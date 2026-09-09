//! Thin Python adapters for the core native scalar evaluator.
//!
//! Standalone wheels expose numbers only. Expression interop is enabled only
//! through the `community` feature, linked into one Symbolica community kernel.
//! Identical revisions in separate extension binaries do not imply shared state.
//! Import initializes the core symbols and scalar caches eagerly on the calling
//! thread. Community hosts defer this work to their module initialize hook.
//! Construction/evaluation remains on the calling thread; no hidden worker or
//! automatic precision escalation is introduced here.

mod precision;

use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarEvaluator, ScalarIntegral,
};
use pyo3::{
    exceptions::{PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAnyMethods, PyComplex, PyComplexMethods, PyModuleMethods},
};
use symbolica::{
    atom::Atom,
    domains::float::{Complex, Float, SingleFloat},
    evaluate::ExpressionEvaluator,
};

type Coefficients = (Py<PyAny>, Py<PyAny>, Py<PyAny>);

#[derive(Clone, Copy)]
enum BackendChoice {
    Auto,
    Native,
    SymJit,
    Expression,
}

impl BackendChoice {
    fn parse(name: &str) -> PyResult<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "native" => Ok(Self::Native),
            "symjit" => Ok(Self::SymJit),
            "expression" | "symbolica" => Ok(Self::Expression),
            _ => Err(PyValueError::new_err(
                "backend must be auto, native, symjit or expression",
            )),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Native => "native",
            Self::SymJit => "symjit",
            Self::Expression => "expression",
        }
    }

    fn resolve(self, arbitrary: bool) -> PyResult<EvaluationBackend> {
        let backend = match self {
            Self::Auto if arbitrary => match oneloop::DEFAULT_BACKEND {
                EvaluationBackend::Native => EvaluationBackend::Native,
                _ => EvaluationBackend::Expression,
            },
            Self::Auto => oneloop::DEFAULT_BACKEND,
            Self::Native => EvaluationBackend::Native,
            Self::SymJit => EvaluationBackend::SymJit,
            Self::Expression => EvaluationBackend::Expression,
        };
        if arbitrary && backend == EvaluationBackend::SymJit {
            return Err(PyValueError::new_err(
                "backend='symjit' supports binary64 only; Decimal/large-integer inputs or prec != 16 require native, expression or auto",
            ));
        }
        Ok(backend)
    }
}

fn backend_name(backend: EvaluationBackend) -> &'static str {
    match backend {
        EvaluationBackend::Native => "native",
        EvaluationBackend::SymJit => "symjit",
        EvaluationBackend::Expression => "expression",
    }
}

enum MachineEvaluator {
    Native(NativeEvaluator<f64>),
    SymJit(ScalarEvaluator),
    Expression(ExpressionEvaluator<Complex<f64>>),
}

impl MachineEvaluator {
    fn new(family: ScalarIntegral, backend: EvaluationBackend, rebuild: bool) -> PyResult<Self> {
        match backend {
            EvaluationBackend::Native => NativeEvaluator::new(family)
                .map(Self::Native)
                .map_err(PyRuntimeError::new_err),
            EvaluationBackend::SymJit => build(family, rebuild).map(Self::SymJit),
            EvaluationBackend::Expression => {
                let parameters = (0..family.arity())
                    .map(|index| {
                        symbolica::symbol!(format!(
                            "__oneloop_python_expression_{}_{}",
                            family.name(),
                            index
                        ))
                        .to_atom()
                    })
                    .collect::<Vec<_>>();
                let master = match family {
                    ScalarIntegral::A0 => oneloop::A0(),
                    ScalarIntegral::B0 => oneloop::B0(),
                    ScalarIntegral::DB0 => oneloop::dB0(),
                    ScalarIntegral::C0 => oneloop::C0(),
                    ScalarIntegral::D0 => oneloop::D0(),
                };
                let calls = [0, -1, -2].map(|tag| {
                    let mut arguments = vec![Atom::num(tag)];
                    arguments.extend(parameters.iter().cloned());
                    master.call(&arguments)
                });
                let exact = oneloop::OneLoopExpressions::new()
                    .evaluator(&calls, &parameters)
                    .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
                Ok(Self::Expression(exact.map_coeff(&|value| {
                    Complex::new(value.re.to_f64(), value.im.to_f64())
                })))
            }
        }
    }

    fn evaluate(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>]) -> PyResult<()> {
        match self {
            Self::Native(evaluator) => evaluator
                .evaluate(input, output)
                .map_err(PyRuntimeError::new_err),
            Self::SymJit(evaluator) => evaluator
                .evaluate(input, output)
                .map_err(PyRuntimeError::new_err),
            Self::Expression(evaluator) => {
                evaluator.evaluate(input, output);
                Ok(())
            }
        }
    }

    fn evaluate_batch(
        &mut self,
        input: &[Complex<f64>],
        output: &mut [Complex<f64>],
        rows: usize,
        arity: usize,
    ) -> PyResult<()> {
        match self {
            Self::Native(evaluator) => evaluator
                .evaluate_batch(input, output, rows)
                .map_err(PyRuntimeError::new_err),
            Self::SymJit(evaluator) => evaluator
                .evaluate_batch(input, output, rows)
                .map_err(PyRuntimeError::new_err),
            Self::Expression(evaluator) => {
                for (point, result) in input.chunks_exact(arity).zip(output.chunks_exact_mut(3)) {
                    evaluator.evaluate(point, result);
                }
                Ok(())
            }
        }
    }
}

fn family(name: &str) -> PyResult<ScalarIntegral> {
    match name.trim().to_ascii_uppercase().as_str() {
        "A0" => Ok(ScalarIntegral::A0),
        "B0" => Ok(ScalarIntegral::B0),
        "DB0" => Ok(ScalarIntegral::DB0),
        "C0" => Ok(ScalarIntegral::C0),
        "D0" => Ok(ScalarIntegral::D0),
        _ => Err(PyValueError::new_err(
            "family must be A0, B0, dB0, C0 or D0",
        )),
    }
}

/// A symbolic family selector belongs to the same kernel as this module.
#[cfg(feature = "community")]
fn symbolic_family(value: &Bound<'_, PyAny>) -> PyResult<ScalarIntegral> {
    // Symbolica's extractor accepts both Symbol and a variable Expression,
    // without accepting strings or casting a composite into a made-up name.
    let symbol = value.extract::<symbolica::atom::Symbol>()?;
    let name = symbol.get_name();
    let short = name
        .strip_prefix("oneloopmaster::")
        .or_else(|| name.strip_prefix("oneloop::"))
        .ok_or_else(|| {
            PyValueError::new_err("expected an oneloopmaster/oneloop scalar master symbol")
        })?;
    match short {
        "A0" | "B0" | "dB0" | "C0" | "D0" => family(short),
        _ => Err(PyValueError::new_err("expected A0, B0, dB0, C0 or D0")),
    }
}

fn number(value: &Bound<'_, PyAny>) -> PyResult<Complex<f64>> {
    if let Ok(value) = value.cast::<PyComplex>() {
        Ok(Complex::new(value.real(), value.imag()))
    } else if let Ok(value) = value.extract::<f64>() {
        Ok(Complex::new(value, 0.0))
    } else {
        Err(PyTypeError::new_err(
            "numeric inputs must be Python real or complex numbers",
        ))
    }
}

fn validate(family: ScalarIntegral, input: &[Complex<f64>]) -> PyResult<()> {
    if input.len() != family.arity() {
        return Err(PyValueError::new_err(format!(
            "{} expects {} ordered arguments, including mu_squared last; received {}",
            family.name(),
            family.arity(),
            input.len()
        )));
    }
    if input.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(PyValueError::new_err("all input components must be finite"));
    }
    let momenta = match family {
        ScalarIntegral::A0 => 0,
        ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
        ScalarIntegral::C0 => 3,
        ScalarIntegral::D0 => 6,
    };
    if input[..momenta].iter().any(|z| z.im != 0.0) {
        return Err(PyValueError::new_err(
            "external squared invariants must be real",
        ));
    }
    if input[momenta..input.len() - 1].iter().any(|z| z.im > 0.0) {
        return Err(PyValueError::new_err(
            "squared masses require nonpositive imaginary parts",
        ));
    }
    let mu = input[input.len() - 1];
    if mu.im != 0.0 || mu.re <= 0.0 {
        return Err(PyValueError::new_err(
            "mu_squared must be positive and real",
        ));
    }
    Ok(())
}

fn coefficients(py: Python<'_>, output: &[Complex<f64>]) -> Coefficients {
    let convert = |z: &Complex<f64>| PyComplex::from_doubles(py, z.re, z.im).unbind().into_any();
    (
        convert(&output[0]),
        convert(&output[1]),
        convert(&output[2]),
    )
}

fn arbitrary_required(arguments: &[Bound<'_, PyAny>], digits: u32) -> PyResult<bool> {
    let mut required = digits != 16;
    for value in arguments {
        // Validate every input type, even after a Decimal has selected this path.
        required |= precision::needs_arbitrary(value)?;
    }
    Ok(required)
}

fn validate_arbitrary(family: ScalarIntegral, input: &[Complex<Float>]) -> PyResult<()> {
    if input.len() != family.arity() {
        return Err(PyValueError::new_err(format!(
            "{} expects {} ordered arguments, including mu_squared last; received {}",
            family.name(),
            family.arity(),
            input.len()
        )));
    }
    if input.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(PyValueError::new_err("all input components must be finite"));
    }
    let momenta = match family {
        ScalarIntegral::A0 => 0,
        ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
        ScalarIntegral::C0 => 3,
        ScalarIntegral::D0 => 6,
    };
    if input[..momenta].iter().any(|z| !z.im.is_zero()) {
        return Err(PyValueError::new_err(
            "external squared invariants must be real",
        ));
    }
    if input[momenta..input.len() - 1]
        .iter()
        .any(|z| !z.im.is_zero() && !z.im.is_negative())
    {
        return Err(PyValueError::new_err(
            "squared masses require nonpositive imaginary parts",
        ));
    }
    let mu = &input[input.len() - 1];
    if !mu.im.is_zero() || mu.re.is_zero() || mu.re.is_negative() {
        return Err(PyValueError::new_err(
            "mu_squared must be positive and real",
        ));
    }
    Ok(())
}

fn arbitrary_coefficients(
    py: Python<'_>,
    output: &[Complex<Float>],
    digits: u32,
) -> PyResult<Coefficients> {
    Ok((
        precision::decimal_output(py, &output[0], digits)?,
        precision::decimal_output(py, &output[1], digits)?,
        precision::decimal_output(py, &output[2], digits)?,
    ))
}

fn arbitrary_single(
    py: Python<'_>,
    evaluator: &mut PrecisionEvaluator,
    arguments: &[Bound<'_, PyAny>],
    digits: u32,
) -> PyResult<Coefficients> {
    let bits = evaluator.binary_precision();
    let input = arguments
        .iter()
        .map(|value| precision::arbitrary_number(value, bits))
        .collect::<PyResult<Vec<_>>>()?;
    validate_arbitrary(evaluator.family(), &input)?;
    let mut output =
        core::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(bits), Float::new(bits)));
    evaluator
        .evaluate(&input, &mut output)
        .map_err(PyRuntimeError::new_err)?;
    arbitrary_coefficients(py, &output, digits)
}

fn single_inputs(
    py: Python<'_>,
    family: ScalarIntegral,
    arguments: &[Bound<'_, PyAny>],
    digits: u32,
    rebuild: bool,
    backend: BackendChoice,
) -> PyResult<Coefficients> {
    if arbitrary_required(arguments, digits)? {
        let mut evaluator =
            PrecisionEvaluator::with_backend(family, digits, backend.resolve(true)?)
                .map_err(PyValueError::new_err)?;
        arbitrary_single(py, &mut evaluator, arguments, digits)
    } else {
        let input = arguments.iter().map(number).collect::<PyResult<Vec<_>>>()?;
        single(py, family, &input, rebuild, backend.resolve(false)?)
    }
}

fn build(family: ScalarIntegral, rebuild: bool) -> PyResult<ScalarEvaluator> {
    if rebuild {
        ScalarEvaluator::rebuild(family)
    } else {
        ScalarEvaluator::cached(family)
    }
    .map_err(PyRuntimeError::new_err)
}

fn single(
    py: Python<'_>,
    family: ScalarIntegral,
    input: &[Complex<f64>],
    rebuild: bool,
    backend: EvaluationBackend,
) -> PyResult<Coefficients> {
    validate(family, input)?;
    let mut output = [Complex::new(0.0, 0.0); 3];
    if rebuild {
        MachineEvaluator::new(family, backend, true)?.evaluate(input, &mut output)?;
    } else {
        oneloop::evaluate_with_backend(family, input, &mut output, backend)
            .map_err(PyRuntimeError::new_err)?;
    }
    Ok(coefficients(py, &output))
}

/// Return whether eager core startup finished successfully; never initialize it.
#[pyfunction]
fn is_initialized() -> bool {
    oneloop::is_initialized()
}

/// Reusable numeric evaluator. prec counts decimal digits; Decimal inputs select
/// arbitrary precision even at the default prec=16. Ordinary default inputs use
/// the default binary64 backend unless explicitly selected. Arguments follow the Rust scalar API,
/// with mu_squared last. Outputs are (finite, simple_pole, double_pole).
/// Use the creating thread; respect Symbolica's license and stack requirements.
#[cfg_attr(
    feature = "community",
    pyclass(name = "Evaluator", unsendable, module = "symbolica.community.oneloop")
)]
#[cfg_attr(
    not(feature = "community"),
    pyclass(name = "Evaluator", unsendable, module = "oneloop_native")
)]
struct Evaluator {
    family: ScalarIntegral,
    evaluator: Option<(EvaluationBackend, MachineEvaluator)>,
    prec: u32,
    backend: BackendChoice,
    arbitrary: Option<(u32, PrecisionEvaluator)>,
}

impl Evaluator {
    fn machine(&mut self, backend: EvaluationBackend) -> PyResult<&mut MachineEvaluator> {
        if self
            .evaluator
            .as_ref()
            .is_none_or(|(current, _)| *current != backend)
        {
            self.evaluator = Some((backend, MachineEvaluator::new(self.family, backend, false)?));
        }
        Ok(&mut self.evaluator.as_mut().unwrap().1)
    }

    fn arbitrary(
        &mut self,
        digits: u32,
        backend: EvaluationBackend,
    ) -> PyResult<&mut PrecisionEvaluator> {
        if self
            .arbitrary
            .as_ref()
            .is_none_or(|(current, value)| *current != digits || value.backend() != backend)
        {
            let evaluator = PrecisionEvaluator::with_backend(self.family, digits, backend)
                .map_err(PyValueError::new_err)?;
            self.arbitrary = Some((digits, evaluator));
        }
        Ok(&mut self.arbitrary.as_mut().unwrap().1)
    }
}

#[pymethods]
impl Evaluator {
    #[new]
    #[pyo3(signature = (family, rebuild=false, *, prec=16, backend="auto"))]
    fn new(
        family: &Bound<'_, PyAny>,
        rebuild: bool,
        #[pyo3(from_py_with = precision::parse_precision)] prec: u32,
        backend: &str,
    ) -> PyResult<Self> {
        #[cfg(feature = "community")]
        let family = symbolic_family(family)?;
        #[cfg(not(feature = "community"))]
        let family = self::family(&family.extract::<String>()?)?;
        let backend = BackendChoice::parse(backend)?;
        let resolved = backend.resolve(prec != 16)?;
        Ok(Self {
            family,
            evaluator: if prec == 16 {
                Some((resolved, MachineEvaluator::new(family, resolved, rebuild)?))
            } else {
                None
            },
            prec,
            backend,
            arbitrary: if prec != 16 {
                Some((
                    prec,
                    PrecisionEvaluator::with_backend(family, prec, resolved)
                        .map_err(PyValueError::new_err)?,
                ))
            } else {
                None
            },
        })
    }

    #[getter]
    fn family(&self) -> &'static str {
        self.family.name()
    }

    #[getter]
    fn arity(&self) -> usize {
        self.family.arity()
    }

    #[getter]
    fn prec(&self) -> u32 {
        self.prec
    }

    /// Constructor's requested backend. Per-call overrides do not change it.
    #[getter]
    fn backend(&self) -> &'static str {
        self.backend.name()
    }

    #[pyo3(signature = (arguments, *, prec=None, backend=None))]
    fn evaluate(
        &mut self,
        py: Python<'_>,
        arguments: Vec<Bound<'_, PyAny>>,
        prec: Option<&Bound<'_, PyAny>>,
        backend: Option<&str>,
    ) -> PyResult<Coefficients> {
        let digits = precision::requested_precision(prec, self.prec)?;
        let backend = backend
            .map(BackendChoice::parse)
            .transpose()?
            .unwrap_or(self.backend);
        if arbitrary_required(&arguments, digits)? {
            return arbitrary_single(
                py,
                self.arbitrary(digits, backend.resolve(true)?)?,
                &arguments,
                digits,
            );
        }
        let input = arguments.iter().map(number).collect::<PyResult<Vec<_>>>()?;
        validate(self.family, &input)?;
        let mut output = [Complex::new(0.0, 0.0); 3];
        self.machine(backend.resolve(false)?)?
            .evaluate(&input, &mut output)?;
        Ok(coefficients(py, &output))
    }

    /// Evaluate rows of ordered arguments; return one coefficient tuple per row.
    /// Empty input returns an empty list. No NumPy dependency is required.
    #[pyo3(signature = (rows, *, prec=None, backend=None))]
    fn evaluate_batch(
        &mut self,
        py: Python<'_>,
        rows: Vec<Vec<Bound<'_, PyAny>>>,
        prec: Option<&Bound<'_, PyAny>>,
        backend: Option<&str>,
    ) -> PyResult<Vec<Coefficients>> {
        let digits = precision::requested_precision(prec, self.prec)?;
        let backend = backend
            .map(BackendChoice::parse)
            .transpose()?
            .unwrap_or(self.backend);
        let count = rows.len();
        if count == 0 {
            backend.resolve(digits != 16)?;
            return Ok(Vec::new());
        }
        let input_count = count
            .checked_mul(self.family.arity())
            .ok_or_else(|| PyValueError::new_err("input batch size overflow"))?;
        let output_count = count
            .checked_mul(3)
            .ok_or_else(|| PyValueError::new_err("output batch size overflow"))?;
        let mut arbitrary = digits != 16;
        for row in &rows {
            arbitrary |= arbitrary_required(row, digits)?;
        }
        if arbitrary {
            let family = self.family;
            let evaluator = self.arbitrary(digits, backend.resolve(true)?)?;
            let bits = evaluator.binary_precision();
            let mut input = Vec::with_capacity(input_count);
            for row in &rows {
                let start = input.len();
                for value in row {
                    input.push(precision::arbitrary_number(value, bits)?);
                }
                validate_arbitrary(family, &input[start..])?;
            }
            let mut output = (0..output_count)
                .map(|_| Complex::new(Float::new(bits), Float::new(bits)))
                .collect::<Vec<_>>();
            evaluator
                .evaluate_batch(&input, &mut output, count)
                .map_err(PyRuntimeError::new_err)?;
            return output
                .chunks_exact(3)
                .map(|row| arbitrary_coefficients(py, row, digits))
                .collect();
        }
        let mut input = Vec::with_capacity(input_count);
        for row in rows {
            let start = input.len();
            for value in row {
                input.push(number(&value)?);
            }
            validate(self.family, &input[start..])?;
        }
        let mut output = vec![Complex::new(0.0, 0.0); output_count];
        let arity = self.family.arity();
        self.machine(backend.resolve(false)?)?
            .evaluate_batch(&input, &mut output, count, arity)?;
        Ok(output
            .chunks_exact(3)
            .map(|row| coefficients(py, row))
            .collect())
    }

    /// Rebuild current workspaces transactionally. Native rebuild prepares fresh
    /// constants/workspace only; SymJIT rebuild recompiles its evaluator.
    fn rebuild(&mut self) -> PyResult<()> {
        let evaluator = self
            .evaluator
            .as_ref()
            .map(|(backend, _)| {
                MachineEvaluator::new(self.family, *backend, true).map(|value| (*backend, value))
            })
            .transpose()?;
        let arbitrary = self
            .arbitrary
            .as_ref()
            .map(|(digits, value)| {
                PrecisionEvaluator::with_backend(self.family, *digits, value.backend())
                    .map(|value| (*digits, value))
            })
            .transpose()
            .map_err(PyRuntimeError::new_err)?;
        self.evaluator = evaluator;
        self.arbitrary = arbitrary;
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!(
            "Evaluator('{}', arity={}, prec={}, backend='{}', coefficient_order=(0,-1,-2))",
            self.family.name(),
            self.family.arity(),
            self.prec,
            self.backend.name(),
        )
    }
}

macro_rules! scalar {
    ($rust:ident, $python:literal, $kind:ident, [$($argument:ident),+]) => {
        #[pyfunction(name=$python, signature=($($argument,)+ mu_squared=None, *, rebuild=false, prec=16, backend="auto"))]
        #[doc = "Return (finite, simple_pole, double_pole). prec counts decimal digits; omitted mu_squared is exactly 1. Decimal inputs/results preserve arbitrary precision. backend is auto, native, symjit (binary64 only), or expression."]
        // Preserve the conventional ordered scalar-integral Python signature.
        #[allow(clippy::too_many_arguments)]
        fn $rust(py: Python<'_>, $($argument: &Bound<'_, PyAny>,)+ mu_squared: Option<&Bound<'_, PyAny>>, rebuild: bool, #[pyo3(from_py_with = precision::parse_precision)] prec: u32, backend: &str) -> PyResult<Coefficients> {
            let default_scale = 1u8.into_pyobject(py)?.into_any();
            let mut arguments = vec![$($argument.clone()),+];
            arguments.push(mu_squared.cloned().unwrap_or(default_scale));
            single_inputs(py, ScalarIntegral::$kind, &arguments, prec, rebuild, BackendChoice::parse(backend)?)
        }
    };
}

scalar!(a0, "A0", A0, [mass_squared]);
scalar!(
    b0,
    "B0",
    B0,
    [momentum_squared, mass_0_squared, mass_1_squared]
);
scalar!(
    db0,
    "dB0",
    DB0,
    [momentum_squared, mass_0_squared, mass_1_squared]
);
scalar!(
    c0,
    "C0",
    C0,
    [
        p1_squared,
        p2_squared,
        p3_squared,
        mass_0_squared,
        mass_1_squared,
        mass_2_squared
    ]
);
scalar!(
    d0,
    "D0",
    D0,
    [
        p1_squared,
        p2_squared,
        p3_squared,
        p4_squared,
        s12,
        s23,
        mass_0_squared,
        mass_1_squared,
        mass_2_squared,
        mass_3_squared
    ]
);

mod inspection;

fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    inspection::register(module)?;
    // Registering Python types/functions does not initialize Symbolica symbols.
    module.add_class::<Evaluator>()?;
    module.add_class::<precision::DecimalComplex>()?;
    module.add_function(wrap_pyfunction!(a0, module)?)?;
    module.add_function(wrap_pyfunction!(b0, module)?)?;
    module.add_function(wrap_pyfunction!(db0, module)?)?;
    module.add_function(wrap_pyfunction!(c0, module)?)?;
    module.add_function(wrap_pyfunction!(d0, module)?)?;
    module.add_function(wrap_pyfunction!(is_initialized, module)?)?;
    for (lower, upper) in [
        ("a0", "A0"),
        ("b0", "B0"),
        ("db0", "dB0"),
        ("c0", "C0"),
        ("d0", "D0"),
    ] {
        module.add(lower, module.getattr(upper)?)?;
    }
    module.add("COEFFICIENT_ORDER", (0, -1, -2))?;
    module.add("DEFAULT_BACKEND", backend_name(oneloop::DEFAULT_BACKEND))?;
    module.add(
        "SYMBOLICA_REVISION",
        "fb845d34bda8ccf1fedef6544d3aa46dc24944e3",
    )?;
    module.add("EXPRESSION_INTEROP", cfg!(feature = "community"))?;
    #[cfg(feature = "community")]
    community::register(module)?;
    Ok(())
}

#[cfg(all(feature = "extension-module", not(feature = "community")))]
#[pymodule]
fn oneloop_native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    oneloop::initialize().map_err(PyRuntimeError::new_err)?;
    register(module)
}

#[cfg(feature = "community")]
mod community {
    use super::*;
    use symbolica::{
        api::python::{
            ConvertibleToExpression, PythonExpression, PythonExpressionEvaluator,
            SymbolicaCommunityModule,
        },
        atom::Atom,
    };

    pub struct CommunityModule;

    impl SymbolicaCommunityModule for CommunityModule {
        fn get_name() -> String {
            "oneloop".to_owned()
        }
        fn register_module(module: &Bound<'_, PyModule>) -> PyResult<()> {
            super::register(module)
        }
        fn initialize(_py: Python<'_>) -> PyResult<()> {
            oneloop::initialize().map_err(PyRuntimeError::new_err)
        }
    }

    /// Return compact master expressions sharing this community kernel's state.
    /// Keep native definitions by compiling combinations with compile_native.
    #[pyfunction]
    fn master_coefficients(master: PythonExpression) -> PyResult<Vec<PythonExpression>> {
        let (family, input) =
            oneloop::master_arguments(master.expr.as_view()).map_err(PyValueError::new_err)?;
        let master = match family {
            ScalarIntegral::A0 => oneloop::A0(),
            ScalarIntegral::B0 => oneloop::B0(),
            ScalarIntegral::DB0 => oneloop::dB0(),
            ScalarIntegral::C0 => oneloop::C0(),
            ScalarIntegral::D0 => oneloop::D0(),
        };
        Ok([0, -1, -2]
            .into_iter()
            .map(|tag| {
                let arguments = std::iter::once(Atom::num(tag))
                    .chain(input.iter().cloned())
                    .collect::<Vec<_>>();
                master.call(arguments.as_slice()).into()
            })
            .collect())
    }

    /// Compile mixed Symbolica expressions with all transparent OneLOop definitions.
    /// Returns the host Symbolica Evaluator; use its complex-valued evaluation API.
    #[pyfunction]
    fn compile_native(
        expressions: Vec<ConvertibleToExpression>,
        parameters: Vec<PythonExpression>,
    ) -> PyResult<PythonExpressionEvaluator> {
        let expressions = expressions
            .into_iter()
            .map(|e| e.to_expression().expr)
            .collect::<Vec<_>>();
        let parameters = parameters.into_iter().map(|p| p.expr).collect::<Vec<_>>();
        let exact = oneloop::OneLoopExpressions::new()
            .evaluator(&expressions, &parameters)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(PythonExpressionEvaluator {
            rational_constants: exact.get_constants().to_vec(),
            eval_complex: exact.map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64())),
            eval_real: None,
            jit_real: None,
            jit_complex: None,
            eval_double_float: None,
            eval_double_float_complex: None,
            eval_arb_prec: None,
            eval_arb_prec_complex: None,
            jit_compile: false,
            jit_settings: oneloop::jit_settings(),
        })
    }

    pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
        module.add_function(wrap_pyfunction!(master_coefficients, module)?)?;
        module.add_function(wrap_pyfunction!(compile_native, module)?)?;
        Ok(())
    }
}

#[cfg(feature = "community")]
pub use community::CommunityModule;
