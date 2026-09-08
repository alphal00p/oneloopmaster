//! Thin Python adapters for the core native scalar evaluator.
//!
//! Standalone wheels expose numbers only. Expression interop is enabled only
//! through the `community` feature, linked into one Symbolica community kernel.
//! Identical revisions in separate extension binaries do not imply shared state.
//! Import initializes the core symbols and scalar caches eagerly on the calling
//! thread. Community hosts defer this work to their module initialize hook.
//! Construction/evaluation remains on the calling thread; no hidden worker or
//! automatic precision escalation is introduced here.

use oneloop::{ScalarEvaluator, ScalarIntegral};
use pyo3::{
    exceptions::{PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAnyMethods, PyComplex, PyComplexMethods, PyModuleMethods},
};
use symbolica::domains::float::Complex;

type Coefficients = (Py<PyComplex>, Py<PyComplex>, Py<PyComplex>);

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
    let convert = |z: &Complex<f64>| PyComplex::from_doubles(py, z.re, z.im).unbind();
    (
        convert(&output[0]),
        convert(&output[1]),
        convert(&output[2]),
    )
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
) -> PyResult<Coefficients> {
    validate(family, input)?;
    let mut output = [Complex::new(0.0, 0.0); 3];
    if rebuild {
        build(family, true)?
            .evaluate(input, &mut output)
            .map_err(PyRuntimeError::new_err)?;
    } else {
        oneloop::evaluate(family, input, &mut output).map_err(PyRuntimeError::new_err)?;
    }
    Ok(coefficients(py, &output))
}

/// Return whether eager core startup finished successfully; never initialize it.
#[pyfunction]
fn is_initialized() -> bool {
    oneloop::is_initialized()
}

/// Reusable binary64-complex evaluator. Arguments follow the Rust scalar API,
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
    evaluator: ScalarEvaluator,
}

#[pymethods]
impl Evaluator {
    #[new]
    #[pyo3(signature = (family, rebuild=false))]
    fn new(family: &str, rebuild: bool) -> PyResult<Self> {
        let family = self::family(family)?;
        Ok(Self {
            family,
            evaluator: build(family, rebuild)?,
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

    fn evaluate(
        &mut self,
        py: Python<'_>,
        arguments: Vec<Bound<'_, PyAny>>,
    ) -> PyResult<Coefficients> {
        let input = arguments.iter().map(number).collect::<PyResult<Vec<_>>>()?;
        validate(self.family, &input)?;
        let mut output = [Complex::new(0.0, 0.0); 3];
        self.evaluator
            .evaluate(&input, &mut output)
            .map_err(PyRuntimeError::new_err)?;
        Ok(coefficients(py, &output))
    }

    /// Evaluate rows of ordered arguments; return one coefficient tuple per row.
    /// Empty input returns an empty list. No NumPy dependency is required.
    fn evaluate_batch(
        &mut self,
        py: Python<'_>,
        rows: Vec<Vec<Bound<'_, PyAny>>>,
    ) -> PyResult<Vec<Coefficients>> {
        let count = rows.len();
        if count == 0 {
            return Ok(Vec::new());
        }
        let input_count = count
            .checked_mul(self.family.arity())
            .ok_or_else(|| PyValueError::new_err("input batch size overflow"))?;
        let output_count = count
            .checked_mul(3)
            .ok_or_else(|| PyValueError::new_err("output batch size overflow"))?;
        let mut input = Vec::with_capacity(input_count);
        for row in rows {
            let start = input.len();
            for value in row {
                input.push(number(&value)?);
            }
            validate(self.family, &input[start..])?;
        }
        let mut output = vec![Complex::new(0.0, 0.0); output_count];
        self.evaluator
            .evaluate_batch(&input, &mut output, count)
            .map_err(PyRuntimeError::new_err)?;
        Ok(output
            .chunks_exact(3)
            .map(|row| coefficients(py, row))
            .collect())
    }

    /// Force a fresh native evaluator build; replace this instance only on success.
    fn rebuild(&mut self) -> PyResult<()> {
        self.evaluator = build(self.family, true)?;
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!(
            "Evaluator('{}', arity={}, coefficient_order=(0,-1,-2))",
            self.family.name(),
            self.family.arity()
        )
    }
}

macro_rules! scalar {
    ($rust:ident, $python:literal, $kind:ident, [$($argument:ident),+]) => {
        #[pyfunction(name=$python, signature=($($argument,)+ mu_squared=1.0, *, rebuild=false))]
        #[doc = "Return (finite, simple_pole, double_pole); masses, momenta and mu_squared are squared quantities."]
        // Preserve the conventional ordered scalar-integral Python signature.
        #[allow(clippy::too_many_arguments)]
        fn $rust(py: Python<'_>, $($argument: &Bound<'_, PyAny>,)+ mu_squared: f64, rebuild: bool) -> PyResult<Coefficients> {
            let mut input = vec![$(number($argument)?),+];
            input.push(Complex::new(mu_squared, 0.0));
            single(py, ScalarIntegral::$kind, &input, rebuild)
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

fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    // Registering Python types/functions does not initialize Symbolica symbols.
    module.add_class::<Evaluator>()?;
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
    fn master_coefficients(
        family: &str,
        arguments: Vec<ConvertibleToExpression>,
    ) -> PyResult<Vec<PythonExpression>> {
        let family = super::family(family)?;
        if arguments.len() != family.arity() {
            return Err(PyValueError::new_err(format!(
                "{} expects {} arguments including mu_squared",
                family.name(),
                family.arity()
            )));
        }
        let master = match family {
            ScalarIntegral::A0 => oneloop::A0(),
            ScalarIntegral::B0 => oneloop::B0(),
            ScalarIntegral::DB0 => oneloop::dB0(),
            ScalarIntegral::C0 => oneloop::C0(),
            ScalarIntegral::D0 => oneloop::D0(),
        };
        let input = arguments
            .into_iter()
            .map(|a| a.to_expression().expr)
            .collect::<Vec<_>>();
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
