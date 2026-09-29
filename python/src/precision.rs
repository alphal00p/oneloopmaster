//! Lossless Python Decimal boundaries for the native arbitrary-precision path.
use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
    types::{PyAnyMethods, PyBool, PyComplex, PyComplexMethods, PyFloat, PyInt},
};
use symbolica::domains::float::{Complex, Float, SingleFloat};

/// Python-facing precision counts decimal significant digits, never binary bits.
pub(crate) fn parse_precision(value: &Bound<'_, PyAny>) -> PyResult<u32> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyTypeError::new_err(
            "prec must be a positive integer number of decimal digits",
        ));
    }
    let digits = value.extract::<u32>().map_err(|_| {
        PyValueError::new_err("prec must be a positive decimal-digit count fitting u32")
    })?;
    if digits == 0 {
        return Err(PyValueError::new_err("prec must be positive"));
    }
    Ok(digits)
}

pub(crate) fn requested_precision(value: Option<&Bound<'_, PyAny>>, default: u32) -> PyResult<u32> {
    value.map(parse_precision).unwrap_or(Ok(default))
}

fn decimal_class(py: Python<'_>) -> PyResult<Bound<'_, PyAny>> {
    py.import("decimal")?.getattr("Decimal")
}

fn decimal_component(value: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    if value.is_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err(
            "boolean values are not numeric inputs",
        ));
    }
    let decimal = decimal_class(value.py())?.call1((value,))?;
    if !decimal.call_method0("is_finite")?.extract::<bool>()? {
        return Err(PyValueError::new_err("Decimal components must be finite"));
    }
    Ok(decimal.unbind())
}

/// A complex number with Python Decimal components.
///
/// Use decimal strings or Decimal inputs to preserve digits beyond binary64.
/// The constructor retains the supplied components without rounding to the
/// ambient Decimal context. ``complex(z)`` explicitly converts to binary64 and
/// can lose precision. Nonfinite components and booleans are rejected.
///
/// Examples
/// --------
/// >>> from decimal import Decimal
/// >>> from symbolica.community.hep import oneloop
/// >>> z = oneloop.DecimalComplex("1.25", "-0.5")
/// >>> assert z.real == Decimal("1.25")
/// >>> assert z.imag == Decimal("-0.5")
#[cfg_attr(
    feature = "community",
    pyclass(frozen, module = "symbolica.community.hep.oneloop")
)]
#[cfg_attr(not(feature = "community"), pyclass(frozen, module = "oneloop_native"))]
pub(crate) struct DecimalComplex {
    real: Py<PyAny>,
    imag: Py<PyAny>,
}

#[pymethods]
impl DecimalComplex {
    /// Construct finite real and imaginary Decimal components; omitted imaginary part is zero.
    ///
    /// Examples
    /// --------
    /// >>> from decimal import Decimal
    /// >>> from symbolica.community.hep import oneloop
    /// >>> z = oneloop.DecimalComplex("1.25", "-0.5")
    /// >>> assert complex(z) == 1.25 - 0.5j
    ///
    /// Parameters
    /// ----------
    /// real : int, float, str or Decimal
    ///     Real component. A float preserves its binary approximation; use a
    ///     string for an exact decimal input.
    /// imag : int, float, str, Decimal or None, optional
    ///     Imaginary component; None means zero.
    #[new]
    #[pyo3(signature = (real, imag=None))]
    fn new(real: &Bound<'_, PyAny>, imag: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        oneloop::record_usage();
        let py = real.py();
        Ok(Self {
            real: decimal_component(real)?,
            imag: if let Some(imag) = imag {
                decimal_component(imag)?
            } else {
                decimal_class(py)?.call1(("0",))?.unbind()
            },
        })
    }

    /// Real component as a Decimal, without conversion to a Python float.
    ///
    /// Examples
    /// --------
    /// >>> from decimal import Decimal
    /// >>> from symbolica.community.hep import oneloop
    /// >>> z = oneloop.DecimalComplex("1.25", "-0.5")
    /// >>> assert z.real == Decimal("1.25")
    #[getter]
    fn real(&self, py: Python<'_>) -> Py<PyAny> {
        self.real.clone_ref(py)
    }

    /// Imaginary component as a Decimal, without conversion to a Python float.
    ///
    /// Examples
    /// --------
    /// >>> from decimal import Decimal
    /// >>> from symbolica.community.hep import oneloop
    /// >>> z = oneloop.DecimalComplex("1.25", "-0.5")
    /// >>> assert z.imag == Decimal("-0.5")
    #[getter]
    fn imag(&self, py: Python<'_>) -> Py<PyAny> {
        self.imag.clone_ref(py)
    }

    /// Convert both components to a built-in complex number; precision beyond binary64 is lost.
    ///
    /// Examples
    /// --------
    /// >>> from decimal import Decimal
    /// >>> from symbolica.community.hep import oneloop
    /// >>> z = oneloop.DecimalComplex("1.25", "-0.5")
    /// >>> assert complex(z) == complex(1.25, -0.5)
    fn __complex__(&self, py: Python<'_>) -> PyResult<Py<PyComplex>> {
        Ok(PyComplex::from_doubles(
            py,
            self.real.bind(py).extract::<f64>()?,
            self.imag.bind(py).extract::<f64>()?,
        )
        .unbind())
    }

    /// Show both Decimal components for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from decimal import Decimal
    /// >>> from symbolica.community.hep import oneloop
    /// >>> z = oneloop.DecimalComplex("1.25", "-0.5")
    /// >>> text = repr(z)
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "DecimalComplex({}, {})",
            self.real.bind(py).repr()?,
            self.imag.bind(py).repr()?
        ))
    }
}

pub(crate) fn needs_arbitrary(value: &Bound<'_, PyAny>) -> PyResult<bool> {
    if value.is_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err(
            "boolean values are not numeric inputs",
        ));
    }
    if value.is_instance_of::<PyFloat>() || value.is_instance_of::<PyComplex>() {
        return Ok(false);
    }
    if value.is_instance_of::<PyInt>() {
        return Ok(value.extract::<i64>().map_or(true, |integer| {
            !(-(1i64 << 53)..=(1i64 << 53)).contains(&integer)
        }));
    }
    if value.is_instance_of::<DecimalComplex>() || value.is_instance(&decimal_class(value.py())?)? {
        return Ok(true);
    }
    #[cfg(feature = "community")]
    if value
        .extract::<symbolica::api::python::PythonExpression>()
        .is_ok()
    {
        // Never lower exact rational or high-precision Symbolica coefficients
        // through binary64 before the requested-precision evaluation.
        return Ok(true);
    }
    Err(PyTypeError::new_err(
        "numeric inputs must be real/complex numbers, Decimal or DecimalComplex",
    ))
}

pub(crate) fn arbitrary_number(value: &Bound<'_, PyAny>, bits: u32) -> PyResult<Complex<Float>> {
    #[cfg(feature = "community")]
    if let Ok(value) = value.extract::<symbolica::api::python::PythonExpression>() {
        use symbolica::atom::{Atom, AtomCore};
        return value
            .expr
            .evaluate_with_prec::<Atom, Complex<Float>>(&std::collections::HashMap::new(), bits)
            .map_err(|error| {
                PyValueError::new_err(format!(
                    "numeric Symbolica input must evaluate without free variables: {error}"
                ))
            });
    }
    if let Ok(value) = value.cast::<PyComplex>() {
        return Ok(Complex::new(
            Float::with_val(bits, value.real()),
            Float::with_val(bits, value.imag()),
        ));
    }
    if value.is_instance_of::<PyFloat>() {
        return Ok(Complex::new(
            Float::with_val(bits, value.extract::<f64>()?),
            Float::new(bits),
        ));
    }
    if let Ok(value) = value.extract::<PyRef<'_, DecimalComplex>>() {
        return Ok(Complex::new(
            decimal_float(value.real.bind(value.py()), bits)?,
            decimal_float(value.imag.bind(value.py()), bits)?,
        ));
    }
    // needs_arbitrary has already rejected booleans and unsupported types.
    // Decimal and arbitrary-sized integers are parsed from their exact text.
    Ok(Complex::new(decimal_float(value, bits)?, Float::new(bits)))
}

fn decimal_float(value: &Bound<'_, PyAny>, bits: u32) -> PyResult<Float> {
    // Decimal(int) is exact and avoids Python's security limit on int->str
    // conversion for integers with more than 4,300 decimal digits.
    let decimal;
    let printable = if value.is_instance_of::<PyInt>() {
        decimal = decimal_class(value.py())?.call1((value,))?;
        &decimal
    } else {
        value
    };
    let text = printable.str()?.extract::<String>()?;
    let result = Float::parse(&text, Some(bits)).map_err(PyValueError::new_err)?;
    if !result.is_finite() {
        return Err(PyValueError::new_err("all input components must be finite"));
    }
    if result.is_zero() && !printable.call_method0("is_zero")?.extract::<bool>()? {
        return Err(PyValueError::new_err(
            "nonzero Decimal input is outside the native float exponent range",
        ));
    }
    Ok(result)
}

pub(crate) fn decimal_output(
    py: Python<'_>,
    value: &Complex<Float>,
    requested: u32,
) -> PyResult<Py<PyAny>> {
    let decimal = decimal_class(py)?;
    let component = |value: &Float| -> PyResult<Py<PyAny>> {
        // Do not pad Numerica's returned precision with fictitious digits.
        let available = (f64::from(value.prec()) * std::f64::consts::LOG10_2).floor() as u32;
        let digits = requested.min(available.max(1)) as usize;
        let text = if value.as_raw().is_nan() {
            "NaN".to_owned()
        } else if !value.is_finite() {
            if value.is_negative() {
                "-Infinity"
            } else {
                "Infinity"
            }
            .to_owned()
        } else {
            value.as_raw().to_string_radix(10, Some(digits))
        };
        // Decimal(str) is exact and ignores decimal.getcontext().prec.
        Ok(decimal.call1((text,))?.unbind())
    };
    Ok(Py::new(
        py,
        DecimalComplex {
            real: component(&value.re)?,
            imag: component(&value.im)?,
        },
    )?
    .into_any())
}
