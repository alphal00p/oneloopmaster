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

/// A complex number whose real and imaginary components are exact Python Decimals.
/// Construction and returned components do not round through the ambient context.
/// Converting with complex(value) is an explicit, lossy binary64 operation.
#[cfg_attr(
    feature = "community",
    pyclass(frozen, module = "symbolica.community.oneloop")
)]
#[cfg_attr(not(feature = "community"), pyclass(frozen, module = "oneloop_native"))]
pub(crate) struct DecimalComplex {
    real: Py<PyAny>,
    imag: Py<PyAny>,
}

#[pymethods]
impl DecimalComplex {
    #[new]
    #[pyo3(signature = (real, imag=None))]
    fn new(real: &Bound<'_, PyAny>, imag: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
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

    #[getter]
    fn real(&self, py: Python<'_>) -> Py<PyAny> {
        self.real.clone_ref(py)
    }

    #[getter]
    fn imag(&self, py: Python<'_>) -> Py<PyAny> {
        self.imag.clone_ref(py)
    }

    fn __complex__(&self, py: Python<'_>) -> PyResult<Py<PyComplex>> {
        Ok(PyComplex::from_doubles(
            py,
            self.real.bind(py).extract::<f64>()?,
            self.imag.bind(py).extract::<f64>()?,
        )
        .unbind())
    }

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
    Err(PyTypeError::new_err(
        "numeric inputs must be real/complex numbers, Decimal or DecimalComplex",
    ))
}

pub(crate) fn arbitrary_number(value: &Bound<'_, PyAny>, bits: u32) -> PyResult<Complex<Float>> {
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
