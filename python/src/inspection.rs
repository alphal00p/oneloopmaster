//! Expression inspection. Standalone text never crosses Symbolica kernel state.
use pyo3::{exceptions::PyValueError, prelude::*, types::PyTuple};
use symbolica::atom::Atom;

fn options(max_nodes: usize, max_depth: usize) -> oneloop::ExpressionOptions {
    oneloop::ExpressionOptions {
        max_nodes,
        max_depth,
    }
}

fn selected(coefficient: Option<i32>) -> PyResult<Option<usize>> {
    match coefficient {
        None => Ok(None),
        Some(0) => Ok(Some(0)),
        Some(-1) => Ok(Some(1)),
        Some(-2) => Ok(Some(2)),
        Some(_) => Err(PyValueError::new_err("coefficient must be 0, -1 or -2")),
    }
}

#[cfg(not(feature = "community"))]
fn declare(names: &[String], positive: bool) -> PyResult<()> {
    use symbolica::atom::{NamespacedSymbol, Symbol, SymbolAttribute, SymbolBuilder};
    for name in names {
        if !name.split("::").all(|part| {
            let mut chars = part.chars();
            chars.next().is_some_and(|c| c.is_alphabetic() || c == '_')
                && chars.all(|c| c.is_alphanumeric() || c == '_')
        }) {
            return Err(PyValueError::new_err(
                "assumption declarations must be variable names, optionally namespaced",
            ));
        }
        let name = if name.contains("::") {
            name.clone()
        } else {
            format!("oneloop_input::{name}")
        };
        let namespaced = NamespacedSymbol::try_parse(name).expect("namespace was provided");
        if let Some(existing) = Symbol::get_symbol(namespaced.clone()) {
            if if positive {
                existing.is_positive()
            } else {
                existing.is_real()
            } {
                continue;
            }
            return Err(PyValueError::new_err(format!(
                "{} already exists without the requested {} attribute; use a new name instead of retagging it",
                existing.get_name(),
                if positive { "Positive" } else { "Real" }
            )));
        }
        let attributes = if positive {
            vec![SymbolAttribute::Positive]
        } else {
            vec![SymbolAttribute::Real]
        };
        SymbolBuilder::new(namespaced)
            .with_attributes(attributes)
            .build()
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
    }
    Ok(())
}

/// Return the complete native coefficients as parseable Symbolica strings.
///
/// Arguments are expression strings; mu_squared is the last argument. Symbols
/// without explicit namespaces belong to oneloop_input. `real` and `positive`
/// declare new variables before parsing; existing symbols cannot be retagged.
/// Native attribute syntax such as `oneloop_input::{real}::p` is also accepted.
/// The returned text includes namespaces and attributes. It is not an Expression
/// object from a separately imported Symbolica extension.
#[cfg(not(feature = "community"))]
#[pyfunction(signature = (family, arguments, *, coefficient=None, max_nodes=1_000_000, max_depth=512, real=None, positive=None))]
// Keep the public Python keyword controls explicit in its generated signature.
#[allow(clippy::too_many_arguments)]
fn get_expression(
    py: Python<'_>,
    family: &str,
    arguments: Vec<String>,
    coefficient: Option<i32>,
    max_nodes: usize,
    max_depth: usize,
    real: Option<Vec<String>>,
    positive: Option<Vec<String>>,
) -> PyResult<Py<PyAny>> {
    use symbolica::{atom::AtomCore, parser::ParseSettings, printer::PrintOptions};
    let family = super::family(family)?;
    let coefficient = selected(coefficient)?;
    if arguments.len() != family.arity() {
        return Err(PyValueError::new_err(format!(
            "{} expects {} arguments including mu_squared",
            family.name(),
            family.arity()
        )));
    }
    // Positive implies Real: declare it first if a name appears in both lists.
    declare(positive.as_deref().unwrap_or_default(), true)?;
    declare(real.as_deref().unwrap_or_default(), false)?;
    let input = arguments
        .into_iter()
        .map(|a| {
            Atom::parse(a, "oneloop_input", ParseSettings::default()).map_err(PyValueError::new_err)
        })
        .collect::<PyResult<Vec<_>>>()?;
    let series =
        oneloop::get_expression_with_options(family, &input, options(max_nodes, max_depth))
            .map_err(PyValueError::new_err)?;
    let strings = series
        .coefficients()
        .iter()
        .map(|a| a.printer(PrintOptions::full()).to_string())
        .collect::<Vec<_>>();
    if let Some(index) = coefficient {
        Ok(strings[index]
            .clone()
            .into_pyobject(py)?
            .unbind()
            .into_any())
    } else {
        Ok(PyTuple::new(py, strings)?.unbind().into_any())
    }
}

/// Return complete native Expression objects in the shared community kernel.
///
/// Use host S(..., is_real=True/is_positive=True) to declare assumptions before
/// constructing arguments. No external/standalone Expression objects are used.
#[cfg(feature = "community")]
#[pyfunction(signature = (family, arguments, *, coefficient=None, max_nodes=1_000_000, max_depth=512))]
fn get_expression(
    py: Python<'_>,
    family: &str,
    arguments: Vec<symbolica::api::python::ConvertibleToExpression>,
    coefficient: Option<i32>,
    max_nodes: usize,
    max_depth: usize,
) -> PyResult<Py<PyAny>> {
    use symbolica::api::python::PythonExpression;
    let family = super::family(family)?;
    let coefficient = selected(coefficient)?;
    let input = arguments
        .into_iter()
        .map(|a| a.to_expression().expr)
        .collect::<Vec<Atom>>();
    let series =
        oneloop::get_expression_with_options(family, &input, options(max_nodes, max_depth))
            .map_err(PyValueError::new_err)?;
    let expressions = series
        .into_coefficients()
        .into_iter()
        .map(|a| Py::new(py, PythonExpression::from(a)))
        .collect::<PyResult<Vec<_>>>()?;
    if let Some(index) = coefficient {
        Ok(expressions[index].clone_ref(py).into_any())
    } else {
        Ok(PyTuple::new(py, expressions)?.unbind().into_any())
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(get_expression, module)?)?;
    Ok(())
}
