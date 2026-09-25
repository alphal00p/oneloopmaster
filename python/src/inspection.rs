//! Expression inspection in the host's single Symbolica kernel.
//!
//! A separate numeric-only extension must never borrow foreign Expression
//! objects. These functions are therefore exposed only by the community host.
#[cfg(feature = "community")]
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyList, PyTuple},
};
#[cfg(feature = "community")]
use symbolica::{
    api::python::{PythonExpression, PythonReplacement, PythonTransformer},
    atom::AliasedAtom,
    id::Replacement,
    transformer::Transformer,
};

/// A complete Symbolica expression with transparent common-subexpression bindings.
/// No numerical callback or OneLOop helper is needed to interpret its definitions.
#[cfg(feature = "community")]
#[pyclass(name = "SharedExpression", module = "symbolica.community.oneloop")]
struct PythonSharedExpression {
    expression: AliasedAtom,
}

#[cfg(feature = "community")]
#[pymethods]
impl PythonSharedExpression {
    #[getter]
    fn root(&self) -> PythonExpression {
        self.expression.get_root().clone().into()
    }

    /// Every pair is (alias Expression, complete defining Expression).
    #[getter]
    fn definitions(&self) -> Vec<(PythonExpression, PythonExpression)> {
        let mut definitions = self.expression.get_aliases().iter().collect::<Vec<_>>();
        definitions.sort_by(|a, b| a.0.cmp(b.0));
        definitions
            .into_iter()
            .map(|(a, b)| (a.clone().into(), b.clone().into()))
            .collect()
    }

    #[getter]
    fn byte_size(&self) -> usize {
        self.expression.get_byte_size()
    }

    #[getter]
    fn num_definitions(&self) -> usize {
        self.expression.get_aliases().len()
    }

    fn __repr__(&self) -> String {
        format!(
            "SharedExpression(root={}, definitions={}, bytes={})",
            self.expression.get_root(),
            self.expression.get_aliases().len(),
            self.expression.get_byte_size()
        )
    }

    /// Substitute every binding; potentially much larger than the shared form.
    #[pyo3(signature = (*, max_nodes=100_000_000, max_depth=4096))]
    fn to_expression(&self, max_nodes: usize, max_depth: usize) -> PyResult<PythonExpression> {
        oneloop::select_branch_shared(
            &self.expression,
            &[],
            oneloop::ExpressionOptions {
                max_nodes,
                max_depth,
            },
        )
        .map(PythonExpression::from)
        .map_err(PyValueError::new_err)
    }
}

#[cfg(feature = "community")]
fn selected(coefficient: Option<i32>) -> PyResult<Option<usize>> {
    match coefficient {
        None => Ok(None),
        Some(0) => Ok(Some(0)),
        Some(-1) => Ok(Some(1)),
        Some(-2) => Ok(Some(2)),
        Some(_) => Err(PyValueError::new_err("coefficient must be 0, -1 or -2")),
    }
}

/// Expand a complete master call to its three Laurent coefficients.
///
/// The input is a Symbolica Expression such as B0(psq,m2,m2,mu_squared),
/// with physical arguments only. Attributes belong to the input's Symbolica
/// symbols; this function never parses strings or declares assumptions.
/// Optional branch_rules select conditions during expansion, before discarded
/// branches are built. Probe values are never substituted into returned bodies.
#[cfg(feature = "community")]
#[pyfunction(signature = (master, *, coefficient=None, max_nodes=10_000_000, max_depth=512, branch_rules=None, shared=false))]
fn get_expression(
    py: Python<'_>,
    master: PythonExpression,
    coefficient: Option<i32>,
    max_nodes: usize,
    max_depth: usize,
    branch_rules: Option<Vec<PythonReplacement>>,
    shared: bool,
) -> PyResult<Py<PyAny>> {
    let coefficient = selected(coefficient)?;
    let options = oneloop::ExpressionOptions {
        max_nodes,
        max_depth,
    };
    if shared {
        if branch_rules.is_some() {
            return Err(PyValueError::new_err(
                "shared=True constructs all branches; apply select_branch afterwards instead of supplying branch_rules",
            ));
        }
        let series = oneloop::get_expression_shared_with_options(master.expr.as_view(), options)
            .map_err(PyValueError::new_err)?;
        let expressions = series
            .into_coefficients()
            .into_iter()
            .map(|expression| Py::new(py, PythonSharedExpression { expression }))
            .collect::<PyResult<Vec<_>>>()?;
        return if let Some(index) = coefficient {
            Ok(expressions[index].clone_ref(py).into_any())
        } else {
            Ok(PyTuple::new(py, expressions)?.unbind().into_any())
        };
    }
    let series = if let Some(rules) = branch_rules {
        oneloop::get_expression_on_branch_with_options(
            master.expr.as_view(),
            &native_rules(rules)?,
            options,
        )
    } else {
        oneloop::get_expression_with_options(master.expr.as_view(), options)
    }
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

/// Select native if branches using replacements exclusively in conditions.
///
/// Accept one Expression or a tuple/list of Expressions. Return the same
/// container shape. Kept branches remain parametric; unresolved predicates
/// retain their original symbolic form. Every nested if is visited.
#[cfg(feature = "community")]
#[pyfunction(signature = (expression, replacement_rules, *, max_nodes=100_000_000, max_depth=4096))]
fn select_branch(
    expression: &Bound<'_, PyAny>,
    replacement_rules: Vec<PythonReplacement>,
    max_nodes: usize,
    max_depth: usize,
) -> PyResult<Py<PyAny>> {
    let py = expression.py();
    let replacement_rules = native_rules(replacement_rules)?;
    let select_one = |value: &Bound<'_, PyAny>| -> PyResult<Py<PythonExpression>> {
        if let Ok(value) = value.extract::<PyRef<'_, PythonSharedExpression>>() {
            let selected = oneloop::select_branch_shared(
                &value.expression,
                &replacement_rules,
                oneloop::ExpressionOptions {
                    max_nodes,
                    max_depth,
                },
            )
            .map_err(PyValueError::new_err)?;
            return Py::new(py, PythonExpression::from(selected));
        }
        let value = value.extract::<PythonExpression>()?;
        let selected = oneloop::select_branch(value.expr.as_view(), &replacement_rules);
        Py::new(py, PythonExpression::from(selected))
    };
    if let Ok(values) = expression.cast::<PyTuple>() {
        let selected = values
            .iter()
            .map(|v| select_one(&v))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyTuple::new(py, selected)?.unbind().into_any())
    } else if let Ok(values) = expression.cast::<PyList>() {
        let selected = values
            .iter()
            .map(|v| select_one(&v))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, selected)?.unbind().into_any())
    } else {
        Ok(select_one(expression)?.into_any())
    }
}

#[cfg(feature = "community")]
fn native_rules(replacement_rules: Vec<PythonReplacement>) -> PyResult<Vec<Replacement>> {
    // PythonReplacement's field is private. Symbolica's public Transformer
    // bridge exposes the native rules without reimplementing matching or
    // discarding rule conditions, callbacks, or settings.
    let mut transformer =
        PythonTransformer::new().replace_multiple(replacement_rules, false, false, false)?;
    let Some(Transformer::ReplaceAllMultiple(replacement_rules, _)) = transformer.chain.pop()
    else {
        return Err(PyValueError::new_err(
            "Symbolica did not construct a replacement transformer",
        ));
    };
    Ok(replacement_rules)
}

pub(super) fn register(module: &pyo3::Bound<'_, pyo3::types::PyModule>) -> pyo3::PyResult<()> {
    #[cfg(feature = "community")]
    {
        module.add_function(wrap_pyfunction!(get_expression, module)?)?;
        module.add_function(wrap_pyfunction!(select_branch, module)?)?;
        module.add_class::<PythonSharedExpression>()?;
    }
    #[cfg(not(feature = "community"))]
    let _ = module;
    Ok(())
}
