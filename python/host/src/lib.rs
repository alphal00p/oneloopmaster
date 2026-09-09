//! Minimal development host: one extension and one Symbolica state.
use pyo3::{prelude::*, types::PyModule};
use symbolica::api::python::{SymbolicaCommunityModule, create_symbolica_module};

#[pymodule]
fn core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    create_symbolica_module(module)?;
    let oneloop = PyModule::new(module.py(), "oneloop_native")?;
    oneloop_native::CommunityModule::register_module(&oneloop)?;
    module.add_submodule(&oneloop)?;
    module
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item("symbolica.community.oneloop_native", &oneloop)?;
    // Eagerly initialize before any Python expression or numeric call is served.
    oneloop_native::CommunityModule::initialize(module.py())?;
    Ok(())
}
