//! Eager registration and backend preparation at Symbolica state startup.
use std::sync::OnceLock;

static INITIALIZED: OnceLock<Result<(), String>> = OnceLock::new();

/// Trigger Symbolica's initializer inventory before entering any OneLOop
/// OnceLock. Its initializer may call back into the same public constructors.
pub(crate) fn ensure_symbolica_state() {
    if INITIALIZED.get().is_some() {
        return;
    }
    let _ = symbolica::symbol!("oneloopmaster::__state_initialized");
}

pub(crate) fn from_symbolica() {
    INITIALIZED.get_or_init(|| {
        // Register every master, native helper and argument symbol together.
        // The transparent definitions are also ready for manual compilation.
        let _ = crate::expressions::shared_definitions();
        crate::backend::initialize_native_all()?;
        crate::evaluators::initialize_all()
    });
}

/// Initialize every OneLOop symbol and all five shared numerical backends.
///
/// The `initialize!` registration also runs automatically with Symbolica's
/// global-state initialization. Call this function at application startup to
/// receive configuration/cache errors before evaluating anything. Python calls
/// it during module import. Work stays on the calling thread; see the README's
/// stack and Symbolica-license requirements.
pub fn initialize() -> Result<(), String> {
    if let Some(status) = INITIALIZED.get() {
        return status.clone();
    }
    crate::evaluators::portable_environment()?;
    ensure_symbolica_state();
    INITIALIZED
        .get()
        .ok_or_else(|| "OneLOop initialization was requested recursively".to_owned())?
        .clone()
}

/// Whether symbol registration and all five shared backends are complete.
/// This query does not itself initialize Symbolica or load evaluators.
pub fn is_initialized() -> bool {
    INITIALIZED.get().is_some_and(Result::is_ok)
}
