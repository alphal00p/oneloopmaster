//! Lightweight symbol registration and optional explicit backend warmup.
use std::sync::OnceLock;

static REGISTERED: OnceLock<()> = OnceLock::new();
static INITIALIZED: OnceLock<Result<(), String>> = OnceLock::new();

/// Trigger Symbolica's initializer inventory before entering any OneLOop
/// OnceLock. Its initializer may call back into the same public constructors.
pub(crate) fn ensure_symbolica_state() {
    if REGISTERED.get().is_some() {
        return;
    }
    let _ = symbolica::symbol!("oneloopmaster::__state_initialized");
}

pub(crate) fn from_symbolica() {
    REGISTERED.get_or_init(|| {
        // Parsed masters need their numeric hooks before callers can create
        // symbols with these names. Keep formulas and evaluators out of startup.
        let _ = (
            crate::A0(),
            crate::B0(),
            crate::dB0(),
            crate::C0(),
            crate::D0(),
        );
        crate::definitions::register_callbacks();
    });
}

/// Initialize every OneLOop symbol and all five shared numerical backends.
///
/// This is an optional eager warmup. Ordinary imports and evaluations prepare
/// only the requested backend on first use. Call this function explicitly to
/// receive configuration/cache errors before evaluating anything.
/// Work stays on the calling thread; see the README's
/// stack and Symbolica-license requirements.
pub fn initialize() -> Result<(), String> {
    if let Some(status) = INITIALIZED.get() {
        return status.clone();
    }
    crate::evaluators::portable_environment()?;
    ensure_symbolica_state();
    INITIALIZED
        .get_or_init(|| {
            let _ = crate::expressions::shared_definitions();
            crate::inspection::prepare_symbols();
            crate::backend::initialize_native_all()?;
            crate::evaluators::initialize_all()
        })
        .clone()
}

/// Whether explicit full warmup via [`initialize`] completed successfully.
/// This query does not itself initialize Symbolica or load evaluators.
pub fn is_initialized() -> bool {
    INITIALIZED.get().is_some_and(Result::is_ok)
}
