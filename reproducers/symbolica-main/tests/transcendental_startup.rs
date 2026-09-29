use std::sync::atomic::{AtomicBool, Ordering};

static INITIALIZED: AtomicBool = AtomicBool::new(false);

symbolica::initialize!(|| {
    // A downstream initializer can legitimately use the registered accessors.
    let _ = symbolica::transcendental::polylog();
    let _ = symbolica::transcendental::tan();
    let _ = symbolica::transcendental::bessel_j();
    INITIALIZED.store(true, Ordering::Release);
});

#[test]
fn cold_accessors_allow_initializer_reentry() {
    use std::{
        process::Command,
        time::{Duration, Instant},
    };
    let mut failures = Vec::new();
    for entry in ["polylog", "tan", "bessel_j"] {
        let mut child = Command::new(std::env::current_exe().unwrap())
            .args(["--ignored", "--exact", "cold_accessor_child", "--nocapture"])
            .env("SYMBOLICA_STARTUP_REGRESSION", entry)
            .spawn()
            .unwrap();
        let start = Instant::now();
        loop {
            if let Some(status) = child.try_wait().unwrap() {
                if !status.success() {
                    failures.push(format!("{entry}: {status}"));
                }
                break;
            }
            if start.elapsed() > Duration::from_secs(10) {
                child.kill().unwrap();
                child.wait().unwrap();
                failures.push(format!("{entry}: initialization deadlocked"));
                break;
            }
            std::thread::sleep(Duration::from_millis(20));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
#[ignore = "fresh-process helper invoked by cold_accessors_allow_initializer_reentry"]
fn cold_accessor_child() {
    // Do not initialize State before requesting the first accessor.
    let entry = std::env::var("SYMBOLICA_STARTUP_REGRESSION").unwrap();
    let first = match entry.as_str() {
        "polylog" => symbolica::transcendental::polylog(),
        "tan" => symbolica::transcendental::tan(),
        "bessel_j" => symbolica::transcendental::bessel_j(),
        _ => unreachable!(),
    };
    assert!(INITIALIZED.load(Ordering::Acquire));
    assert_eq!(
        Some(first),
        symbolica::get_symbol!(format!("symbolica::{entry}"))
    );
}
