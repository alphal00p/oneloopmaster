//! Fresh-process coverage for requesting transcendental symbols before State.
use std::sync::atomic::{AtomicBool, Ordering};

static DEPENDENT_INITIALIZER_RAN: AtomicBool = AtomicBool::new(false);

// Exercise accessor reentry during State initialization. Upstream must enter
// State before locking its symbol caches, including when OneLOop is linked.
symbolica::initialize!(
    || {
        assert_eq!(
            symbolica::transcendental::polylog().get_name(),
            "symbolica::polylog"
        );
        assert_eq!(
            symbolica::transcendental::tan().get_name(),
            "symbolica::tan"
        );
        assert_eq!(
            symbolica::transcendental::bessel_j().get_name(),
            "symbolica::bessel_j"
        );
        DEPENDENT_INITIALIZER_RAN.store(true, Ordering::Release);
    },
    "oneloop"
);

#[test]
fn first_transcendental_call_initializes_linked_inventory() {
    use std::{
        process::{Command, Stdio},
        time::{Duration, Instant},
    };
    for entry in ["polylog", "geometric", "bessel"] {
        let mut child = Command::new(std::env::current_exe().unwrap())
            .args([
                "--ignored",
                "--exact",
                "native_startup_child",
                "--nocapture",
                "--test-threads=1",
            ])
            .env("ONELOOP_NATIVE_STARTUP_ENTRY", entry)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let start = Instant::now();
        let timed_out = loop {
            if child.try_wait().unwrap().is_some() {
                break false;
            }
            if start.elapsed() > Duration::from_secs(180) {
                child.kill().expect("terminate the timed-out startup child");
                break true;
            }
            std::thread::sleep(Duration::from_millis(20));
        };
        let output = child.wait_with_output().unwrap();
        assert!(
            !timed_out && output.status.success(),
            "{entry} startup failed (timeout={timed_out}): {}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
#[ignore = "fresh-process helper; the parent supplies a timeout"]
fn native_startup_child() {
    let entry = std::env::var("ONELOOP_NATIVE_STARTUP_ENTRY")
        .expect("invoke through the startup parent test");
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(move || {
            use symbolica::{atom::Atom, transcendental::TranscendentalFunctions};
            assert!(!oneloop::is_initialized());
            // No initialize(), parse!(), symbol!(), or evaluator precedes this call.
            let expression = match entry.as_str() {
                "polylog" => Atom::num((1, 2)).polylog(2),
                "geometric" => symbolica::transcendental::tan().call(1),
                "bessel" => symbolica::transcendental::bessel_j().call((0, 1)),
                _ => panic!("unknown first-entry route"),
            };
            assert!(!expression.is_zero());
            assert!(DEPENDENT_INITIALIZER_RAN.load(Ordering::Acquire));
            assert!(oneloop::is_initialized());
            for family in [
                oneloop::ScalarIntegral::A0,
                oneloop::ScalarIntegral::B0,
                oneloop::ScalarIntegral::DB0,
                oneloop::ScalarIntegral::C0,
                oneloop::ScalarIntegral::D0,
            ] {
                let _ = oneloop::ScalarEvaluator::cached(family).unwrap();
            }
            oneloop::initialize().unwrap(); // Already complete and idempotent.
        })
        .unwrap()
        .join()
        .unwrap();
}
