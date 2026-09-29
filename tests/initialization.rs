//! Fresh-process coverage for lazy first-entry paths and explicit warmup.
#[path = "support/fixtures.rs"]
mod fixtures;

#[test]
fn first_entries_leave_full_warmup_opt_in() {
    for entry in [
        "symbol",
        "master",
        "context",
        "manual",
        "native",
        "precision",
    ] {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--ignored",
                "--exact",
                "initialization_child",
                "--nocapture",
                "--test-threads=1",
            ])
            .env("ONELOOP_INITIALIZATION_TEST_ENTRY", entry)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{entry} startup failed: {}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        for log in [&output.stdout, &output.stderr] {
            let log = String::from_utf8_lossy(log);
            assert!(!log.contains("Created infinity"), "{entry}: {log}");
            assert!(!log.contains("Created indeterminate"), "{entry}: {log}");
        }
    }
}

#[test]
#[ignore = "fresh-process helper invoked by first_entries_leave_full_warmup_opt_in"]
fn initialization_child() {
    let entry = std::env::var("ONELOOP_INITIALIZATION_TEST_ENTRY")
        .expect("run through the parent initialization regression");
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(move || {
            use oneloop::{ScalarEvaluator, ScalarIntegral};
            use symbolica::prelude::*;
            assert!(!oneloop::is_initialized());
            match entry.as_str() {
                "symbol" => {
                    // No OneLOop accessor has been called. The State initializer
                    // must attach all hooks before this parsed symbol is used.
                    let _ = parse!("oneloopmaster::A0(0,x,1)");
                }
                "master" => {
                    let _ = oneloop::A0();
                }
                "context" => {
                    let _ = oneloop::OneLoopExpressions::new();
                }
                "manual" => {
                    let _ = ScalarEvaluator::cached(ScalarIntegral::A0).unwrap();
                }
                "native" => {
                    let mut output = [Complex::new(0., 0.); 3];
                    oneloop::evaluate(
                        ScalarIntegral::A0,
                        &[Complex::new(2., 0.), Complex::new(1., 0.)],
                        &mut output,
                    )
                    .unwrap();
                    assert!((output[0].re - 2. * (1. - 2_f64.ln())).abs() < 1e-14);
                }
                "precision" => {
                    let _ = oneloop::PrecisionEvaluator::new(ScalarIntegral::A0, 50).unwrap();
                }
                _ => panic!("unknown initialization entry"),
            }
            assert!(!oneloop::is_initialized());
            oneloop::initialize().unwrap();
            assert!(oneloop::is_initialized());
            oneloop::initialize().unwrap(); // Explicit warmup remains idempotent.
            let rows = fixtures::parse(include_str!("data/benchmark.txt"));
            for family in [
                ScalarIntegral::A0,
                ScalarIntegral::B0,
                ScalarIntegral::DB0,
                ScalarIntegral::C0,
                ScalarIntegral::D0,
            ] {
                let row = rows.iter().find(|row| row.family == family).unwrap();
                let mut result = [Complex::new(f64::NAN, f64::NAN); 3];
                oneloop::evaluate(family, &row.args, &mut result).unwrap();
                row.check(&result);
                let mut clone = ScalarEvaluator::cached(family).unwrap();
                clone.evaluate(&row.args, &mut result).unwrap();
                row.check(&result);
                let mut batch = [Complex::new(f64::NAN, f64::NAN); 15];
                oneloop::evaluate_batch(family, &row.args.repeat(5), &mut batch, 5).unwrap();
                for result in batch.chunks_exact(3) {
                    row.check(result);
                }
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
