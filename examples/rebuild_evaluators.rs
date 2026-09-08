//! Rebuild portable SymJIT O2 assets without requiring existing embedded blobs.
//! cargo run --release --no-default-features --example rebuild_evaluators
#[path = "../tests/support/fixtures.rs"]
mod fixtures;
use oneloop::{ScalarEvaluator, ScalarIntegral};
use std::{path::PathBuf, time::Instant};
use symbolica::domains::float::Complex;

fn main() {
    let mut args = std::env::args().skip(1);
    let destination = args.next().unwrap_or_else(|| "assets/evaluators".into());
    let selected = args.next();
    assert!(
        args.next().is_none(),
        "usage: rebuild_evaluators [directory] [family]"
    );
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(move || {
            let directory = PathBuf::from(destination);
            std::fs::create_dir_all(&directory)?;
            let mut fixtures = fixtures::parse(include_str!("../tests/data/parity.txt"));
            fixtures.extend(fixtures::parse(include_str!("../tests/data/benchmark.txt")));
            let mut count = 0;
            let mut failures = Vec::new();
            for (family, file) in [
                (ScalarIntegral::A0, "a0.bin"),
                (ScalarIntegral::B0, "b0.bin"),
                (ScalarIntegral::DB0, "db0.bin"),
                (ScalarIntegral::C0, "c0.bin"),
                (ScalarIntegral::D0, "d0.bin"),
            ] {
                if selected
                    .as_ref()
                    .is_some_and(|s| !s.eq_ignore_ascii_case(family.name()))
                {
                    continue;
                }
                let started = Instant::now();
                eprintln!("building {}", family.name());
                let evaluator = ScalarEvaluator::rebuild(family).map_err(std::io::Error::other)?;
                let data = evaluator.to_bytes().map_err(std::io::Error::other)?;
                let mut restored =
                    ScalarEvaluator::from_bytes(family, &data).map_err(std::io::Error::other)?;
                let previous_failures = failures.len();
                for row in fixtures.iter().filter(|row| row.family == family) {
                    let mut result = [Complex::new(f64::NAN, f64::NAN); 3];
                    restored
                        .evaluate(&row.args, &mut result)
                        .map_err(std::io::Error::other)?;
                    failures.extend(row.failures(&result));
                }
                count += 1;
                if failures.len() != previous_failures {
                    eprintln!(
                        "{}: validation failed; existing asset left untouched",
                        family.name()
                    );
                    continue;
                }
                let path = directory.join(file);
                std::fs::write(&path, &data)?;
                eprintln!(
                    "{}: {} bytes, {:?}, {}",
                    family.name(),
                    data.len(),
                    started.elapsed(),
                    path.display()
                );
            }
            if count == 0 {
                return Err(std::io::Error::other("unknown scalar family"));
            }
            if !failures.is_empty() {
                return Err(std::io::Error::other(failures.join("\n")));
            }
            Ok::<_, std::io::Error>(())
        })
        .expect("start evaluator builder")
        .join()
        .expect("evaluator builder panicked")
        .expect("could not rebuild evaluators");
}
