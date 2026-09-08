//! `cargo run --release --example arbitrary_precision -- 1000`
use oneloop::{PrecisionEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, Float};

fn main() -> Result<(), String> {
    let digits = std::env::args()
        .nth(1)
        .map(|text| text.parse::<u32>().map_err(|error| error.to_string()))
        .transpose()?
        .unwrap_or(32);
    std::thread::Builder::new()
        .name("oneloop-precision-example".into())
        .stack_size(128 * 1024 * 1024)
        .spawn(move || -> Result<(), String> {
            // Construct Symbolica values only after entering the large-stack
            // calling thread. No second active Symbolica worker is introduced.
            let mut evaluator = PrecisionEvaluator::new(ScalarIntegral::B0, digits)?;
            let bits = evaluator.binary_precision();
            let input = ["-1", "1", "1", "1", "-1", "1", "1", "4"]
                .into_iter()
                .map(|text| {
                    Float::parse(text, Some(bits)).map(|real| Complex::new(real, Float::new(bits)))
                })
                .collect::<Result<Vec<_>, _>>()?;
            let mut output = (0..6)
                .map(|_| Complex::new(Float::new(bits), Float::new(bits)))
                .collect::<Vec<_>>();
            evaluator.evaluate_batch(&input, &mut output, 2)?;
            println!("Requested {digits} decimal digits; {bits} working bits.");
            for (row, values) in output.chunks_exact(3).enumerate() {
                println!("B0 row {row} [finite, simple pole, double pole]:");
                for value in values {
                    println!("  {} + ({}) i", value.re, value.im);
                }
            }
            Ok(())
        })
        .map_err(|error| error.to_string())?
        .join()
        .map_err(|_| "precision example thread failed".to_owned())?
}
