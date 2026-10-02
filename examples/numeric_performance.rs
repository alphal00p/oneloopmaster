//! Compare native and interpreted throughput in all three numeric domains.
//! Usage: numeric_performance [iterations=100] [repetitions=5] [native|expression]
//! Construction, input conversion and fixture validation are outside timing.
#[path = "../tests/support/fixtures.rs"]
mod fixtures;
use oneloop::{NativeEvaluator, NativeFloat, OneLoopExpressions, ScalarIntegral};
use std::{hint::black_box, time::Instant};
use symbolica::{
    atom::Atom,
    domains::float::{Complex, DoubleFloat, Float, RealLike},
    evaluate::EvaluationDomain,
};

fn survey<T: NativeFloat + RealLike>(
    name: &str,
    iterations: usize,
    repetitions: usize,
    native: bool,
    convert: impl Fn(f64) -> T,
) where
    Complex<T>: EvaluationDomain,
{
    let fixtures = fixtures::parse(include_str!("../tests/data/benchmark.txt"));
    for (family, symbol) in [
        (ScalarIntegral::A0, oneloop::A0()),
        (ScalarIntegral::B0, oneloop::B0()),
        (ScalarIntegral::DB0, oneloop::dB0()),
        (ScalarIntegral::C0, oneloop::C0()),
        (ScalarIntegral::D0, oneloop::D0()),
    ] {
        let prototype = T::zero_at(T::DEFAULT_BITS);
        let rows = fixtures
            .iter()
            .filter(|r| r.family == family)
            .collect::<Vec<_>>();
        let inputs = rows
            .iter()
            .map(|r| {
                r.args
                    .iter()
                    .map(|z| Complex::new(convert(z.re), convert(z.im)))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut output =
            core::array::from_fn::<_, 3, _>(|_| Complex::new(prototype.zero(), prototype.zero()));
        let mut evaluate: Box<dyn FnMut(&[Complex<T>], &mut [Complex<T>])> = if native {
            let mut evaluator = NativeEvaluator::<T>::new(family).unwrap();
            Box::new(move |input, output| evaluator.evaluate(input, output).unwrap())
        } else {
            let parameters = (0..family.arity())
                .map(|i| {
                    symbolica::symbol!(format!("numeric_survey_{}_{i}", family.name())).to_atom()
                })
                .collect::<Vec<_>>();
            let calls = [0, -1, -2].map(|tag| {
                let mut args = vec![Atom::num(tag)];
                args.extend(parameters.iter().cloned());
                symbol.call(&args)
            });
            let mut evaluator = OneLoopExpressions::new()
                .evaluator(&calls, &parameters)
                .unwrap()
                .map_coeff_with_prec(
                    &|z| {
                        Complex::new(
                            prototype.from_rational(&z.re),
                            prototype.from_rational(&z.im),
                        )
                    },
                    T::DEFAULT_BITS,
                );
            Box::new(move |input, output| evaluator.evaluate(input, output))
        };
        for (row, input) in rows.iter().zip(&inputs) {
            evaluate(input, &mut output);
            row.check(
                &output
                    .clone()
                    .map(|z| Complex::new(z.re.to_f64(), z.im.to_f64())),
            );
        }
        let mut times = Vec::new();
        for _ in 0..repetitions {
            let start = Instant::now();
            for _ in 0..iterations {
                for input in &inputs {
                    evaluate(black_box(input), black_box(&mut output));
                    black_box(&output);
                }
            }
            times.push(start.elapsed().as_secs_f64() * 1e9 / (iterations * inputs.len()) as f64);
        }
        times.sort_by(f64::total_cmp);
        println!("{name}\t{}\t{:.1}", family.name(), times[times.len() / 2]);
    }
}

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let iterations = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(100);
    let repetitions = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(5);
    assert!(iterations > 0 && repetitions > 0);
    let native = match args.get(3).map(String::as_str).unwrap_or("native") {
        "native" => true,
        "expression" => false,
        _ => panic!("backend must be native or expression"),
    };
    println!("domain\tfamily\tmedian_ns_per_point");
    survey::<f64>("f64", iterations, repetitions, native, |v| v);
    survey::<DoubleFloat>(
        "DoubleFloat",
        iterations,
        repetitions,
        native,
        DoubleFloat::from,
    );
    survey::<Float>("Float128", iterations, repetitions, native, |v| {
        Float::with_val(128, v)
    });
}
