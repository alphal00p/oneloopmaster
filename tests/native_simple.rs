//! Algebraically compact native sectors versus the independent mapped formulas.
use oneloop::{EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, Float, Real, SingleFloat};

#[test]
fn compact_massless_sectors_match_mapped_expressions() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let bits = 320;
            let number = |s: &str| Float::parse(s, Some(bits)).unwrap();
            let complex = |r: Float| Complex::new(r, Float::new(bits));
            let zero = || complex(Float::new(bits));
            for family in [ScalarIntegral::B0, ScalarIntegral::C0, ScalarIntegral::D0] {
                let mut native =
                    NativeEvaluator::<Float>::with_binary_precision(family, bits).unwrap();
                let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
                    family,
                    bits,
                    EvaluationBackend::Expression,
                )
                .unwrap();
                let mut cases = Vec::new();
                for mu in ["1e-50", "1", "1e50"] {
                    if family == ScalarIntegral::B0 {
                        cases.push(vec![zero(), zero(), zero(), complex(number(mu))]);
                    } else if family == ScalarIntegral::C0 {
                        for p in [
                            ["0", "0", "0"],
                            ["0", "0", "-2"],
                            ["0", "0", "2"],
                            ["0", "-2", "3"],
                            ["0", "2", "3"],
                            ["0", "-2", "-3"],
                            ["0", "2", "2"],
                            ["0", "-2", "-2"],
                            ["0", "2", "2.000000000000000000000000000001"],
                        ] {
                            for rotation in 0..3 {
                                let mut row = (0..3)
                                    .map(|i| complex(number(p[(i + rotation) % 3])))
                                    .collect::<Vec<_>>();
                                row.extend([zero(), zero(), zero(), complex(number(mu))]);
                                cases.push(row);
                            }
                        }
                    } else {
                        for s in ["-3", "2"] {
                            for t in ["-5", "2", "4"] {
                                let mut row = (0..4).map(|_| zero()).collect::<Vec<_>>();
                                row.extend([complex(number(s)), complex(number(t))]);
                                row.extend((0..4).map(|_| zero()));
                                row.push(complex(number(mu)));
                                cases.push(row);
                            }
                        }
                    }
                }
                let mut actual = std::array::from_fn::<_, 3, _>(|_| zero());
                let mut expected = actual.clone();
                for (index, row) in cases.iter().enumerate() {
                    native.evaluate(row, &mut actual).unwrap();
                    expression.evaluate(row, &mut expected).unwrap();
                    for tag in 0..3 {
                        for (a, e) in [
                            (&actual[tag].re, &expected[tag].re),
                            (&actual[tag].im, &expected[tag].im),
                        ] {
                            assert!(
                                a.is_finite() && e.is_finite(),
                                "{family:?} row{index} tag{tag}: nonfinite"
                            );
                            let scale = if e.is_zero() { number("1") } else { e.norm() };
                            assert!(
                                (a.clone() - e).norm() <= scale * number("1e-60"),
                                "{family:?} row{index} tag{tag}: compact {a:e}, mapped {e:e}"
                            );
                        }
                    }
                }
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
