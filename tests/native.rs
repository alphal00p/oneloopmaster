//! Direct Rust arithmetic versus original fixtures and native expressions.
#[path = "support/fixtures.rs"]
mod fixtures;
use oneloop::{EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, DoubleFloat, Float, Real, RealLike, SingleFloat};

const FAMILIES: [ScalarIntegral; 5] = [
    ScalarIntegral::A0,
    ScalarIntegral::B0,
    ScalarIntegral::DB0,
    ScalarIntegral::C0,
    ScalarIntegral::D0,
];

fn on_stack(action: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(action)
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn native_f64_acceptance_scalar_and_mixed_batches() {
    on_stack(|| {
        let fixtures = fixtures::parse(include_str!("data/parity.txt"));
        assert_eq!(fixtures.len(), 293);
        let mut failures = Vec::new();
        for family in FAMILIES {
            let mut evaluator = NativeEvaluator::<f64>::new(family).unwrap();
            let rows = fixtures
                .iter()
                .filter(|row| row.family == family)
                .collect::<Vec<_>>();
            let mut single = [Complex::new(0., 0.); 3];
            for row in &rows {
                evaluator.evaluate(&row.args, &mut single).unwrap();
                failures.extend(row.failures(&single));
            }
            let input = rows
                .iter()
                .flat_map(|r| r.args.iter().copied())
                .collect::<Vec<_>>();
            for batch in [1, 3, 4, 5, 31, 256, 1024] {
                let mut output = vec![Complex::new(f64::NAN, f64::NAN); rows.len() * 3];
                for (points, values) in input
                    .chunks(family.arity() * batch)
                    .zip(output.chunks_mut(3 * batch))
                {
                    evaluator
                        .evaluate_batch(points, values, points.len() / family.arity())
                        .unwrap();
                }
                for (row, values) in rows.iter().zip(output.chunks_exact(3)) {
                    failures.extend(
                        row.failures(values)
                            .into_iter()
                            .map(|f| format!("batch={batch}: {f}")),
                    );
                }
            }
        }
        assert!(
            failures.is_empty(),
            "{} native comparisons failed:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}

#[test]
fn native_double_float_acceptance() {
    on_stack(|| {
        let fixtures = fixtures::parse(include_str!("data/parity.txt"));
        let mut failures = Vec::new();
        for family in FAMILIES {
            let mut evaluator = NativeEvaluator::<DoubleFloat>::new(family).unwrap();
            for row in fixtures.iter().filter(|row| row.family == family) {
                let input = row
                    .args
                    .iter()
                    .map(|v| Complex::new(v.re.into(), v.im.into()))
                    .collect::<Vec<_>>();
                let mut output =
                    std::array::from_fn::<_, 3, _>(|_| Complex::new(0.0.into(), 0.0.into()));
                evaluator.evaluate(&input, &mut output).unwrap();
                let converted = output.map(|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
                failures.extend(row.failures(&converted));
            }
        }
        assert!(
            failures.is_empty(),
            "{} DoubleFloat comparisons failed:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}

#[test]
fn native_expanded_fortran_fixtures_at_256_bits() {
    on_stack(|| {
        let fixtures = fixtures::parse(include_str!("data/scalar_audit.txt"));
        assert_eq!(fixtures.len(), 342);
        let mut failures = Vec::new();
        for family in FAMILIES {
            let mut evaluator =
                NativeEvaluator::<Float>::with_binary_precision(family, 256).unwrap();
            for row in fixtures.iter().filter(|row| row.family == family) {
                // The oracle table is binary64. Preserve precisely those inputs
                // for this comparison; this is not a 256-bit-accuracy oracle.
                let input = row
                    .args
                    .iter()
                    .map(|z| Complex::new(Float::with_val(256, z.re), Float::with_val(256, z.im)))
                    .collect::<Vec<_>>();
                let mut output = std::array::from_fn::<_, 3, _>(|_| {
                    Complex::new(Float::new(256), Float::new(256))
                });
                evaluator.evaluate(&input, &mut output).unwrap();
                let converted = output.map(|z| Complex::new(z.re.to_f64(), z.im.to_f64()));
                failures.extend(row.failures(&converted));
            }
        }
        assert!(
            failures.is_empty(),
            "{} expanded native comparisons failed:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}

#[test]
#[ignore = "exploratory binary64 stress audit; report failures without accepting them"]
fn native_expanded_binary64_diagnostic() {
    on_stack(|| {
        let fixtures = fixtures::parse(include_str!("data/scalar_audit.txt"));
        let mut failures = Vec::new();
        for family in FAMILIES {
            let mut evaluator = NativeEvaluator::<f64>::new(family).unwrap();
            for row in fixtures.iter().filter(|row| row.family == family) {
                let mut output = [Complex::new(0., 0.); 3];
                evaluator.evaluate(&row.args, &mut output).unwrap();
                failures.extend(row.failures(&output));
            }
        }
        assert!(
            failures.is_empty(),
            "{} expanded native comparisons failed:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}

#[test]
fn native_arbitrary_precision_matches_expressions_at_1000_digits() {
    on_stack(|| {
        for (family, values) in [
            (ScalarIntegral::A0, vec![("2", "-0.125"), ("4", "0")]),
            (
                ScalarIntegral::B0,
                vec![("-1", "0"), ("1", "-0.1"), ("2", "-0.2"), ("4", "0")],
            ),
            (
                ScalarIntegral::DB0,
                vec![("-1", "0"), ("1", "-0.1"), ("2", "-0.2"), ("4", "0")],
            ),
            (
                ScalarIntegral::C0,
                vec![
                    ("-1", "0"),
                    ("-2", "0"),
                    ("-3", "0"),
                    ("1", "-0.1"),
                    ("2", "-0.2"),
                    ("3", "-0.3"),
                    ("4", "0"),
                ],
            ),
            (
                ScalarIntegral::D0,
                vec![
                    ("-1", "0"),
                    ("-2", "0"),
                    ("-3", "0"),
                    ("-4", "0"),
                    ("-5", "0"),
                    ("-6", "0"),
                    ("1", "-0.1"),
                    ("2", "-0.2"),
                    ("3", "-0.3"),
                    ("4", "-0.4"),
                    ("4", "0"),
                ],
            ),
        ] {
            let mut reference =
                PrecisionEvaluator::with_backend(family, 1120, EvaluationBackend::Expression)
                    .unwrap();
            let bits = 3456;
            let mut native = NativeEvaluator::<Float>::with_binary_precision(family, bits).unwrap();
            let input = values
                .iter()
                .map(|(r, i)| {
                    Complex::new(
                        Float::parse(r, Some(reference.binary_precision())).unwrap(),
                        Float::parse(i, Some(reference.binary_precision())).unwrap(),
                    )
                })
                .collect::<Vec<_>>();
            let zero = || Complex::new(Float::new(bits), Float::new(bits));
            let mut actual = std::array::from_fn::<_, 3, _>(|_| zero());
            let mut expected = actual.clone();
            native.evaluate(&input, &mut actual).unwrap();
            reference.evaluate(&input, &mut expected).unwrap();
            let tolerance = Float::parse("1e-1000", Some(bits)).unwrap();
            for (tag, (a, e)) in actual.iter().zip(&expected).enumerate() {
                for (component, (a, e)) in [(&a.re, &e.re), (&a.im, &e.im)].into_iter().enumerate()
                {
                    assert!(
                        a.is_finite(),
                        "{} tag{tag} component{component} nonfinite",
                        family.name()
                    );
                    let error = (a.clone() - e).norm();
                    let scale = if e.is_zero() {
                        Float::with_val(bits, 1)
                    } else {
                        e.norm()
                    };
                    assert!(
                        error < scale * &tolerance,
                        "{} tag{tag} component{component}: error={error:e}",
                        family.name()
                    );
                }
            }
        }
    });
}

#[test]
fn native_shape_and_domain_errors_do_not_modify_outputs() {
    let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::B0).unwrap();
    let mut output = [Complex::new(42., -17.); 6];
    let sentinel = output;
    let good = [
        Complex::new(-1., 0.),
        Complex::new(1., 0.),
        Complex::new(2., 0.),
        Complex::new(4., 0.),
    ];
    assert!(evaluator.evaluate(&good[..3], &mut output[..3]).is_err());
    assert!(evaluator.evaluate_batch(&[], &mut [], usize::MAX).is_err());
    evaluator.evaluate_batch(&[], &mut [], 0).unwrap();
    let mut batch = good.to_vec();
    let mut invalid = good;
    invalid[1].im = f64::from_bits(1);
    batch.extend(invalid);
    assert!(evaluator.evaluate_batch(&batch, &mut output, 2).is_err());
    assert_eq!(output, sentinel);
}
