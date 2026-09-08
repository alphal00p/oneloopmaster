//! Regression inputs come from the original Fortran corpus. Precision tests
//! use the same inputs in multiple native Symbolica numeric domains; there is
//! no numerical OneLOop fallback or altered comparison tolerance.
use oneloop::{a0, b0, db0};
use symbolica::prelude::*;

fn exact(
    series: &oneloop::LaurentSeries,
    parameters: &[Atom],
) -> symbolica::evaluate::ExpressionEvaluator<Complex<Rational>> {
    let views = series.coefficients().each_ref().map(Atom::as_view);
    Atom::evaluator_multiple(&views, parameters)
        .build()
        .unwrap()
}

#[test]
fn lower_lip_is_independent_of_signed_zero() {
    let m = parse!("axis_m");
    let mu = Atom::num(1);
    let formula = a0(&m, &mu);
    let mut evaluator =
        exact(&formula, &[m]).map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    for imaginary in [0.0, -0.0, -1e-12, -1e-24] {
        let mut result = [Complex::new(0., 0.); 3];
        evaluator.evaluate(&[Complex::new(-2., imaginary)], &mut result);
        assert!((result[0].re + 0.6137056388801094).abs() < 1e-10);
        assert!((result[0].im + 2. * std::f64::consts::PI).abs() < 1e-10);
    }
    // Exact construction and runtime substitution must agree on the same lip.
    let literal = a0(&Atom::num(-2), &mu);
    let mut evaluator =
        exact(&literal, &[]).map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let mut result = [Complex::new(0., 0.); 3];
    evaluator.evaluate(&[], &mut result);
    assert!((result[0].im + 2. * std::f64::consts::PI).abs() < 1e-12);
}

#[test]
fn derivative_coincident_roots_match_the_original_convention() {
    let params = [
        parse!("limit_s"),
        parse!("limit_m0"),
        parse!("limit_m1"),
        parse!("limit_mu"),
    ];
    let formula = db0(&params[0], &params[1], &params[2], &params[3]);
    let mut evaluator =
        exact(&formula, &params).map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let expected = 3. * 2_f64.ln() - 2.;
    for (m0, m1) in [(1., 4.), (4., 1.)] {
        for mu_squared in [1e-12, 1., 1e12] {
            let mut result = [Complex::new(0., 0.); 3];
            let input = [1., m0, m1, mu_squared].map(|v| Complex::new(v, 0.));
            evaluator.evaluate(&input, &mut result);
            assert!(
                (result[0].re - expected).hypot(result[0].im) < 1e-12,
                "{result:?}"
            );
            assert_eq!(result[1], Complex::new(0., 0.));
            assert_eq!(result[2], Complex::new(0., 0.));
        }
    }
    let mut result = [Complex::new(0., 0.); 3];
    evaluator.evaluate(
        &[-1., -1., -4., 1.].map(|v| Complex::new(v, 0.)),
        &mut result,
    );
    assert!((result[0].re + expected).hypot(result[0].im) < 1e-12);
    // OneLOop prescribes a finite value at coincident normal-threshold roots.
    // This is API parity, not a claim that the one-sided derivative is finite.
    evaluator.evaluate(&[4., 1., 1., 1.].map(|v| Complex::new(v, 0.)), &mut result);
    assert!((result[0].re + 0.5).hypot(result[0].im) < 1e-12);
}

#[test]
fn all_exploratory_bubbles_converge_at_high_precision() {
    let params = [
        parse!("prec_s"),
        parse!("prec_m0"),
        parse!("prec_m1"),
        parse!("prec_mu"),
    ];
    let originals = [
        exact(&b0(&params[0], &params[1], &params[2], &params[3]), &params),
        exact(
            &db0(&params[0], &params[1], &params[2], &params[3]),
            &params,
        ),
    ];
    let mut previous = Vec::<[Complex<Float>; 3]>::new();
    for bits in [256, 384] {
        let mut evaluators = originals.clone().map(|e| {
            e.map_coeff_with_prec(
                &|c| {
                    Complex::new(
                        c.re.to_multi_prec_float(bits),
                        c.im.to_multi_prec_float(bits),
                    )
                },
                bits,
            )
        });
        let mut count = 0;
        for (line_index, line) in include_str!("data/scalar_audit.txt").lines().enumerate() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let row = line
                .split_whitespace()
                .map(|v| v.parse::<f64>().unwrap())
                .collect::<Vec<_>>();
            let family = if row[0] == 2. {
                0
            } else if row[0] == -2. {
                1
            } else {
                continue;
            };
            assert_eq!(row.len(), 13);
            let input = [
                (row[2], 0.),
                (row[3], row[4]),
                (row[5], row[6]),
                (row[1], 0.),
            ]
            .map(|(re, im)| Complex::new(Float::with_val(bits, re), Float::with_val(bits, im)));
            let mut result = core::array::from_fn::<_, 3, _>(|_| {
                Complex::new(Float::new(bits), Float::new(bits))
            });
            evaluators[family].evaluate(&input, &mut result);
            let scale = row[2..7]
                .iter()
                .fold(0_f64, |s, v| s.max(v.abs()))
                .max(f64::MIN_POSITIVE);
            let normalization = if family == 0 { 1. } else { scale };
            for (k, (value, reference)) in result.iter().zip(row[7..].chunks_exact(2)).enumerate() {
                let error = (value.re.to_f64() - reference[0])
                    .hypot(value.im.to_f64() - reference[1])
                    * normalization;
                let tolerance = 2e-10 + 2e-8 * reference[0].hypot(reference[1]) * normalization;
                assert!(
                    error.is_finite() && error < tolerance,
                    "line {} coefficient {k} at {bits} bits: {value}; reference {reference:?}",
                    line_index + 1
                );
                if bits == 384 {
                    let difference = value.clone() - &previous[count][k];
                    let convergence =
                        difference.re.to_f64().hypot(difference.im.to_f64()) * normalization;
                    assert!(
                        convergence.is_finite() && convergence < 1e-32,
                        "line {} coefficient {k} failed precision convergence: {convergence:e}",
                        line_index + 1
                    );
                }
            }
            if bits == 256 {
                previous.push(result);
            }
            count += 1;
        }
        assert_eq!(count, 99, "all 49 B0 and 50 dB0 points must be tested");
    }
}

#[test]
fn derivative_negative_real_mass_limits_match_original_ratio_lips() {
    let parameters = [
        parse!("negative_limit_s"),
        parse!("negative_limit_m0"),
        parse!("negative_limit_m1"),
        parse!("negative_limit_mu"),
    ];
    let formula = db0(
        &parameters[0],
        &parameters[1],
        &parameters[2],
        &parameters[3],
    );
    let mut evaluator =
        exact(&formula, &parameters).map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let pi = std::f64::consts::PI;
    let cases = [
        (
            0.,
            -1.,
            2.,
            Complex::new(1. / 18. + 2. * 2_f64.ln() / 27., -2. * pi / 27.),
            0.,
        ),
        (
            0.,
            2.,
            -1.,
            Complex::new(1. / 18. + 2. * 2_f64.ln() / 27., -2. * pi / 27.),
            0.,
        ),
        (
            0.,
            -4.,
            1.,
            Complex::new(-3. / 50. - 4. * 4_f64.ln() / 125., -4. * pi / 125.),
            0.,
        ),
        (
            0.,
            1.,
            -4.,
            Complex::new(-3. / 50. - 4. * 4_f64.ln() / 125., -4. * pi / 125.),
            0.,
        ),
        (
            -2.,
            -2.,
            0.,
            Complex::new(0.5 - 2_f64.ln() / 4., pi / 4.),
            0.25,
        ),
        (
            -2.,
            0.,
            -2.,
            Complex::new(0.5 - 2_f64.ln() / 4., pi / 4.),
            0.25,
        ),
    ];
    for (s, m0, m1, expected, pole) in cases {
        for imaginary in [0., -0.] {
            let inputs = [
                Complex::new(s, 0.),
                Complex::new(m0, imaginary),
                Complex::new(m1, imaginary),
                Complex::new(1., 0.),
            ];
            let mut actual = [Complex::new(0., 0.); 3];
            evaluator.evaluate(&inputs, &mut actual);
            let difference = actual[0] - expected;
            assert!(
                difference.re.hypot(difference.im) < 1e-12,
                "s={s}, m0={m0}, m1={m1}: {actual:?}"
            );
            assert_eq!(actual[1], Complex::new(pole, 0.));
            assert_eq!(actual[2], Complex::new(0., 0.));
        }
    }
}

#[test]
fn derivative_normal_threshold_prescription_is_mass_symmetric() {
    let parameters = [
        parse!("ordered_threshold_s"),
        parse!("ordered_threshold_m0"),
        parse!("ordered_threshold_m1"),
        parse!("ordered_threshold_mu"),
    ];
    let formula = db0(
        &parameters[0],
        &parameters[1],
        &parameters[2],
        &parameters[3],
    );
    let mut evaluator =
        exact(&formula, &parameters).map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let pi = std::f64::consts::PI;
    for (s, m0, m1, expected) in [
        (
            9.,
            1.,
            4.,
            Complex::new(2_f64.ln() / 27. - 2. / 9., -pi / 27.),
        ),
        (
            -9.,
            -1.,
            -4.,
            Complex::new(2. / 9. - 2_f64.ln() / 27., -pi / 27.),
        ),
        (
            16.,
            1.,
            9.,
            Complex::new(3_f64.ln() / 32. - 1. / 8., -pi / 32.),
        ),
        (
            -16.,
            -1.,
            -9.,
            Complex::new(1. / 8. - 3_f64.ln() / 32., -pi / 32.),
        ),
        (-1., -1., -4., Complex::new(2. - 3. * 2_f64.ln(), 0.)),
    ] {
        for (first, second) in [(m0, m1), (m1, m0)] {
            let mut actual = [Complex::new(0., 0.); 3];
            evaluator.evaluate(
                &[s, first, second, 1.].map(|value| Complex::new(value, 0.)),
                &mut actual,
            );
            let difference = actual[0] - expected;
            assert!(
                difference.re.hypot(difference.im) < 1e-12,
                "s={s}, m0={first}, m1={second}: {actual:?}, expected {expected:?}"
            );
            assert_eq!(actual[1], Complex::new(0., 0.));
            assert_eq!(actual[2], Complex::new(0., 0.));
        }
    }
    let complex_masses = [Complex::new(3., -4.), Complex::new(0., -2.)];
    let expected = Complex::new(0.01793720660451692, 0.04896093131577151);
    for masses in [complex_masses, [complex_masses[1], complex_masses[0]]] {
        let mut actual = [Complex::new(0., 0.); 3];
        evaluator.evaluate(
            &[
                Complex::new(1., 0.),
                masses[0],
                masses[1],
                Complex::new(1., 0.),
            ],
            &mut actual,
        );
        let difference = actual[0] - expected;
        assert!(
            difference.re.hypot(difference.im) < 1e-12,
            "complex repeated roots {masses:?}: {actual:?}"
        );
        assert_eq!(actual[1], Complex::new(0., 0.));
        assert_eq!(actual[2], Complex::new(0., 0.));
    }
}
