//! Original-Fortran regressions and independent analytic triangle limits.
use oneloop::OneLoopExpressions;
use symbolica::prelude::*;

#[test]
fn accepted_triangle_limits_and_real_mass_lips() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(check_limits)
        .unwrap()
        .join()
        .unwrap();
}

fn check_limits() {
    let parameters = [
        "tc_p1", "tc_p2", "tc_p3", "tc_m1", "tc_m2", "tc_m3", "tc_mu",
    ]
    .map(|name| symbol!(name).to_atom());
    let context = OneLoopExpressions::new();
    let series = context.c0(
        [&parameters[0], &parameters[1], &parameters[2]],
        [&parameters[3], &parameters[4], &parameters[5]],
        &parameters[6],
    );
    let exact = context
        .evaluator(series.coefficients(), &parameters)
        .unwrap();
    let zero = Complex::new(0., 0.);
    let mut cases = vec![
        (
            [-2., 0., -1.],
            [zero, zero, Complex::new(1., 0.)],
            [Complex::new(-1.0626935403834068, 0.), zero, zero],
        ),
        (
            [-6., -1., -2.],
            [zero, zero, Complex::new(1., 0.)],
            [Complex::new(-0.5320783565121019, 0.), zero, zero],
        ),
        (
            [-3., 0., -2.],
            [zero, zero, Complex::new(1., 0.)],
            [Complex::new(-0.8116621803719681, 0.), zero, zero],
        ),
        (
            [1., 0., 2.],
            [zero, zero, Complex::new(1., 0.)],
            [
                Complex::new(0.8224670334243456, -2.177586090303606),
                zero,
                zero,
            ],
        ),
        (
            [2., 1., 1.],
            [zero, zero, Complex::new(0., -1.)],
            [
                Complex::new(0.915965594177171, -0.6168502750681577),
                zero,
                zero,
            ],
        ),
        (
            [-2., -1., -1.],
            [zero, zero, Complex::new(0., -1.)],
            [
                Complex::new(-0.915965594177292, -0.6168502750679563),
                zero,
                zero,
            ],
        ),
        (
            [2., 0., 1.],
            [zero, zero, Complex::new(-1., 0.)],
            [Complex::new(1.0626935403834068, 0.), zero, zero],
        ),
        (
            [-1., 0., -2.],
            [zero, zero, Complex::new(-1., 0.)],
            [
                Complex::new(-0.8224670334243456, -2.1775860903035778),
                zero,
                zero,
            ],
        ),
        (
            [-1., -3., 2.],
            [zero, Complex::new(-1., 0.), Complex::new(2., 0.)],
            [
                Complex::new(-0.5329515763938921, 0.06829261308509249),
                Complex::new(0.23397032752589253, 0.320637457540466),
                zero,
            ],
        ),
        (
            [-1., 0., 4.],
            [zero, Complex::new(-1., 0.), Complex::new(4., 0.)],
            [
                Complex::new(-0.5895708228381082, 0.),
                Complex::new(0.13862943611198908, 0.3141592653589793),
                zero,
            ],
        ),
    ];
    for (p, finite) in [
        ([-1., 3., 4.], -0.13850921108485956),
        ([0., 6., 4.], -0.15080867256915886),
        ([1., 9., 4.], -0.16612872594673547),
    ] {
        let m = [
            Complex::new(1., 0.),
            Complex::new(4., 0.),
            Complex::new(9., 0.),
        ];
        for shift in 0..3 {
            cases.push((
                core::array::from_fn(|i| p[(i + shift) % 3]),
                core::array::from_fn(|i| m[(i + shift) % 3]),
                [Complex::new(finite, 0.), zero, zero],
            ));
        }
    }
    // This Euclidean point is finite even though the original Fortran reports
    // an error. Its projective denominator has A=B=2 and A+B=A*B, giving
    // C0=-(log(A)+(A-1)*log(B))/A=-log(2), with both Laurent poles zero.
    let analytic_index = cases.len();
    cases.push((
        [-4., -1., -1.],
        [zero, zero, Complex::new(1., 0.)],
        [Complex::new(-2_f64.ln(), 0.), zero, zero],
    ));
    let mut double = exact
        .clone()
        .map_coeff(&|value| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut failures = Vec::new();
    for (index, (p, m, expected)) in cases.iter().enumerate() {
        let input = [
            Complex::new(p[0], 0.),
            Complex::new(p[1], 0.),
            Complex::new(p[2], 0.),
            m[0],
            m[1],
            m[2],
            Complex::new(1., 0.),
        ];
        let mut actual = [zero; 3];
        double.evaluate(&input, &mut actual);
        for coefficient in 0..3 {
            let delta = actual[coefficient] - expected[coefficient];
            let error = delta.re.hypot(delta.im);
            if !error.is_finite() || error >= 5e-10 {
                failures.push(format!(
                    "case {index} coefficient {coefficient} at f64: {actual:?}, expected {expected:?}"
                ));
            }
        }
    }
    // The same exact expressions must retain their lips beyond double precision.
    for bits in [256, 384] {
        let mut precise = exact.clone().map_coeff_with_prec(
            &|value| {
                Complex::new(
                    value.re.to_multi_prec_float(bits),
                    value.im.to_multi_prec_float(bits),
                )
            },
            bits,
        );
        for (index, (p, m, expected)) in cases.iter().enumerate() {
            let input = [
                Complex::new(p[0], 0.),
                Complex::new(p[1], 0.),
                Complex::new(p[2], 0.),
                m[0],
                m[1],
                m[2],
                Complex::new(1., 0.),
            ]
            .map(|value| {
                Complex::new(
                    Float::with_val(bits, value.re),
                    Float::with_val(bits, value.im),
                )
            });
            let mut actual = core::array::from_fn::<_, 3, _>(|_| {
                Complex::new(Float::new(bits), Float::new(bits))
            });
            precise.evaluate(&input, &mut actual);
            if index == analytic_index {
                let expected = -Float::with_val(bits, 2).log();
                let error = (actual[0].re.clone() - expected)
                    .to_f64()
                    .hypot(actual[0].im.to_f64());
                if !error.is_finite() || error >= 1e-60 {
                    failures.push(format!(
                        "independent -log(2) limit at {bits} bits: {:?}, error {error}",
                        actual[0]
                    ));
                }
            }
            for coefficient in 0..3 {
                let error = (actual[coefficient].re.to_f64() - expected[coefficient].re)
                    .hypot(actual[coefficient].im.to_f64() - expected[coefficient].im);
                if !error.is_finite() || error >= 5e-10 {
                    failures.push(format!(
                        "case {index} coefficient {coefficient} at {bits} bits: {actual:?}, expected {expected:?}"
                    ));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} triangle failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
