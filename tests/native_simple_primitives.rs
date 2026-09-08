//! State-free analytic controls for the exactly selected simple sectors.
#[allow(dead_code)]
#[path = "../src/native/primitives.rs"]
mod primitives;
#[path = "../src/native/simple.rs"]
mod simple;

use primitives::{C, NativeFloat};
use symbolica::domains::float::{Float, FloatLike, Real};

const REFERENCE_BITS: u32 = 4096;

fn real(value: i64) -> Float {
    Float::with_val(REFERENCE_BITS, value)
}

// Independent reference: take the logarithm of the ratio directly, at higher
// precision. Do not call simple::momentum_log or the native complex primitives.
fn logarithm(s: &Float, mu: &Float) -> C<Float> {
    C::new(
        (s.norm() / mu).log(),
        if *s > real(0) { -real(1).pi() } else { real(0) },
    )
}

fn triangle_reference(momenta: &[Float; 3], mu: &Float) -> [C<Float>; 3] {
    let scales = momenta
        .iter()
        .filter(|s| !s.is_fully_zero())
        .collect::<Vec<_>>();
    let zero = || C::new(real(0), real(0));
    if scales.is_empty() {
        return [zero(), zero(), zero()];
    }
    let s = scales[0];
    let ls = logarithm(s, mu);
    if scales.len() == 1 {
        // Public normalization cancels the raw -pi²/(12s) exactly.
        return [
            C::new(
                (ls.re.clone() * &ls.re - ls.im.clone() * &ls.im) / (real(2) * s),
                ls.re.clone() * &ls.im / s,
            ),
            C::new(-ls.re / s, -ls.im / s),
            C::new(s.inv(), real(0)),
        ];
    }
    assert_eq!(scales.len(), 2);
    let t = scales[1];
    if s == t {
        return [
            C::new(ls.re / s, ls.im / s),
            C::new(-s.inv(), real(0)),
            zero(),
        ];
    }
    let lt = logarithm(t, mu);
    let difference = s.clone() - t;
    // A difference of squared logarithms, independently of the production
    // factorization into the pole times the sum of logarithms.
    [
        C::new(
            (ls.re.clone() * &ls.re - ls.im.clone() * &ls.im - lt.re.clone() * &lt.re
                + lt.im.clone() * &lt.im)
                / (real(2) * &difference),
            (ls.re.clone() * &ls.im - lt.re.clone() * &lt.im) / &difference,
        ),
        C::new((lt.re - ls.re) / &difference, (lt.im - ls.im) / difference),
        zero(),
    ]
}

fn box_reference(s: &Float, t: &Float, mu: &Float) -> [C<Float>; 3] {
    let ls = logarithm(s, mu);
    let lt = logarithm(t, mu);
    let difference = ls.clone() - &lt;
    let factor = (s.clone() * t).inv();
    let pi_squared = real(1).pi().pow(2);
    // Keep the raw squared-log formula and add the public normalization
    // separately, rather than duplicating the simplified 2 Ls Lt - pi².
    let raw = ls.clone() * &ls + lt.clone() * &lt
        - difference.clone() * difference
        - C::new(real(4) * &pi_squared / real(3), real(0));
    let finite = raw * C::new(factor.clone(), real(0))
        + C::new(real(4) * &factor * pi_squared / real(12), real(0));
    let pole = (ls + lt) * C::new(-real(2) * &factor, real(0));
    [finite, pole, C::new(real(4) * factor, real(0))]
}

fn check<T: NativeFloat>(
    actual: &[C<T>; 3],
    expected: &[C<Float>; 3],
    convert: &impl Fn(&T) -> Float,
    tolerance: &Float,
) {
    for (tag, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        for (component, (actual, expected)) in
            [(&actual.re, &expected.re), (&actual.im, &expected.im)]
                .into_iter()
                .enumerate()
        {
            assert!(
                actual.is_finite(),
                "tag {tag} component {component} is nonfinite"
            );
            let error = (convert(actual) - expected).norm();
            let scale = if expected.is_fully_zero() {
                real(1)
            } else {
                expected.norm()
            };
            assert!(
                error <= scale * tolerance,
                "tag {tag} component {component}: actual={actual:e}, expected={expected:e}, error={error:e}"
            );
        }
    }
}

fn sectors<T: NativeFloat>(bits: u32, convert: impl Fn(&T) -> Float, tolerance: &str) {
    let prototype = T::zero_at(bits);
    let constants = simple::Constants::new(&prototype);
    let tolerance = Float::parse(tolerance, Some(REFERENCE_BITS)).unwrap();
    let z = || C::new(prototype.zero(), prototype.zero());
    let permutations = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for momenta in [
        [0, 0, 0],
        [0, 0, -1],
        [0, 0, 1],
        [0, 0, -2],
        [0, 0, 2],
        [0, -1, -2],
        [0, 1, 2],
        [0, -1, 2],
        [0, 1, -2],
        [0, -2, -2],
        [0, 2, 2],
    ] {
        for mu in [1, 2, 4] {
            let expected = triangle_reference(&momenta.map(real), &real(mu));
            for order in permutations {
                for negative_zero in [false, true] {
                    let mut input = core::array::from_fn::<_, 7, _>(|_| z());
                    for index in 0..3 {
                        input[index].re = prototype.from_i64(momenta[order[index]]);
                    }
                    input[6].re = prototype.from_i64(mu);
                    if negative_zero {
                        for value in &mut input {
                            value.im = -value.im.clone();
                            if value.re.is_fully_zero() {
                                value.re = -value.re.clone();
                            }
                        }
                    }
                    let actual = simple::massless_triangle(&input, &constants).unwrap();
                    check(&actual, &expected, &convert, &tolerance);
                }
            }
        }
    }
    for (s, t) in [(-1, -2), (1, 2), (-1, 2), (1, -2), (-2, -2), (2, 2)] {
        for mu in [1, 2, 4] {
            let expected = box_reference(&real(s), &real(t), &real(mu));
            for channels in [[s, t], [t, s]] {
                for negative_zero in [false, true] {
                    let mut input = core::array::from_fn::<_, 11, _>(|_| z());
                    input[4].re = prototype.from_i64(channels[0]);
                    input[5].re = prototype.from_i64(channels[1]);
                    input[10].re = prototype.from_i64(mu);
                    if negative_zero {
                        for value in &mut input {
                            value.im = -value.im.clone();
                            if value.re.is_fully_zero() {
                                value.re = -value.re.clone();
                            }
                        }
                    }
                    let actual = simple::massless_onshell_box(&input, &constants).unwrap();
                    check(&actual, &expected, &convert, &tolerance);
                }
            }
        }
    }
}

#[test]
fn simple_massless_coefficients_binary64() {
    sectors::<f64>(53, |value| Float::with_val(REFERENCE_BITS, *value), "2e-14");
}

#[test]
fn simple_massless_coefficients_thousand_decimal_digits() {
    sectors::<Float>(3456, |value| value.at_precision(REFERENCE_BITS), "1e-1000");
}

fn exact_selection<T: NativeFloat>(bits: u32, tiny: T) {
    let prototype = T::zero_at(bits);
    let constants = simple::Constants::new(&prototype);
    let z = || C::new(prototype.zero(), prototype.zero());
    let mut triangle = core::array::from_fn::<_, 7, _>(|_| z());
    triangle[0].re = prototype.one();
    triangle[1].re = prototype.from_i64(-2);
    triangle[6].re = prototype.one();
    assert!(simple::massless_triangle(&triangle, &constants).is_some());
    triangle[2].re = tiny.clone();
    assert!(simple::massless_triangle(&triangle, &constants).is_none());
    triangle[2] = z();
    for mass in 3..6 {
        triangle[mass].im = -tiny.clone();
        assert!(simple::massless_triangle(&triangle, &constants).is_none());
        triangle[mass] = z();
    }
    let mut box_input = core::array::from_fn::<_, 11, _>(|_| z());
    box_input[4].re = prototype.one();
    box_input[5].re = prototype.from_i64(-2);
    box_input[10].re = prototype.one();
    assert!(simple::massless_onshell_box(&box_input, &constants).is_some());
    for external in 0..4 {
        box_input[external].re = tiny.clone();
        assert!(simple::massless_onshell_box(&box_input, &constants).is_none());
        box_input[external] = z();
    }
    for mass in 6..10 {
        box_input[mass].im = -tiny.clone();
        assert!(simple::massless_onshell_box(&box_input, &constants).is_none());
        box_input[mass] = z();
    }
    for channel in [4, 5] {
        let saved = box_input[channel].clone();
        box_input[channel] = z();
        assert!(simple::massless_onshell_box(&box_input, &constants).is_none());
        box_input[channel] = saved;
    }
}

#[test]
fn simple_sector_selection_uses_exact_zeros() {
    exact_selection::<f64>(53, f64::from_bits(1));
    exact_selection::<Float>(3456, Float::parse("1e-1000", Some(3456)).unwrap());
}
