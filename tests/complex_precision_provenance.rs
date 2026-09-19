//! Independent regressions for exact constants and mixed-precision phases.
use symbolica::domains::float::{Complex, Float, FloatLike, Real};

#[test]
fn a_small_uncertain_component_does_not_spoil_a_known_large_phase() {
    for (bits, tiny, tolerance) in [(512, "1e-170", "1e-140"), (3456, "1e-1100", "1e-1000")] {
        let n = |text| Float::parse(text, Some(bits)).unwrap();
        let tiny = Float::parse(tiny, Some(2)).unwrap();
        let pi = n("0").pi();
        for imaginary_negative in [false, true] {
            let sign = if imaginary_negative { -1 } else { 1 };
            let y = tiny.clone() * sign;
            let actual = y.atan2(&n("-1"));
            let expected = pi.clone() * sign;
            assert!((actual.clone() - expected).norm() < n(tolerance));
            assert!(
                actual.prec() >= bits - 8,
                "lost phase bits: {}",
                actual.prec()
            );
            for real_negative in [false, true] {
                let x = if real_negative {
                    -tiny.clone()
                } else {
                    tiny.clone()
                };
                let actual = (n("1.25") * sign).atan2(&x);
                let expected = pi.clone() * sign / 2;
                assert!((actual.clone() - expected).norm() < n(tolerance));
                assert!(
                    actual.prec() >= bits - 8,
                    "lost phase bits: {}",
                    actual.prec()
                );
            }
        }
        // A small angle itself is known to only two relative bits. Do not
        // claim the large component can recover the missing input digits.
        let actual = tiny.atan2(&n("1"));
        assert!(actual.prec() <= 3);
    }
}

#[test]
fn exact_integer_power_identity_uses_the_strong_component() {
    for (bits, tiny, tolerance) in [(512, "1e-170", "1e-140"), (3456, "1e-1100", "1e-1000")] {
        let n = |text| Float::parse(text, Some(bits)).unwrap();
        let tiny = Float::parse(tiny, Some(2)).unwrap();
        for real_negative in [false, true] {
            for imaginary_negative in [false, true] {
                let z = Complex::new(
                    if real_negative {
                        -tiny.clone()
                    } else {
                        tiny.clone()
                    },
                    if imaginary_negative {
                        -n("1.25")
                    } else {
                        n("1.25")
                    },
                );
                let actual = z.pow(2);
                assert!((actual.re.clone() - n("-1.5625")).norm() < n(tolerance));
                assert!(actual.re.prec() >= bits - 8);
                assert_eq!(z.pow(0).re, n("1"));
            }
        }
    }
}

#[test]
fn integer_powers_preserve_real_axis_signed_lips() {
    for bits in [512, 3456] {
        for real in [-2, 2] {
            for imaginary_negative in [false, true] {
                let z = Complex::new(
                    Float::with_val(bits, real),
                    Float::with_val(bits, if imaginary_negative { -0.0 } else { 0.0 }),
                );
                for exponent in 1..=5 {
                    let actual = z.pow(exponent);
                    assert!(actual.im.is_fully_zero());
                    let expected_negative = imaginary_negative ^ (real < 0 && exponent % 2 == 0);
                    assert_eq!(
                        actual.im.is_sign_negative(),
                        expected_negative,
                        "base={real}, exponent={exponent}, imaginary_negative={imaginary_negative}"
                    );
                    assert_eq!(
                        actual.re,
                        Float::with_val(bits, i64::from(real).pow(exponent as u32))
                    );
                }
            }
        }
    }
}

#[test]
fn phase_preserves_axis_and_signed_zero_conventions() {
    for (real_bits, imaginary_bits) in [(2, 512), (512, 2), (512, 512)] {
        for real in ["-1", "-0", "0", "1"] {
            for imaginary in ["-1", "-0", "0", "1"] {
                let x = Float::parse(real, Some(real_bits)).unwrap();
                let y = Float::parse(imaginary, Some(imaginary_bits)).unwrap();
                if !x.is_fully_zero() && !y.is_fully_zero() {
                    continue;
                }
                let mut expected = y.clone();
                expected.set_prec(512);
                let expected = expected.as_raw().clone().atan2(x.as_raw());
                let actual = y.atan2(&x);
                let error = (actual.clone() - Float::from(expected)).norm();
                assert!(
                    error < Float::parse("1e-140", Some(512)).unwrap(),
                    "atan2({imaginary}, {real})"
                );
                if actual.is_fully_zero() {
                    assert_eq!(actual.is_sign_negative(), y.is_sign_negative());
                }
            }
        }
    }
}

#[test]
fn half_powers_use_the_precision_safe_square_root() {
    for (bits, tiny, tolerance) in [(512, "1e-170", "1e-140"), (3456, "1e-1100", "1e-1000")] {
        let n = |text| Float::parse(text, Some(bits)).unwrap();
        let tiny = Float::parse(tiny, Some(2)).unwrap();
        for real_negative in [false, true] {
            for imaginary_negative in [false, true] {
                for tiny_real in [false, true] {
                    let re = if tiny_real { tiny.clone() } else { n("2") };
                    let im = if tiny_real { n("1.25") } else { tiny.clone() };
                    let z = Complex::new(
                        if real_negative { -re } else { re },
                        if imaginary_negative { -im } else { im },
                    );
                    let mut reference = z.clone();
                    reference.re.set_prec(bits + 128);
                    reference.im.set_prec(bits + 128);
                    let reference = reference.sqrt();
                    for exponent in ["-1.5", "-0.5", "0.5", "1.5"] {
                        let actual = z.powf(&Complex::new(n(exponent), n("0")));
                        let expected = if exponent.ends_with("1.5") {
                            let mut strong_z = z.clone();
                            strong_z.re.set_prec(bits + 128);
                            strong_z.im.set_prec(bits + 128);
                            reference.clone() * strong_z
                        } else {
                            reference.clone()
                        };
                        let expected = if exponent.starts_with('-') {
                            expected.inv()
                        } else {
                            expected
                        };
                        assert!(
                            (actual.re.clone() - expected.re).norm() < n(tolerance),
                            "exponent={exponent}"
                        );
                        assert!(
                            (actual.im.clone() - expected.im).norm() < n(tolerance),
                            "exponent={exponent}"
                        );
                    }
                }
            }
        }
        for re in ["-2", "2"] {
            for im in ["-0", "0"] {
                let z = Complex::new(n(re), n(im));
                for exponent in ["-1.5", "-0.5", "0.5", "1.5"] {
                    let actual = z.powf(&Complex::new(n(exponent), n("0")));
                    let expected = if exponent.ends_with("1.5") {
                        z.clone() * z.sqrt()
                    } else {
                        z.sqrt()
                    };
                    let expected = if exponent.starts_with('-') {
                        expected.inv()
                    } else {
                        expected
                    };
                    // Equivalent reciprocal-first products may differ by a
                    // few final rounding bits; the branch lips must agree.
                    assert!((actual.re.clone() - expected.re.clone()).norm() < n(tolerance));
                    assert!((actual.im.clone() - expected.im.clone()).norm() < n(tolerance));
                    assert_eq!(actual.re.is_sign_negative(), expected.re.is_sign_negative());
                    assert_eq!(actual.im.is_sign_negative(), expected.im.is_sign_negative());
                }
            }
        }
    }
}
