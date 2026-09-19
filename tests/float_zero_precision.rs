//! Adding exact zero must not advertise recovered arbitrary-precision bits.
use symbolica::domains::float::Float;

fn same(value: &Float, expected: &Float) {
    assert_eq!(value.as_raw(), expected.as_raw());
    assert_eq!(
        value.prec(),
        expected.prec(),
        "zero inflated a nonzero value"
    );
    assert_eq!(value.is_sign_negative(), expected.is_sign_negative());
}

#[test]
fn zero_addition_and_subtraction_preserve_nonzero_precision() {
    for bits in [2, 53, 512, 3456] {
        for zero_bits in [2, 53, 512, 3456] {
            for text in ["1.5", "-1.5", "0.125", "-1e-1000"] {
                let value = Float::parse(text, Some(bits)).unwrap();
                let negative = -value.clone();
                for zero_sign in [0.0, -0.0] {
                    let zero = Float::with_val(zero_bits, zero_sign);
                    for result in [
                        value.clone() + &zero,
                        value.clone() + zero.clone(),
                        zero.clone() + &value,
                        zero.clone() + value.clone(),
                        value.clone() - &zero,
                        value.clone() - zero.clone(),
                    ] {
                        same(&result, &value);
                    }
                    same(&(zero.clone() - &value), &negative);
                    same(&(zero.clone() - value.clone()), &negative);

                    let mut result = value.clone();
                    result += &zero;
                    same(&result, &value);
                    result += zero.clone();
                    same(&result, &value);
                    result -= &zero;
                    same(&result, &value);
                    result -= zero.clone();
                    same(&result, &value);
                    for borrowed in [false, true] {
                        let mut result = zero.clone();
                        if borrowed {
                            result += &value;
                        } else {
                            result += value.clone();
                        }
                        same(&result, &value);
                        let mut result = zero.clone();
                        if borrowed {
                            result -= &value;
                        } else {
                            result -= value.clone();
                        }
                        same(&result, &negative);
                    }
                }
            }
        }
    }
}

#[test]
fn signed_zero_and_nonfinite_addition_keep_numerical_behavior() {
    for left in [0.0_f64, -0.0] {
        for right in [0.0_f64, -0.0] {
            let a = Float::with_val(2, left);
            let b = Float::with_val(512, right);
            let sum = a.clone() + &b;
            let difference = a.clone() - &b;
            assert_eq!(sum.is_sign_negative(), (left + right).is_sign_negative());
            assert_eq!(
                difference.is_sign_negative(),
                (left - right).is_sign_negative()
            );
            assert_eq!(sum.prec(), 512);
            assert_eq!(difference.prec(), 512);
            let mut sum = a.clone();
            sum += &b;
            let mut difference = a;
            difference -= &b;
            assert_eq!(sum.is_sign_negative(), (left + right).is_sign_negative());
            assert_eq!(
                difference.is_sign_negative(),
                (left - right).is_sign_negative()
            );
        }
    }
    for value in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        let a = Float::with_val(2, value);
        let zero = Float::new(512);
        for result in [a.clone() + &zero, zero.clone() + &a, a.clone() - &zero] {
            assert!(!result.is_finite());
            assert_eq!(result.as_raw().is_nan(), value.is_nan());
            if !value.is_nan() {
                assert_eq!(result.is_sign_negative(), value.is_sign_negative());
            }
        }
    }
}
