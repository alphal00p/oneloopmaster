//! Fixed arbitrary precision: exact decimal inputs and independently known limits.
use oneloop::{A0, B0, C0, D0, PrecisionEvaluator, ScalarIntegral, dB0};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::{
        backend::float::Constant,
        float::{Complex, Float, Real, SingleFloat},
    },
    evaluate::ExpressionEvaluator,
};

type C = Complex<Float>;

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

fn decimal(value: &str, bits: u32) -> Float {
    Float::parse(value, Some(bits)).unwrap()
}

fn complex(re: &str, im: &str, bits: u32) -> C {
    Complex::new(decimal(re, bits), decimal(im, bits))
}

fn real(value: i32, bits: u32) -> C {
    Complex::new(Float::with_val(bits, value), Float::new(bits))
}

fn zeros(bits: u32) -> [C; 3] {
    std::array::from_fn(|_| real(0, bits))
}

fn bare_hook(family: ScalarIntegral, bits: u32) -> ExpressionEvaluator<C> {
    let symbol = match family {
        ScalarIntegral::A0 => A0(),
        ScalarIntegral::B0 => B0(),
        ScalarIntegral::DB0 => dB0(),
        ScalarIntegral::C0 => C0(),
        ScalarIntegral::D0 => D0(),
    };
    let parameters = (0..family.arity())
        .map(|index| {
            symbolica::symbol!(format!("precision_api_{}_{index}", family.name())).to_atom()
        })
        .collect::<Vec<_>>();
    let calls = [0, -1, -2].map(|tag| {
        let mut arguments = vec![Atom::num(tag)];
        arguments.extend(parameters.iter().cloned());
        symbol.call(&arguments)
    });
    // Deliberately no FunctionMap: this exercises the public Float hook.
    let exact = Atom::evaluator_multiple(
        &calls.iter().map(Atom::as_view).collect::<Vec<_>>(),
        &parameters,
    )
    .direct_translation(true)
    .build()
    .unwrap();
    let converter = Float::new(bits);
    exact.map_coeff_with_prec(
        &|value| {
            Complex::new(
                converter.from_rational(&value.re),
                converter.from_rational(&value.im),
            )
        },
        bits,
    )
}

fn absolute(value: Float) -> Float {
    if value.is_negative() { -value } else { value }
}

fn assert_close(actual: &C, expected: &C, digits: u32, label: &str) {
    assert!(
        actual.re.is_finite() && actual.im.is_finite(),
        "{label}: nonfinite output {actual}"
    );
    let bits = expected.re.prec().max(expected.im.prec());
    let tolerance = |component: &Float| {
        let scale = if component.is_zero() {
            Float::with_val(bits, 1)
        } else {
            absolute(component.clone())
        };
        decimal(&format!("1e-{digits}"), bits) * scale
    };
    let re_tolerance = tolerance(&expected.re);
    let im_tolerance = tolerance(&expected.im);
    let re_error = absolute(actual.re.clone() - expected.re.clone());
    let im_error = absolute(actual.im.clone() - expected.im.clone());
    assert!(
        re_error < re_tolerance && im_error < im_tolerance,
        "{label}: real error={re_error:e}, imaginary error={im_error:e}, real tolerance={re_tolerance:e}, imaginary tolerance={im_tolerance:e}"
    );
}

fn elementary_control(family: ScalarIntegral, bits: u32) -> (Vec<C>, [C; 3]) {
    let n = |value| Float::with_val(bits, value);
    let z = || Float::new(bits);
    let mut expected = zeros(bits);
    let arguments = match family {
        ScalarIntegral::A0 => {
            // The physical lower-lip logarithm fixes the negative-mass sheet.
            expected[0] = Complex::new(
                n(-2) * (n(1) - n(2).log()),
                n(-2) * Float::with_val(bits, Constant::Pi),
            );
            expected[1] = real(-2, bits);
            vec![-2, 1]
        }
        ScalarIntegral::B0 => {
            expected[0] = Complex::new(-(n(3) / n(5)).log(), z());
            expected[1] = real(1, bits);
            vec![0, 3, 3, 5]
        }
        ScalarIntegral::DB0 => {
            // A regular derivative pseudothreshold, not a numerical derivative.
            expected[0] = Complex::new(n(3) * n(2).log() - n(2), z());
            vec![1, 1, 4, 1]
        }
        ScalarIntegral::C0 => {
            expected[0] = Complex::new(-n(1) / n(2), z());
            vec![0, 0, 0, 1, 1, 1, 1]
        }
        ScalarIntegral::D0 => {
            expected[0] = Complex::new(n(1) / n(6), z());
            vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
        }
    };
    (
        arguments
            .into_iter()
            .map(|value| real(value, bits))
            .collect(),
        expected,
    )
}

#[test]
fn all_scalar_families_have_elementary_limits_at_16_32_and_1000_digits() {
    on_stack(|| {
        for digits in [16, 32, 1000] {
            for family in [
                ScalarIntegral::A0,
                ScalarIntegral::B0,
                ScalarIntegral::DB0,
                ScalarIntegral::C0,
                ScalarIntegral::D0,
            ] {
                eprintln!("elementary {} at {digits} decimal digits", family.name());
                let mut evaluator = PrecisionEvaluator::new(family, digits).unwrap();
                assert_eq!(evaluator.family(), family);
                assert_eq!(evaluator.precision(), digits);
                let bits = evaluator.binary_precision();
                let (arguments, expected) = elementary_control(family, bits + 256);
                let mut output = zeros(bits);
                evaluator.evaluate(&arguments, &mut output).unwrap();
                for tag in 0..3 {
                    assert_close(
                        &output[tag],
                        &expected[tag],
                        digits,
                        &format!("{} at {digits} digits, coefficient {tag}", family.name()),
                    );
                }
            }
        }
    });
}

#[test]
fn decimal_inputs_survive_large_integers_tiny_exponents_and_lower_lips() {
    on_stack(|| {
        for digits in [32, 1000] {
            let mut evaluator = PrecisionEvaluator::new(ScalarIntegral::A0, digits).unwrap();
            let bits = evaluator.binary_precision();
            let mut arguments = vec![
                complex("9007199254740993", "0", bits),
                real(1, bits),
                complex("1e-10000", "0", bits),
                real(1, bits),
                complex("-2", "-1e-10000", bits),
                real(1, bits),
            ];
            let mut output = vec![real(0, bits); 9];
            evaluator
                .evaluate_batch(&arguments, &mut output, 3)
                .unwrap();
            assert_eq!(output[1], arguments[0]);
            assert_ne!(output[1].re, decimal("9007199254740992", bits));
            assert_eq!(output[4], arguments[2]);
            assert!(!output[4].re.is_zero());
            assert_eq!(output[7], arguments[4]);
            assert!(!output[7].im.is_zero());
            assert!(output[7].im.is_negative());

            // Independent A0 formula, with a genuinely negative imaginary part
            // selecting the lower principal logarithm even below f64's range.
            for row in 0..3 {
                let m = &arguments[2 * row];
                let expected = m * (real(1, bits) - m.log());
                assert_close(&output[3 * row], &expected, digits, "exact-decimal A0");
                if row == 1 {
                    assert_close(
                        &(output[3 * row].clone() / expected),
                        &real(1, bits),
                        digits,
                        "relative accuracy of tiny-exponent A0",
                    );
                }
                assert_eq!(output[3 * row + 2], real(0, bits));
            }

            // Reusing and cloning an evaluator must not keep a previous point.
            let mut cloned = evaluator.clone();
            arguments[0] = real(2, bits);
            let mut first = zeros(bits);
            let mut second = zeros(bits);
            evaluator.evaluate(&arguments[..2], &mut first).unwrap();
            cloned.evaluate(&arguments[..2], &mut second).unwrap();
            assert_eq!(first, second);
        }
    });
}

#[test]
fn shape_and_exact_domain_validation_are_transactional() {
    assert!(PrecisionEvaluator::new(ScalarIntegral::A0, 0).is_err());
    assert!(PrecisionEvaluator::new(ScalarIntegral::A0, u32::MAX).is_err());
    assert!(PrecisionEvaluator::with_binary_precision(ScalarIntegral::A0, 0).is_err());
    on_stack(|| {
        let mut evaluator =
            PrecisionEvaluator::with_binary_precision(ScalarIntegral::B0, 128).unwrap();
        assert_eq!(evaluator.binary_precision(), 128);
        assert_eq!(evaluator.precision(), 38);
        let bits = 256;
        let point = vec![real(0, bits), real(1, bits), real(1, bits), real(1, bits)];
        let sentinel = vec![real(42, bits); 6];
        let mut output = sentinel.clone();
        assert!(evaluator.evaluate(&point[..3], &mut output[..3]).is_err());
        assert!(evaluator.evaluate(&point, &mut output[..2]).is_err());
        assert!(evaluator.evaluate_batch(&[], &mut [], usize::MAX).is_err());
        evaluator.evaluate_batch(&[], &mut [], 0).unwrap();
        assert_eq!(output, sentinel);

        for (index, invalid) in [
            (0, complex("0", "1e-10000", bits)),
            (1, complex("1", "1e-10000", bits)),
            (2, complex("1", "1e-10000", bits)),
            (3, complex("1", "-1e-10000", bits)),
            (3, real(0, bits)),
            (3, real(-1, bits)),
            (1, complex("NaN", "0", bits)),
            (2, complex("1", "inf", bits)),
        ] {
            let mut batch = point.clone();
            let mut invalid_point = point.clone();
            invalid_point[index] = invalid;
            batch.extend(invalid_point);
            assert!(evaluator.evaluate_batch(&batch, &mut output, 2).is_err());
            assert_eq!(
                output, sentinel,
                "invalid second row changed earlier outputs"
            );
        }
        // Negative squared masses and a tiny positive real scale are valid.
        let valid = vec![
            real(-1, bits),
            complex("-2", "-1e-10000", bits),
            real(1, bits),
            complex("1e-10000", "-0", bits),
        ];
        evaluator.evaluate(&valid, &mut output[..3]).unwrap();
    });
}

#[test]
fn bare_master_float_callbacks_preserve_1000_decimal_digits() {
    on_stack(|| {
        for family in [
            ScalarIntegral::A0,
            ScalarIntegral::B0,
            ScalarIntegral::DB0,
            ScalarIntegral::C0,
            ScalarIntegral::D0,
        ] {
            let mut native = PrecisionEvaluator::new(family, 1000).unwrap();
            let bits = native.binary_precision();
            let mut hook = bare_hook(family, bits);
            let (mut arguments, _) = elementary_control(family, bits);
            if family == ScalarIntegral::A0 {
                arguments[0] = complex("-2", "-1e-10000", bits);
            }
            let mut actual = zeros(bits);
            let mut expected = zeros(bits);
            hook.evaluate(&arguments, &mut actual);
            native.evaluate(&arguments, &mut expected).unwrap();
            for tag in 0..3 {
                assert_close(&actual[tag], &expected[tag], 1000, "bare Float hook");
            }
            if family == ScalarIntegral::A0 {
                assert_eq!(actual[1].im, arguments[0].im);
            }
            if family == ScalarIntegral::C0 {
                // A second point exercises genuinely complex dilogarithms,
                // unlike the elementary vacuum control above.
                let arguments = [
                    real(-1, bits),
                    real(-2, bits),
                    real(-3, bits),
                    complex("1", "-0.1", bits),
                    complex("2", "-0.2", bits),
                    complex("3", "-0.3", bits),
                    real(1, bits),
                ];
                hook.evaluate(&arguments, &mut actual);
                native.evaluate(&arguments, &mut expected).unwrap();
                for tag in 0..3 {
                    assert_close(&actual[tag], &expected[tag], 1000, "complex C0 Float hook");
                }
            }
        }
    });
}

#[test]
fn complex_triangle_and_box_are_stable_from_1000_to_1120_digits() {
    on_stack(|| {
        for (family, momenta, masses) in [
            (
                ScalarIntegral::C0,
                vec!["-1", "-2", "-3"],
                vec![("1", "-0.1"), ("2", "-0.2"), ("3", "-0.3")],
            ),
            (
                ScalarIntegral::D0,
                vec!["-1", "-2", "-3", "-4", "-5", "-6"],
                vec![("1", "-0.1"), ("2", "-0.2"), ("3", "-0.3"), ("4", "-0.4")],
            ),
        ] {
            let mut outputs = Vec::new();
            for digits in [1000, 1120] {
                eprintln!("complex {} at {digits} decimal digits", family.name());
                let mut evaluator = PrecisionEvaluator::new(family, digits).unwrap();
                let bits = evaluator.binary_precision();
                let mut arguments = momenta
                    .iter()
                    .map(|p| complex(p, "0", bits))
                    .collect::<Vec<_>>();
                arguments.extend(masses.iter().map(|(re, im)| complex(re, im, bits)));
                arguments.push(real(1, bits));
                let mut output = zeros(bits);
                evaluator.evaluate(&arguments, &mut output).unwrap();
                outputs.push(output);
            }
            // This is a precision-stability check, not an independent oracle.
            // The elementary tests above provide independent analytic targets.
            for (actual, reference) in outputs[0].iter().zip(&outputs[1]) {
                assert_close(actual, reference, 1000, family.name());
            }
        }
    });
}
