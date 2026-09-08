//! Direct-native contract checks. Deliberately never initialize Symbolica State.
use oneloop::{NativeEvaluator, NativeFloat, ScalarIntegral};
use symbolica::domains::float::{Complex, DoubleFloat, Float, FloatLike};

type C = Complex<f64>;

fn valid_point(family: ScalarIntegral) -> Vec<C> {
    let momenta = match family {
        ScalarIntegral::A0 => 0,
        ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
        ScalarIntegral::C0 => 3,
        ScalarIntegral::D0 => 6,
    };
    let mut point = vec![C::new(1.0, 0.0); family.arity()];
    for momentum in &mut point[..momenta] {
        momentum.re = -1.0;
    }
    point
}

#[test]
fn native_shape_domain_and_batch_errors_are_atomic() {
    let sentinel = [C::new(901.25, -704.5); 3];
    for family in [
        ScalarIntegral::A0,
        ScalarIntegral::B0,
        ScalarIntegral::DB0,
        ScalarIntegral::C0,
        ScalarIntegral::D0,
    ] {
        let mut evaluator = NativeEvaluator::<f64>::new(family).unwrap();
        assert_eq!(evaluator.family(), family);
        assert_eq!(evaluator.binary_precision(), 53);
        let point = valid_point(family);
        let mut output = sentinel;
        assert!(
            evaluator
                .evaluate(&point[..point.len() - 1], &mut output)
                .is_err()
        );
        assert_eq!(output, sentinel);
        assert!(evaluator.evaluate(&point, &mut output[..2]).is_err());
        assert_eq!(output, sentinel);
        assert!(
            evaluator
                .evaluate_batch(&[], &mut output, usize::MAX)
                .is_err()
        );
        assert_eq!(output, sentinel);
        evaluator.evaluate_batch(&[], &mut [], 0).unwrap();
        assert!(evaluator.evaluate_batch(&[], &mut output, 0).is_err());
        assert_eq!(output, sentinel);

        let momenta = family.arity()
            - match family {
                ScalarIntegral::A0 => 2,
                ScalarIntegral::B0 | ScalarIntegral::DB0 => 3,
                ScalarIntegral::C0 => 4,
                ScalarIntegral::D0 => 5,
            };
        let mut invalid = Vec::new();
        for component in [
            C::new(f64::NAN, 0.0),
            C::new(1.0, f64::INFINITY),
            C::new(1.0, f64::from_bits(1)),
        ] {
            let mut row = point.clone();
            row[momenta] = component;
            invalid.push(row);
        }
        if momenta != 0 {
            let mut row = point.clone();
            row[0].im = -f64::from_bits(1);
            invalid.push(row);
        }
        for mu in [
            C::new(0.0, 0.0),
            C::new(-0.0, -0.0),
            C::new(-1.0, 0.0),
            C::new(1.0, f64::from_bits(1)),
        ] {
            let mut row = point.clone();
            *row.last_mut().unwrap() = mu;
            invalid.push(row);
        }
        for row in invalid {
            let mut batch = point.clone();
            batch.extend(row);
            let mut output = sentinel.repeat(2);
            assert!(evaluator.evaluate_batch(&batch, &mut output, 2).is_err());
            assert_eq!(
                output,
                sentinel.repeat(2),
                "{family:?}: invalid second row changed first output"
            );
        }
    }
    assert!(!oneloop::is_initialized());
}

#[test]
fn native_precision_and_clone_contract() {
    assert!(NativeEvaluator::<f64>::with_binary_precision(ScalarIntegral::A0, 0).is_err());
    assert!(NativeEvaluator::<f64>::with_binary_precision(ScalarIntegral::A0, 106).is_err());
    assert!(NativeEvaluator::<DoubleFloat>::with_binary_precision(ScalarIntegral::A0, 53).is_err());
    assert!(NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::A0, 0).is_err());

    let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::A0).unwrap();
    let mut clone = evaluator.clone();
    let mut scaled = [C::new(0.0, 0.0); 3];
    evaluator
        .evaluate(&[C::new(2.0, 0.0), C::new(4.0, 0.0)], &mut scaled)
        .unwrap();
    assert!((scaled[0].re - 2.0 * (1.0 + 2.0_f64.ln())).abs() < 2e-15);
    assert_eq!(scaled[1], C::new(2.0, 0.0));
    assert_eq!(scaled[2], C::new(0.0, 0.0));
    drop(evaluator);
    let point = [C::new(2.0, -0.0), C::new(1.0, -0.0)];
    let mut batch = [C::new(f64::NAN, f64::NAN); 15];
    clone
        .evaluate_batch(&point.repeat(5), &mut batch, 5)
        .unwrap();
    for row in batch.chunks_exact(3) {
        assert!((row[0].re - 2.0 * (1.0 - 2.0_f64.ln())).abs() < 2e-15);
        assert_eq!(row[1], C::new(2.0, 0.0));
    }
    assert!(!oneloop::is_initialized());
}

#[test]
fn native_float_validation_precedes_precision_rounding() {
    let mut evaluator =
        NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::A0, 128).unwrap();
    let mut point = [
        Complex::new(
            Float::with_val(512, 1) + Float::with_val(512, 0.5).pow(200),
            Float::new(512),
        ),
        Complex::new(Float::with_val(512, 1), Float::new(512)),
    ];
    let sentinel = Complex::new(Float::with_val(512, 17), Float::with_val(512, -9));
    let mut output = [sentinel.clone(), sentinel.clone(), sentinel.clone()];
    evaluator.evaluate(&point, &mut output).unwrap();
    // The fixed 128-bit evaluator rounds the supplied 512-bit mass to one.
    assert_eq!(output[1].re, Float::with_val(128, 1));
    let tiny = Float::parse("1e-1000", Some(512)).unwrap();
    point[0].im = tiny.clone();
    output.fill(sentinel.clone());
    assert!(evaluator.evaluate(&point, &mut output).is_err());
    assert!(output.iter().all(|value| value == &sentinel));
    point[0].im = -tiny;
    evaluator.evaluate(&point, &mut output).unwrap();
    assert_eq!(output[1].im, point[0].im.clone().at_precision(128));
    assert!(!oneloop::is_initialized());
}

#[test]
fn native_double_float_keeps_compensated_mass() {
    let mut evaluator = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::A0).unwrap();
    let mass = DoubleFloat::from_compensated_sum(2.0, 1e-20);
    let point = [
        Complex::new(mass, DoubleFloat::from(0.0)),
        Complex::new(DoubleFloat::from(1.0), DoubleFloat::from(0.0)),
    ];
    let mut output = [Complex::new(DoubleFloat::from(0.0), DoubleFloat::from(0.0)); 3];
    evaluator.evaluate(&point, &mut output).unwrap();
    assert_eq!(output[1].re, mass);
    assert_ne!(output[1].re, DoubleFloat::from(2.0));
    assert!(!oneloop::is_initialized());
}
