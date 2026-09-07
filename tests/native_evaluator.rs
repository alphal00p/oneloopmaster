use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    evaluate::{FunctionMap, FunctionRegistrationOptions, InliningPolicy},
    parse,
};

#[test]
#[ignore = "targeted precision diagnostic for the exploratory scalar audit"]
fn audit_near_zero_bubble_precision() {
    use symbolica::prelude::{Complex, Float, RealLike};
    let parameters = [
        parse!("audit_p"),
        parse!("audit_m0"),
        parse!("audit_m1"),
        parse!("audit_mu"),
    ];
    let series = oneloop::db0(
        &parameters[0],
        &parameters[1],
        &parameters[2],
        &parameters[3],
    );
    let exact = series.coefficients()[0]
        .evaluator(&parameters)
        .build()
        .unwrap();
    for bits in [128, 256] {
        let mut evaluator = exact.clone().map_coeff_with_prec(
            &|v| {
                Complex::new(
                    v.re.to_multi_prec_float(bits),
                    v.im.to_multi_prec_float(bits),
                )
            },
            bits,
        );
        let p = Float::with_val(bits, 1) / Float::with_val(bits, 10_000_000_000_i64);
        let input = [
            p,
            Float::with_val(bits, 1),
            Float::with_val(bits, 2),
            Float::with_val(bits, 1),
        ]
        .map(|v| Complex::new(v, Float::new(bits)));
        let value = evaluator.evaluate_single(&input);
        eprintln!("near-zero dB0 at {bits} bits: {value}");
        // This is an explicit precision diagnostic, not an adaptive evaluator:
        // 128 bits still loses digits to cancellation for this representation.
        if bits == 256 {
            assert!((value.re.to_f64() - 0.11370563888166617).abs() < 1e-12);
            assert!(value.im.to_f64().abs() < 1e-24);
        }
    }
}

#[test]
fn polylog_preserves_high_precision_endpoint_neighborhoods() {
    use symbolica::{
        domains::backend::float::Constant,
        prelude::{Complex, Float, FloatLike, Real, RealLike},
    };
    let x = parse!("polylog_endpoint_x");
    let exact = parse!("polylog(2,polylog_endpoint_x)")
        .evaluator(&[x])
        .build()
        .unwrap();
    for bits in [128, 256] {
        let mut evaluator = exact.clone().map_coeff_with_prec(
            &|v| {
                Complex::new(
                    v.re.to_multi_prec_float(bits),
                    v.im.to_multi_prec_float(bits),
                )
            },
            bits,
        );
        let one = Float::with_val(bits, 1);
        let delta = one.clone() / Float::with_val(bits, 2).pow(80);
        let pi = Float::with_val(bits, Constant::Pi);
        let li2 = |e: &mut symbolica::evaluate::ExpressionEvaluator<Complex<Float>>, z: Float| {
            e.evaluate_single(&[Complex::new(z, Float::new(bits))])
        };
        let small = li2(&mut evaluator, delta.clone());
        let below = li2(&mut evaluator, one.clone() - delta.clone());
        // Reflection identity, evaluated without rounding the input to f64.
        let expected =
            pi.clone().pow(2) / 6 - small.re - delta.log() * (one.clone() - delta.clone()).log();
        assert!((below.re - expected).to_f64().abs() < 1e-34);
        assert!(below.im.to_f64().abs() < 1e-34);

        let negative = li2(&mut evaluator, -one.clone() + delta.clone());
        // First-order expansion at -1; the omitted term is O(delta^2) < 1e-48.
        let expected = -pi.clone().pow(2) / 12 + delta.clone() * Float::with_val(bits, 2).log();
        assert!((negative.re - expected).to_f64().abs() < 1e-34);
        assert!(negative.im.to_f64().abs() < 1e-34);

        let above = li2(&mut evaluator, one.clone() + delta.clone());
        assert!(above.re.to_f64().is_finite());
        assert!((above.im + pi * (one + delta).log()).to_f64().abs() < 1e-34);
    }
}

#[test]
fn native_predicates_and_noninlined_functions() {
    let x = parse!("native_regression_x");
    let f = symbolica::symbol!("native_regression_f");
    let body = parse!("pi+polylog(2,native_regression_x)");
    let mut map = FunctionMap::new();
    map.add_function_with_options(
        f,
        vec![symbolica::symbol!("native_regression_x")],
        body.clone(),
        FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
    )
    .unwrap();
    let mut eval = f
        .call((&x,))
        .evaluator(std::slice::from_ref(&x))
        .function_map(map)
        .build()
        .unwrap()
        .map_coeff(&|v| symbolica::prelude::Complex::new(v.re.to_f64(), v.im.to_f64()));
    let mut reference = body
        .evaluator(std::slice::from_ref(&x))
        .build()
        .unwrap()
        .map_coeff(&|v| symbolica::prelude::Complex::new(v.re.to_f64(), v.im.to_f64()));
    let input = [symbolica::prelude::Complex::new(0.25, 0.)];
    assert_eq!(
        eval.evaluate_single(&input),
        reference.evaluate_single(&input)
    );
    let sign = Symbol::IF.call((&x - Symbol::ABS.call((&x,)), Atom::num(-1), Atom::num(1)));
    let mut eval = sign
        .evaluator(&[x])
        .build()
        .unwrap()
        .map_coeff(&|v| symbolica::prelude::Complex::new(v.re.to_f64(), v.im.to_f64()));
    for value in [-1e150, -3.21, -0.123, 0., 0.123, 3.21, 1e150] {
        assert_eq!(
            eval.evaluate_single(&[symbolica::prelude::Complex::new(value, 0.)])
                .re,
            if value < 0. { -1. } else { 1. }
        );
    }
}
