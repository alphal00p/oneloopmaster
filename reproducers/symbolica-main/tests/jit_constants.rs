#![cfg(feature = "native_code_generation")]

use symbolica::{
    atom::{Atom, AtomCore, EvaluationInfo},
    domains::float::{Complex, Float},
    evaluate::JITCompilationSettings,
    parse, symbol,
};

#[test]
fn exact_jit_resolves_registered_and_fixed_argument_constants() {
    let tagged = symbol!(
        "symbolica::jit_constants::tagged",
        eval = EvaluationInfo::constant(|tags, bits| {
            Ok(Float::with_val(bits, if tags[0] == 3 { 11 } else { 17 }).into())
        })
        .with_tags(1)
    );
    let cases = [
        (parse!("pi"), std::f64::consts::PI),
        (parse!("polylog(2,1/4)"), 0.2676526390827326),
        (tagged.call((3,)) + tagged.call((5,)), 28.),
        (parse!("pi + 1/2"), std::f64::consts::PI + 0.5),
    ];
    for (expression, expected) in cases {
        let exact = expression.evaluator(&[] as &[Atom]).build().unwrap();
        for direct in [false, true] {
            let settings = JITCompilationSettings::new()
                .direct_translation(direct)
                .with_option("use_simd", "false")
                .with_option("use_threads", "false");
            let mut real = exact.jit_compile::<f64>(settings.clone()).unwrap();
            let mut output = [0.];
            real.evaluate(&[], &mut output);
            assert!(
                (output[0] - expected).abs() < 1e-14,
                "{expression}, direct={direct}: got {}, expected {expected}",
                output[0]
            );
            let mut complex = exact.jit_compile::<Complex<f64>>(settings).unwrap();
            let mut output = [Complex::new(0., 0.)];
            complex.evaluate(&[], &mut output);
            assert!((output[0].re - expected).abs() < 1e-14);
            assert_eq!(output[0].im, 0.);
        }
    }
}

#[test]
fn exact_jit_resolves_complex_constants_and_reports_domain_errors() {
    let imaginary = symbol!(
        "symbolica::jit_constants::imaginary",
        eval = EvaluationInfo::constant(|_, bits| {
            Ok(Complex::new(
                Float::with_val(bits, 2),
                Float::with_val(bits, 3),
            ))
        })
    );
    let exact = imaginary
        .to_atom()
        .evaluator(&[] as &[Atom])
        .build()
        .unwrap();
    let settings = JITCompilationSettings::new()
        .with_option("use_simd", "false")
        .with_option("use_threads", "false");
    let mut complex = exact.jit_compile::<Complex<f64>>(settings.clone()).unwrap();
    let mut output = [Complex::new(0., 0.)];
    complex.evaluate(&[], &mut output);
    assert_eq!(output, [Complex::new(2., 3.)]);
    assert!(exact.jit_compile::<f64>(settings).is_err());
}

#[test]
fn exact_jit_resolves_constants_lifted_from_helpers_in_supported_domains() {
    use symbolica::evaluate::{FunctionMap, FunctionRegistrationOptions, InliningPolicy};
    let x = symbol!("symbolica::jit_constants::x");
    let f = symbol!("symbolica::jit_constants::f");
    let mut functions = FunctionMap::new();
    functions
        .add_function_with_options(
            f,
            vec![x],
            parse!("pi") + x,
            FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
        )
        .unwrap();
    let exact = f
        .call((x,))
        .evaluator(&[x.to_atom()])
        .function_map(functions)
        .build()
        .unwrap();
    for direct in [false, true] {
        let settings = JITCompilationSettings::new()
            .direct_translation(direct)
            .with_option("use_threads", "false");
        let mut real = exact.jit_compile::<f64>(settings.clone()).unwrap();
        let mut output = [0.];
        real.evaluate(&[2.], &mut output);
        assert_eq!(output, [2. + std::f64::consts::PI]);
        let mut complex = exact.jit_compile::<Complex<f64>>(settings.clone()).unwrap();
        let mut output = [Complex::new(0., 0.)];
        complex.evaluate(&[Complex::new(2., 3.)], &mut output);
        assert_eq!(output, [Complex::new(2. + std::f64::consts::PI, 3.)]);
        let input = wide::f64x4::new([0., 1., 2., 3.]);
        let expected = input + wide::f64x4::splat(std::f64::consts::PI);
        let mut real = exact.jit_compile::<wide::f64x4>(settings.clone()).unwrap();
        let mut output = [wide::f64x4::ZERO];
        real.evaluate(&[input], &mut output);
        assert_eq!(output, [expected]);
    }
}

#[test]
fn exact_jit_propagates_constant_evaluation_errors() {
    let broken = symbol!(
        "symbolica::jit_constants::broken",
        eval = EvaluationInfo::constant(|_, _| Err("constant unavailable".to_owned()))
    );
    let exact = broken.to_atom().evaluator(&[] as &[Atom]).build().unwrap();
    let result = exact.jit_compile::<f64>(JITCompilationSettings::new());
    assert!(matches!(result, Err(message) if message == "constant unavailable"));
}
