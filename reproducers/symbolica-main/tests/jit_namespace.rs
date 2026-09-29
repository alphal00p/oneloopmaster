#![cfg(feature = "native_code_generation")]

use symbolica::{
    atom::AtomCore,
    domains::float::Complex,
    evaluate::{FunctionMap, FunctionRegistrationOptions, InliningPolicy, JITCompilationSettings},
    symbol,
};

#[test]
fn jit_registers_user_helpers_in_symbolica_namespace() {
    let x = symbol!("symbolica::jit_namespace_x");
    let f = symbol!("symbolica::jit_namespace_f");
    let g = symbol!("symbolica::jit_namespace_g");
    assert!(f.is_builtin());
    let mut functions = FunctionMap::new();
    let options = FunctionRegistrationOptions::new().inlining(InliningPolicy::Never);
    functions
        .add_function_with_options(f, vec![x], x.to_atom() + 1, options.clone())
        .unwrap();
    functions
        .add_function_with_options(g, vec![x], f.call((x,)) + 2, options)
        .unwrap();
    let exact = g
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
        real.evaluate(&[4.], &mut output);
        assert_eq!(output, [7.]);
        let mut complex = exact.jit_compile::<Complex<f64>>(settings.clone()).unwrap();
        let mut output = [Complex::new(0., 0.)];
        complex.evaluate(&[Complex::new(4., 2.)], &mut output);
        assert_eq!(output, [Complex::new(7., 2.)]);
        let mut vector = exact.jit_compile::<wide::f64x4>(settings).unwrap();
        let input = wide::f64x4::new([1., 2., 3., 4.]);
        let mut output = [wide::f64x4::ZERO];
        vector.evaluate(&[input], &mut output);
        assert_eq!(output, [input + wide::f64x4::splat(3.)]);
    }
}
