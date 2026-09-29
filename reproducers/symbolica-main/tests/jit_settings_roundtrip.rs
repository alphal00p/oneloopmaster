#![cfg(all(feature = "native_code_generation", feature = "bincode"))]

use symbolica::{
    atom::{AtomCore, Symbol},
    domains::float::Complex,
    evaluate::{
        FunctionMap, FunctionRegistrationOptions, InliningPolicy, JITCompilationSettings,
        JITCompiledEvaluator,
    },
    symbol,
};

#[test]
fn jit_roundtrips_preserve_nested_compilation_policy() {
    let x = symbol!("symbolica::jit_roundtrip_settings::x");
    let f = symbol!("symbolica::jit_roundtrip_settings::f");
    let g = symbol!("symbolica::jit_roundtrip_settings::g");
    let mut functions = FunctionMap::new();
    let options = FunctionRegistrationOptions::new().inlining(InliningPolicy::Never);
    functions
        .add_function_with_options(
            f,
            vec![x],
            x.to_atom() * Symbol::CONJ.call((x,)),
            options.clone(),
        )
        .unwrap();
    functions
        .add_function_with_options(g, vec![x], f.call((x,)), options)
        .unwrap();
    let exact = g
        .call((x,))
        .evaluator(&[x.to_atom()])
        .function_map(functions)
        .build()
        .unwrap();
    let inputs = [
        Complex::new(0.1, 0.2),
        Complex::new(0.3, -0.7),
        Complex::new(2., 3.),
    ];
    for direct in [false, true] {
        for fast_complex in [false, true] {
            let settings = JITCompilationSettings::new()
                .optimization_level(2)
                .direct_translation(direct)
                .with_option("fastmath", "false")
                .with_option("fast_complex", fast_complex.to_string())
                .with_option("use_simd", "false")
                .with_option("use_threads", "false");
            let mut evaluator = exact.jit_compile::<Complex<f64>>(settings).unwrap();
            let mut expected = [Complex::new(0., 0.); 3];
            evaluator.batch_evaluate(&inputs, &mut expected, 3);
            // The two backends currently round some operations differently.
            // This test checks preservation of the selected backend and policy,
            // independently of either backend's absolute numerical accuracy.
            for _ in 0..2 {
                let bytes =
                    bincode::encode_to_vec(&evaluator, bincode::config::standard()).unwrap();
                let (restored, consumed): (JITCompiledEvaluator<Complex<f64>>, usize) =
                    bincode::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
                assert_eq!(consumed, bytes.len());
                evaluator = restored;
                let mut actual = [Complex::new(0., 0.); 3];
                evaluator.batch_evaluate(&inputs, &mut actual, 3);
                assert_eq!(
                    actual, expected,
                    "direct={direct}, fast_complex={fast_complex}"
                );
                for (input, expected) in inputs.iter().zip(expected) {
                    let mut output = [Complex::new(0., 0.)];
                    evaluator.evaluate(&[*input], &mut output);
                    assert_eq!(output, [expected]);
                }
            }
        }
    }
}
