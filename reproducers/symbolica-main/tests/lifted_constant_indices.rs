use symbolica::{
    atom::{Atom, AtomCore, AtomView, EvalFn, EvaluationInfo},
    domains::float::Float,
    evaluate::{FunctionMap, FunctionRegistrationOptions, InliningPolicy},
    parse, symbol,
};

#[test]
fn lifted_constants_preserve_nested_external_function_indices_and_tags() {
    let constant = symbol!(
        "symbolica::constant_indices::constant",
        eval = EvaluationInfo::constant(|tags, prec| {
            let value = if tags[0] == 1 { 3 } else { 5 };
            Ok(Float::with_val(prec, value).into())
        })
        .with_tags(1)
    );
    let scale = symbol!(
        "symbolica::constant_indices::scale",
        eval = EvaluationInfo::new().with_tags(1).register_tagged(
            |tags: &[AtomView]| -> EvalFn<f64> {
                let factor = if tags[0] == 1 { 2.0 } else { 7.0 };
                Box::new(move |args: &[f64]| factor * args[0])
            }
        )
    );
    let x = symbol!("symbolica::constant_indices::x");
    let inner = symbol!("symbolica::constant_indices::inner");
    let outer = symbol!("symbolica::constant_indices::outer");
    let options = FunctionRegistrationOptions::new().inlining(InliningPolicy::Never);
    let parameters = [x.to_atom()];

    // Exercise mixed tables, a table containing only constants, and a table
    // without constants. The two tags must keep their distinct implementations.
    let bodies = [
        parse!("pi")
            + constant.call((1,))
            + constant.call((2,))
            + scale.call((1, x))
            + scale.call((2, x)),
        parse!("pi") + constant.call((1,)) + constant.call((2,)),
        scale.call((1, x)) + scale.call((2, x)),
    ];
    for (case, body) in bodies.into_iter().enumerate() {
        let mut functions = FunctionMap::new();
        functions
            .add_tagged_function_with_options(
                inner,
                vec![Atom::num(11)],
                vec![x],
                body,
                options.clone(),
            )
            .unwrap();
        functions
            .add_function_with_options(
                outer,
                vec![x],
                inner.call((11, x)) + inner.call((11, x.to_atom() + 1)),
                options.clone(),
            )
            .unwrap();
        let evaluator = outer
            .call((x,))
            .evaluator(&parameters)
            .function_map(functions)
            .horner_iterations(0)
            .build()
            .unwrap();

        // Before the fix this follows a stale callback index and panics.
        assert_eq!(evaluator.export_instructions().output_count, 1);
        let mut numerical = evaluator.map_coeff(&|c| c.re.to_f64());
        let expected = |x: f64| match case {
            0 => 2.0 * (std::f64::consts::PI + 8.0) + 9.0 * (2.0 * x + 1.0),
            1 => 2.0 * (std::f64::consts::PI + 8.0),
            _ => 9.0 * (2.0 * x + 1.0),
        };
        for x in [-0.5, 0.0, 2.0] {
            assert!((numerical.evaluate_single(&[x]) - expected(x)).abs() < 1e-12);
        }

        #[cfg(feature = "bincode")]
        {
            let bytes = bincode::encode_to_vec(&numerical, bincode::config::standard()).unwrap();
            let (mut restored, consumed): (symbolica::evaluate::ExpressionEvaluator<f64>, _) =
                bincode::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
            assert_eq!(consumed, bytes.len());
            for x in [-0.5, 0.0, 2.0] {
                assert!((restored.evaluate_single(&[x]) - expected(x)).abs() < 1e-12);
            }
        }

        #[cfg(feature = "native_code_generation")]
        {
            let mut compiled = numerical
                .jit_compile(
                    symbolica::evaluate::JITCompilationSettings::new()
                        .direct_translation(false)
                        .with_option("use_threads", "false")
                        .with_option("use_simd", "false"),
                )
                .unwrap();
            for x in [-0.5, 0.0, 2.0] {
                let mut output = [0.0];
                compiled.evaluate(&[x], &mut output);
                assert!((output[0] - expected(x)).abs() < 1e-12);
            }
        }
    }
}
