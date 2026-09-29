use symbolica::{
    atom::{Atom, AtomCore, EvaluationInfo, Symbol},
    domains::float::{Complex, Float},
    evaluate::{ExpressionEvaluator, FunctionMap, FunctionRegistrationOptions, InliningPolicy},
    parse, symbol,
};

fn check_constants(exact: &ExpressionEvaluator<Complex<symbolica::domains::rational::Rational>>) {
    let exported = exact.export_instructions();
    assert!(!exported.constant_functions.is_empty());
    for bits in [53, 128, 512] {
        let convert = |c: &Complex<symbolica::domains::rational::Rational>| {
            Complex::new(
                c.re.to_multi_prec_float(bits),
                c.im.to_multi_prec_float(bits),
            )
        };
        let numeric = exact
            .clone()
            .map_coeff_with_prec(&convert, bits)
            .export_instructions();
        let mut reconstructed: Vec<_> = exported.constants.iter().map(convert).collect();
        let mut indices = std::collections::HashSet::new();
        for constant in &exported.constant_functions {
            assert!(indices.insert(constant.index));
            assert!(constant.index < reconstructed.len());
            let tags = constant
                .tags
                .iter()
                .map(|tag| parse!(tag))
                .collect::<Vec<_>>();
            let tags = tags.iter().map(Atom::as_view).collect::<Vec<_>>();
            let info = constant.symbol.get_evaluation_info().unwrap();
            let value = if info.has_constant_evaluator() {
                assert!(constant.fixed_args.is_empty());
                info.evaluate_constant(&tags, bits).unwrap()
            } else {
                let callback = info.get_evaluator::<Complex<Float>>(&tags).unwrap();
                callback(&constant.fixed_args.iter().map(convert).collect::<Vec<_>>())
            };
            reconstructed[constant.index] = value;
        }
        assert_eq!(reconstructed, numeric.constants);
        // Exporting after coefficient conversion must retain the same definitions.
        assert_eq!(
            exported.constant_functions.len(),
            numeric.constant_functions.len()
        );
        for (a, b) in exported
            .constant_functions
            .iter()
            .zip(&numeric.constant_functions)
        {
            assert_eq!(
                (a.index, a.symbol, &a.tags, &a.fixed_args),
                (b.index, b.symbol, &b.tags, &b.fixed_args)
            );
        }
    }
}

#[test]
fn exact_exports_retain_registered_constants_tags_and_fixed_arguments() {
    let tagged = symbol!(
        "symbolica::exported_constants::tagged",
        eval = EvaluationInfo::constant(|tags, bits| {
            Ok(Float::with_val(bits, if tags[0] == 3 { 11 } else { 17 }).into())
        })
        .with_tags(1)
    );
    let expression = parse!("pi + polylog(2,1/4)") + tagged.call((3,)) + tagged.call((5,));
    let exact = expression.evaluator(&[] as &[Atom]).build().unwrap();
    check_constants(&exact);
    let exported = exact.export_instructions();
    assert!(
        exported
            .constant_functions
            .iter()
            .any(|c| c.symbol == Symbol::PI)
    );
    let polylog = exported
        .constant_functions
        .iter()
        .find(|c| c.symbol == symbolica::transcendental::polylog())
        .unwrap();
    assert_eq!(polylog.tags, ["2"]);
    assert_eq!(polylog.fixed_args, [Complex::new((1, 4).into(), 0.into())]);
    let mut tags = exported
        .constant_functions
        .iter()
        .filter(|c| c.symbol == tagged)
        .map(|c| c.tags.clone())
        .collect::<Vec<_>>();
    tags.sort();
    assert_eq!(tags, [vec!["3".to_owned()], vec!["5".to_owned()]]);

    #[cfg(feature = "bincode")]
    {
        let bytes = bincode::encode_to_vec(&exact, bincode::config::standard()).unwrap();
        let (restored, consumed) =
            bincode::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
        assert_eq!(consumed, bytes.len());
        check_constants(&restored);
    }
}

#[test]
fn nested_exports_keep_lifted_constant_definitions() {
    let x = symbol!("symbolica::exported_constants::x");
    let f = symbol!("symbolica::exported_constants::f");
    let mut map = FunctionMap::new();
    let options = FunctionRegistrationOptions::new().inlining(InliningPolicy::Never);
    map.add_function_with_options(
        f,
        vec![x],
        parse!("pi + polylog(2,1/4)") + x,
        options.clone(),
    )
    .unwrap();
    let exact = f
        .call((x,))
        .evaluator(&[x.to_atom()])
        .function_map(map)
        .build()
        .unwrap();
    check_constants(&exact);
    assert!(!exact.export_instructions().sub_evaluators.is_empty());
}
