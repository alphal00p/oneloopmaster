//! Complete expressions must support selection *after* construction.
use oneloop::{C0, ExpressionOptions, get_expression_with_options, select_branch};
use symbolica::prelude::*;

#[test]
fn complete_massive_triangle_then_select() {
    std::thread::Builder::new().stack_size(128 * 1024 * 1024).spawn(|| {
        let a = symbol!("full_c0_mass2").to_atom();
        let b = symbol!("full_c0_mass2B").to_atom();
        let mu = symbol!("full_c0_mu2").to_atom();
        let master = C0().call((0,-&a,&a,&a,&a,&b,&mu));
        let start = std::time::Instant::now();
        let all = get_expression_with_options(&master, ExpressionOptions { max_nodes: 100_000_000, max_depth: 512 }).unwrap();
        eprintln!("full C0 construction: {:?}", start.elapsed());
        let mut nodes = 0;
        let mut branches = 0;
        for coefficient in all.coefficients() {
            coefficient.visitor(&mut |v| {
                nodes += 1;
                branches += usize::from(v.get_symbol() == Some(Symbol::IF));
                if let Some(s) = v.get_symbol() { assert!(!s.get_name().contains("::__olo_")); }
                true
            });
        }
        assert!(branches > 0);
        eprintln!("full C0: {nodes} nodes, {branches} branches");
        let rules = [Replacement::new(a.clone(),2), Replacement::new(b.clone(),1), Replacement::new(mu.clone(),1)];
        let start = std::time::Instant::now();
        let selected = all.coefficients().each_ref().map(|e| select_branch(e,&rules));
        eprintln!("C0 late selection: {:?}; {} characters",start.elapsed(),selected[0].to_string().len());
        for coefficient in &selected {
            coefficient.visitor(&mut |v| { assert_ne!(v.get_symbol(),Some(Symbol::IF)); true });
        }
        assert_eq!(selected[1..], [Atom::num(0),Atom::num(0)]);
        let variables = selected[0].get_all_symbols(false);
        assert!(variables.contains(&a.get_symbol().unwrap()));
        assert!(variables.contains(&b.get_symbol().unwrap()));
        let inputs = [(a.clone(),2),(b.clone(),1),(mu.clone(),1)].into_iter().map(|(x,n)| (x,Complex::new(Float::with_val(256,n),Float::new(256)))).collect::<std::collections::HashMap<_,_>>();
        let value: Complex<Float> = selected[0].evaluate_with_prec(&inputs,256).unwrap();
        let expected = Float::parse("-0.31341398580589155641085172932611532693192520090780374821622716949068511831637853",Some(256)).unwrap();
        eprintln!("C0 selected value: {value}");
        assert!((value.re-expected).to_f64().abs()<1e-65);
        assert!(value.im.to_f64().abs()<1e-65);
        // These also exercise real-axis lips and rounding-sensitive sign
        // predicates: deciding after stripping their grouping used to fail.
        for (av,bv,scale) in [(-3,1,1),(-1,-1,1),(2,-1,1),(2,0,7),(0,1,7),(0,0,1),(2,8,1),(2,12,3)] {
            let rules = [Replacement::new(a.clone(),av),Replacement::new(b.clone(),bv),Replacement::new(mu.clone(),scale)];
            let inputs = [(a.clone(),av),(b.clone(),bv),(mu.clone(),scale)].into_iter().map(|(x,n)| (x,Complex::new(Float::with_val(256,n),Float::new(256)))).collect::<std::collections::HashMap<_,_>>();
            for expression in all.coefficients() {
                let selected = select_branch(expression,&rules);
                selected.visitor(&mut |v| { assert_ne!(v.get_symbol(),Some(Symbol::IF)); true });
                let before: Complex<Float> = expression.evaluate_with_prec(&inputs,256).unwrap();
                let after: Complex<Float> = selected.evaluate_with_prec(&inputs,256).unwrap();
                let error = (before.re-after.re).to_f64().hypot((before.im-after.im).to_f64());
                assert!(error<1e-65,"a={av}, b={bv}, scale={scale}, error={error}");
            }
        }
    }).unwrap().join().unwrap();
}

#[test]
fn all_generic_masters_have_complete_shared_expressions_and_select_late() {
    use oneloop::{
        EvaluationBackend, PrecisionEvaluator, ScalarIntegral as F, get_expression_shared,
        select_branch_shared,
    };
    std::thread::Builder::new().stack_size(128*1024*1024).spawn(|| {
        for (family, symbol, momenta) in [(F::A0,oneloop::A0(),0),(F::B0,oneloop::B0(),1),(F::DB0,oneloop::dB0(),1),(F::C0,oneloop::C0(),3),(F::D0,oneloop::D0(),6)] {
            let parameters = (0..family.arity()).map(|i| symbol!(format!("full_{}_parameter_{i}",family.name())).to_atom()).collect::<Vec<_>>();
            let start = std::time::Instant::now();
            let all = get_expression_shared(symbol.call(&parameters)).unwrap();
            eprintln!("{} full shared construction: {:?}",family.name(),start.elapsed());
            let mut branches = 0;
            for coefficient in all.coefficients() {
                for expression in std::iter::once(coefficient.get_root()).chain(coefficient.get_aliases().values()) {
                    expression.visitor(&mut |v| {
                        if let Some(s) = v.get_symbol() { assert!(!s.get_name().contains("::__olo_"), "unresolved {s}"); }
                        branches += usize::from(v.get_symbol() == Some(Symbol::IF));
                        true
                    });
                }
            }
            assert!(branches>0);
            let euclidean = (0..family.arity()).map(|i| Atom::num(if i<momenta { -(i as i32)-1 } else { i as i32+1 })).collect::<Vec<_>>();
            let mut complex = euclidean.clone();
            for value in &mut complex[momenta..family.arity()-1] { *value = &*value - parse!("1i/4"); }
            let infrared: Vec<Atom> = match family {
                F::A0 => vec![0,2], F::B0|F::DB0 => vec![1,0,0,2],
                F::C0 => vec![0,0,-3,0,0,0,2],
                F::D0 => vec![0,0,0,0,-3,-5,0,0,0,0,2],
            }.into_iter().map(Atom::num).collect();
            let mut reference = PrecisionEvaluator::with_binary_precision_and_backend(family,384,EvaluationBackend::Native).unwrap();
            for (region,point) in [euclidean,complex,infrared].into_iter().enumerate() {
                let rules = parameters.iter().zip(&point).map(|(x,v)| {
                    let rule = Replacement::new(x.clone(),v.clone());
                    // Match Python's ordinary Replacement defaults as well as Rust's.
                    if region == 1 { rule.when(symbolica::id::Condition::True).rhs_cache_size(100) } else { rule }
                }).collect::<Vec<_>>();
                let start = std::time::Instant::now();
                let selected = all.coefficients().each_ref().map(|coefficient| select_branch_shared(coefficient,&rules,ExpressionOptions { max_nodes: 100_000_000,max_depth:4096 }).unwrap());
                eprintln!("{} region {region} late selection: {:?}; {} characters",family.name(),start.elapsed(),selected[0].to_string().len());
                let values = point.iter().map(|x| x.evaluate_with_prec::<Atom,_>(&std::collections::HashMap::new(),384).unwrap()).collect::<Vec<Complex<Float>>>();
                let inputs = parameters.iter().cloned().zip(values.iter().cloned()).collect::<std::collections::HashMap<_,_>>();
                let mut expected = std::array::from_fn::<_,3,_>(|_| Complex::new(Float::new(384),Float::new(384)));
                reference.evaluate(&values,&mut expected).unwrap();
                for (coefficient,(expression,expected)) in selected.iter().zip(&expected).enumerate() {
                    expression.visitor(&mut |v| { assert_ne!(v.get_symbol(),Some(Symbol::IF)); true });
                    let actual: Complex<Float> = expression.evaluate_with_prec(&inputs,384).unwrap();
                    let error = (actual.re.clone()-&expected.re).to_f64().hypot((actual.im.clone()-&expected.im).to_f64());
                    assert!(error<1e-45,"{} region {region}, coefficient {coefficient}: {actual} != {expected}; error={error}",family.name());
                }
            }
        }
    }).unwrap().join().unwrap();
}

#[test]
fn shared_bindings_compile_natively_and_native_patterns_keep_their_scope() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let x = symbol!("shared_native_x").to_atom();
            let mu = symbol!("shared_native_mu").to_atom();
            let master = oneloop::A0().call((&x, &mu));
            let graph = oneloop::get_expression_shared(&master).unwrap();
            let ordinary = oneloop::get_expression(&master).unwrap();
            let mut evaluator = graph.coefficients()[0]
                .evaluator(&[x.clone(), mu.clone()])
                .direct_translation(true)
                .horner_iterations(0)
                .build()
                .unwrap()
                .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
            let value = evaluator.evaluate_single(&[Complex::new(2., 0.), Complex::new(1., 0.)]);
            assert!((value.re - 2. * (1. - 2f64.ln())).abs() < 1e-14);
            assert!(value.im.abs() < 1e-14);
            for rules in [
                vec![
                    Replacement::new(&x / &mu, 2),
                    Replacement::new(x.clone(), 4),
                ],
                vec![
                    Replacement::new(x.clone(), parse!("1e-1000")),
                    Replacement::new(mu.clone(), 1),
                ],
            ] {
                let selected = oneloop::select_branch_shared(
                    &graph.coefficients()[0],
                    &rules,
                    ExpressionOptions::default(),
                )
                .unwrap();
                assert_eq!(selected, select_branch(&ordinary.coefficients()[0], &rules));
                selected.visitor(&mut |v| {
                    assert_ne!(v.get_symbol(), Some(Symbol::IF));
                    true
                });
            }
            let mut limited = Replacement::new(x.clone(), 1);
            limited.match_settings = symbolica::id::MatchSettings::default().max_level(0);
            let rules = [limited];
            let selected = oneloop::select_branch_shared(
                &graph.coefficients()[0],
                &rules,
                ExpressionOptions::default(),
            )
            .unwrap();
            let ordinary = select_branch(&ordinary.coefficients()[0], &rules);
            let complete = [Replacement::new(x, 2), Replacement::new(mu, 1)];
            assert_eq!(
                select_branch(selected, &complete),
                select_branch(ordinary, &complete)
            );
        })
        .unwrap()
        .join()
        .unwrap();
}
