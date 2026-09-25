//! Early branch selection preserves the same parametric regions as late selection.
use oneloop::{
    B0, C0, ExpressionOptions, get_expression, get_expression_on_branch,
    get_expression_on_branch_with_options, select_branch,
};
use symbolica::domains::float::{Complex, Float};
use symbolica::prelude::*;

#[test]
fn massive_triangle_probe() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let a = symbol!("branch_expansion_a").to_atom();
            let b = symbol!("branch_expansion_b").to_atom();
            let mu = symbol!("branch_expansion_mu").to_atom();
            let master = C0().call((0, -&a, &a, &a, &a, &b, &mu));
            let rules = [
                Replacement::new(a.clone(), 2),
                Replacement::new(b.clone(), 1),
                Replacement::new(mu.clone(), 1),
            ];
            let start = std::time::Instant::now();
            let series = get_expression_on_branch_with_options(
                &master,
                &rules,
                ExpressionOptions {
                    max_nodes: 100_000_000,
                    max_depth: 512,
                },
            )
            .unwrap();
            eprintln!("expansion time: {:?}", start.elapsed());
            for coefficient in series.coefficients() {
                assert_eq!(&select_branch(coefficient, &rules), coefficient);
                coefficient.visitor(&mut |node| {
                    assert_ne!(node.get_symbol(), Some(Symbol::IF));
                    if let Some(symbol) = node.get_symbol() {
                        assert!(!symbol.get_name().contains("::__olo_"));
                    }
                    true
                });
            }
            let finite = &series.coefficients()[0];
            let symbols = finite.get_all_symbols(false);
            assert!(symbols.contains(&a.get_symbol().unwrap()));
            assert!(symbols.contains(&b.get_symbol().unwrap()));
            assert!(
                get_expression_on_branch_with_options(
                    &master,
                    &rules,
                    ExpressionOptions {
                        max_nodes: 10,
                        max_depth: 512,
                    }
                )
                .unwrap_err()
                .contains("max_nodes")
            );
            eprintln!("expression length: {}", finite.to_string().len());
            let precise_inputs = [
                (a.clone(), Complex::new(Float::with_val(256, 2), Float::new(256))),
                (b.clone(), Complex::new(Float::with_val(256, 1), Float::new(256))),
                (mu.clone(), Complex::new(Float::with_val(256, 1), Float::new(256))),
            ].into_iter().collect::<std::collections::HashMap<_, _>>();
            let precise: Complex<Float> = finite.evaluate_with_prec(&precise_inputs, 256).unwrap();
            let expected = Float::parse(
                "-0.31341398580589155641085172932611532693192520090780374821622716949068511831637853",
                Some(256),
            ).unwrap();
            eprintln!("256-bit value: {precise}");
            assert!((precise.re - expected).to_f64().abs() < 1e-65);
            assert!(precise.im.to_f64().abs() < 1e-65);
            let value: Complex<f64> = finite
                .evaluate(
                    &[
                        (a, Complex::new(2., 0.)),
                        (b, Complex::new(1., 0.)),
                        (mu, Complex::new(1., 0.)),
                    ]
                    .into_iter()
                    .collect::<std::collections::HashMap<_, _>>(),
                )
                .unwrap();
            eprintln!("value: {value:?}");
            assert!((value.re + 0.313_413_985_805_891_6).abs() < 1e-12);
            assert!(value.im.abs() < 1e-12);
            assert_eq!(series.coefficients()[1], Atom::num(0));
            assert_eq!(series.coefficients()[2], Atom::num(0));
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn early_selection_matches_late_selection_and_caches_do_not_cross_probes() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let s = symbol!("early_b0_s"; Real).to_atom();
            let m = symbol!("early_b0_m"; Positive).to_atom();
            let master = B0().call((&s, &m, &m, 1));
            let all = get_expression(&master).unwrap();
            for sample in [
                Atom::num(0),
                Atom::num(3),
                Atom::num(4),
                Atom::num(5),
                Atom::num(-2),
                Atom::num((1, 2)),
                Atom::num(Float::with_val(128, 3.23)),
            ] {
                let rules = [
                    Replacement::new(m.clone(), 1),
                    Replacement::new(s.clone(), sample),
                ];
                let early = get_expression_on_branch(&master, &rules).unwrap();
                let late = all
                    .coefficients()
                    .each_ref()
                    .map(|c| select_branch(c, &rules));
                assert_eq!(early.coefficients(), &late);
            }
            for rules in [vec![], vec![Replacement::new(m.clone(), 1)]] {
                let early = get_expression_on_branch(&master, &rules).unwrap();
                let late = all
                    .coefficients()
                    .each_ref()
                    .map(|c| select_branch(c, &rules));
                // Path-aware expansion can remove additional redundant guards
                // while the probes are incomplete. Compare the resulting
                // regions, not the spelling of those unresolved guards.
                for sample in [0, 3, 5, -2] {
                    let complete = [
                        Replacement::new(m.clone(), 1),
                        Replacement::new(s.clone(), sample),
                    ];
                    let early_selected = early
                        .coefficients()
                        .each_ref()
                        .map(|c| select_branch(c, &complete));
                    let late_selected = late.each_ref().map(|c| select_branch(c, &complete));
                    assert_eq!(early_selected, late_selected);
                }
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
