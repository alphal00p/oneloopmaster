//! Exact-expression inspection; these are construction tests, not parity claims.
use oneloop::{
    EvaluationBackend, ExpressionOptions, PrecisionEvaluator, ScalarIntegral,
    get_expression_for_family as get_expression,
    get_expression_for_family_with_options as get_expression_with_options,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    domains::float::{Complex, Float, SingleFloat},
    parser::ParseSettings,
    printer::PrintOptions,
};

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

fn count(expression: &Atom, symbol: Symbol) -> usize {
    let mut count = 0;
    expression.visitor(&mut |view| {
        if view.get_symbol() == Some(symbol) {
            count += 1;
        }
        true
    });
    count
}

fn assert_complete(expression: &Atom) {
    expression.visitor(&mut |view| {
        if let Some(symbol) = view.get_symbol() {
            assert!(
                !symbol.get_name().contains("::__olo_"),
                "unexpanded {symbol}"
            );
            assert!(!matches!(
                symbol.get_name(),
                "oneloopmaster::A0"
                    | "oneloopmaster::B0"
                    | "oneloopmaster::dB0"
                    | "oneloopmaster::C0"
                    | "oneloopmaster::D0"
            ));
        }
        true
    });
}

#[test]
fn all_five_families_have_complete_exact_vacuum_bodies() {
    on_stack(|| {
        for family in [
            ScalarIntegral::A0,
            ScalarIntegral::B0,
            ScalarIntegral::DB0,
            ScalarIntegral::C0,
            ScalarIntegral::D0,
        ] {
            let momenta = match family {
                ScalarIntegral::A0 => 0,
                ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
                ScalarIntegral::C0 => 3,
                ScalarIntegral::D0 => 6,
            };
            let input = (0..family.arity())
                .map(|i| Atom::num(if i < momenta { 0 } else { 1 }))
                .collect::<Vec<_>>();
            let result = get_expression(family, &input).unwrap();
            let expected = match family {
                ScalarIntegral::A0 => [Atom::num(1), Atom::num(1), Atom::num(0)],
                ScalarIntegral::B0 => [Atom::num(0), Atom::num(1), Atom::num(0)],
                ScalarIntegral::DB0 => [Atom::num((1, 6)), Atom::num(0), Atom::num(0)],
                ScalarIntegral::C0 => [Atom::num((-1, 2)), Atom::num(0), Atom::num(0)],
                ScalarIntegral::D0 => [Atom::num((1, 6)), Atom::num(0), Atom::num(0)],
            };
            assert_eq!(result.coefficients(), &expected, "{}", family.name());
            for coefficient in result.coefficients() {
                assert_complete(coefficient);
            }
        }
    });
}

#[test]
fn native_real_attributes_simplify_but_do_not_discard_zero_branches() {
    on_stack(|| {
        let untyped = symbolica::symbol!("inspect_untyped_mass").to_atom();
        let positive = symbolica::symbol!("inspect_positive_mass"; Positive).to_atom();
        let scale =
            symbolica::symbol!("inspect_positive_scale"; Positive; tags = ["inspection::scale"])
                .to_atom();
        let general = get_expression(ScalarIntegral::A0, &[untyped, scale.clone()]).unwrap();
        let simple = ScalarIntegral::A0
            .get_expression(&[positive.clone(), scale])
            .unwrap();
        assert!(count(&general.coefficients()[0], Symbol::CONJ) > 0);
        assert_eq!(count(&simple.coefficients()[0], Symbol::CONJ), 0);
        assert!(
            count(&simple.coefficients()[0], Symbol::IF) > 0,
            "positive is not used as a nonzero proof"
        );
        assert_eq!(simple.coefficients()[1], positive);
        assert_complete(&simple.coefficients()[0]);
        for expression in simple.coefficients() {
            let printed = expression.printer(PrintOptions::full()).to_string();
            assert_eq!(
                Atom::parse(printed, "unused", ParseSettings::default()).unwrap(),
                *expression
            );
        }
        let tagged =
            symbolica::symbol!("inspect_tagged_not_real", tags = ["inspection::real"]).to_atom();
        let tagged = get_expression(ScalarIntegral::A0, &[tagged, Atom::num(1)]).unwrap();
        assert!(
            count(&tagged.coefficients()[0], Symbol::CONJ) > 0,
            "a tag is not a Real attribute"
        );
    });
}

#[test]
fn one_mass_triangle_accepts_real_expressions_and_prunes_specialized_regions() {
    on_stack(|| {
        let [x, y, z] = ["inspect_c0_x", "inspect_c0_y", "inspect_c0_z"]
            .map(|name| symbolica::symbol!(name; Real).to_atom());
        let m = symbolica::symbol!("inspect_c0_m_squared"; Positive).to_atom();
        let mu = symbolica::symbol!("inspect_c0_mu_squared"; Positive).to_atom();
        let generic = [
            x,
            y,
            z.clone(),
            Atom::num(0),
            m.clone(),
            Atom::num(0),
            mu.clone(),
        ];
        assert!(
            get_expression_with_options(
                ScalarIntegral::C0,
                &generic,
                ExpressionOptions {
                    max_nodes: 100,
                    max_depth: 512,
                }
            )
            .unwrap_err()
            .contains("max_nodes")
        );
        // Retain a symbolic real invariant and mass/scale, not just a fully
        // numerical vacuum. The middle mass rotates this to the finite
        // one-mass sector with canonical momenta (z, 0, 0).
        let specialized = [
            Atom::num(0),
            Atom::num(0),
            z,
            Atom::num(0),
            m,
            Atom::num(0),
            mu,
        ];
        let start = std::time::Instant::now();
        let result = get_expression(ScalarIntegral::C0, &specialized).unwrap();
        let mut nodes = 0;
        for expression in result.coefficients() {
            expression.visitor(&mut |_| {
                nodes += 1;
                true
            });
            assert_complete(expression);
            let printed = expression.printer(PrintOptions::full()).to_string();
            assert_eq!(
                Atom::parse(printed, "unused", ParseSettings::default()).unwrap(),
                *expression
            );
        }
        eprintln!(
            "complete real one-mass C0: {nodes} output nodes, {:?}",
            start.elapsed()
        );
        assert!(nodes < 2_000, "the analytic sector must remain compact");
    });
}

#[test]
fn compact_real_one_mass_triangle_matches_independent_mapped_regions() {
    on_stack(|| {
        const BITS: u32 = 384;
        // The original nested-root representation can lose about 196 working
        // bits at s/m=1e-30. Its independent reference therefore receives more
        // precision; the compact expression and 1e-60 target are unchanged.
        const REFERENCE_BITS: u32 = 512;
        let s = symbolica::symbol!("inspect_dense_s"; Real).to_atom();
        let m = symbolica::symbol!("inspect_dense_m"; Positive).to_atom();
        let mu = symbolica::symbol!("inspect_dense_mu"; Positive).to_atom();
        let parameters = [s.clone(), m.clone(), mu.clone()];
        let zero = Atom::num(0);
        let body = get_expression(
            ScalarIntegral::C0,
            &[zero.clone(), zero.clone(), s, zero.clone(), m, zero, mu],
        )
        .unwrap();
        let exact = Atom::evaluator_multiple(
            &body
                .coefficients()
                .iter()
                .map(Atom::as_view)
                .collect::<Vec<_>>(),
            &parameters,
        )
        .direct_translation(true)
        .build()
        .unwrap();
        let converter = Float::new(BITS);
        let mut compact = exact.map_coeff_with_prec(
            &|value| {
                Complex::new(
                    converter.from_rational(&value.re),
                    converter.from_rational(&value.im),
                )
            },
            BITS,
        );
        // Deliberately use the original FunctionMap expression backend, not
        // get_expression or the generated Rust specialization, as the control.
        let mut reference = PrecisionEvaluator::with_binary_precision_and_backend(
            ScalarIntegral::C0,
            REFERENCE_BITS,
            EvaluationBackend::Expression,
        )
        .unwrap();
        let number = |text: &str| Float::parse(text, Some(BITS)).unwrap();
        let real = |value: Float| Complex::new(value, Float::new(BITS));
        let blank = || std::array::from_fn::<_, 3, _>(|_| real(number("0")));
        let tolerance = number("1e-60");
        let absolute = |value: Float| if value.is_negative() { -value } else { value };
        let mut rows = 0;
        for mass in ["0.125", "1", "32"] {
            for scale in ["0.03125", "1", "1024"] {
                for ratio in [
                    "-1000000", "-100", "-4", "-1.00001", "-1", "-0.99999", "-0.01", "-1e-30", "0",
                    "1e-30", "0.1", "1", "2", "4", "100", "1000000",
                ] {
                    for mass_is_zero in [false, true] {
                        // Use the nonzero base mass to set s even when checking
                        // the explicit mass=0 branch, so it is not all vacuum.
                        let momentum = number(ratio) * number(mass);
                        let mass = if mass_is_zero {
                            number("0")
                        } else {
                            number(mass)
                        };
                        let point = [
                            real(momentum.clone()),
                            real(mass.clone()),
                            real(number(scale)),
                        ];
                        let mut actual = blank();
                        compact.evaluate(&point, &mut actual);
                        for mass_index in 0..3 {
                            let mut input = vec![real(number("0")); 7];
                            input[(mass_index + 1) % 3] = real(momentum.clone());
                            input[mass_index + 3] = real(mass.clone());
                            input[6] = real(number(scale));
                            let mut expected = blank();
                            reference.evaluate(&input, &mut expected).unwrap();
                            for (index, (a, b)) in actual.iter().zip(&expected).enumerate() {
                                for (part, a, b) in
                                    [("real", &a.re, &b.re), ("imaginary", &a.im, &b.im)]
                                {
                                    assert!(
                                        a.is_finite() && b.is_finite(),
                                        "nonfinite ratio={ratio}, mass={mass}, scale={scale}, mass_index={mass_index}, coefficient={index} {part}: {a} vs {b}"
                                    );
                                    let error = absolute(a.clone() - b.clone());
                                    let scale = Float::with_val(BITS, 1) + absolute(b.clone());
                                    assert!(
                                        error < tolerance.clone() * scale,
                                        "ratio={ratio}, mass={mass}, mu_squared={}, mass_index={mass_index}, coefficient={index} {part}: {a} vs {b}; error={error:e}",
                                        point[2]
                                    );
                                }
                            }
                            rows += 1;
                        }
                    }
                }
            }
        }
        eprintln!(
            "compact real one-mass C0: {rows} rows, compact {BITS} bits / original reference {REFERENCE_BITS} bits, all three coefficients"
        );
    });
}

#[test]
fn near_zero_reference_converges_without_relaxing_accuracy() {
    on_stack(|| {
        let s = symbolica::symbol!("inspect_convergence_s"; Real).to_atom();
        let m = symbolica::symbol!("inspect_convergence_m"; Positive).to_atom();
        let mu = symbolica::symbol!("inspect_convergence_mu"; Positive).to_atom();
        let parameters = [s.clone(), m.clone(), mu.clone()];
        let zero = Atom::num(0);
        let body = get_expression(
            ScalarIntegral::C0,
            &[zero.clone(), s, zero.clone(), m, zero.clone(), zero, mu],
        )
        .unwrap();
        let exact = Atom::evaluator_multiple(
            &body
                .coefficients()
                .iter()
                .map(Atom::as_view)
                .collect::<Vec<_>>(),
            &parameters,
        )
        .direct_translation(true)
        .build()
        .unwrap();
        let mut values = Vec::new();
        for bits in [384, 512, 768] {
            let number = |text| Float::parse(text, Some(bits)).unwrap();
            let real = |value| Complex::new(value, Float::new(bits));
            let point = [
                real(number("-1.25e-31")),
                real(number("0.125")),
                real(number("0.03125")),
            ];
            let converter = Float::new(bits);
            let mut compact = exact.clone().map_coeff_with_prec(
                &|value| {
                    Complex::new(
                        converter.from_rational(&value.re),
                        converter.from_rational(&value.im),
                    )
                },
                bits,
            );
            let mut actual = std::array::from_fn::<_, 3, _>(|_| real(number("0")));
            compact.evaluate(&point, &mut actual);
            let mut reference = PrecisionEvaluator::with_binary_precision_and_backend(
                ScalarIntegral::C0,
                bits,
                EvaluationBackend::Expression,
            )
            .unwrap();
            let input = [
                real(number("0")),
                point[0].clone(),
                real(number("0")),
                point[1].clone(),
                real(number("0")),
                real(number("0")),
                point[2].clone(),
            ];
            let mut expected = std::array::from_fn::<_, 3, _>(|_| real(number("0")));
            reference.evaluate(&input, &mut expected).unwrap();
            assert!(actual[0].re.is_finite() && expected[0].re.is_finite());
            eprintln!(
                "near-zero C0 at {bits} bits: compact precision={}, original precision={}, difference={:e}",
                actual[0].re.prec(),
                expected[0].re.prec(),
                actual[0].re.clone() - expected[0].re.clone()
            );
            values.push((actual[0].re.clone(), expected[0].re.clone()));
        }
        let absolute = |value: Float| if value.is_negative() { -value } else { value };
        let bound = Float::parse("1e-60", Some(768)).unwrap()
            * (Float::with_val(768, 1) + absolute(values[2].1.clone()));
        // Require convergence of the independent reference and the original
        // 384-bit compact result. Do not require the lower-precision reference
        // to remain inaccurate: future engine improvements are welcome.
        assert!(absolute(values[1].1.clone() - values[2].1.clone()) < bound);
        assert!(absolute(values[0].0.clone() - values[2].1.clone()) < bound);
    });
}

#[test]
fn exact_massless_regions_and_resource_errors_are_explicit() {
    on_stack(|| {
        let a = get_expression(ScalarIntegral::A0, &[Atom::num(0), Atom::num(1)]).unwrap();
        assert!(a.coefficients().iter().all(|c| c == &Atom::num(0)));
        let input = [Atom::num(-2), Atom::num(0), Atom::num(0), Atom::num(1)];
        let b = get_expression(ScalarIntegral::B0, &input).unwrap();
        for coefficient in b.coefficients() {
            assert_complete(coefficient);
        }
        assert_eq!(b.coefficients()[1], 1);
        assert_eq!(b.coefficients()[2], 0);
        assert!(
            get_expression(ScalarIntegral::D0, &input)
                .unwrap_err()
                .contains("expects 11")
        );
        assert!(
            get_expression_with_options(
                ScalarIntegral::B0,
                &input,
                ExpressionOptions {
                    max_nodes: 1,
                    max_depth: 512
                }
            )
            .unwrap_err()
            .contains("max_nodes")
        );
        assert!(
            get_expression_with_options(
                ScalarIntegral::B0,
                &input,
                ExpressionOptions {
                    max_nodes: 1_000_000,
                    max_depth: 1
                }
            )
            .unwrap_err()
            .contains("max_depth")
        );
        assert!(
            get_expression_with_options(
                ScalarIntegral::B0,
                &input,
                ExpressionOptions {
                    max_nodes: 0,
                    max_depth: 512
                }
            )
            .is_err()
        );
        // A real squared expression can still vanish. Keep its scaleless branch.
        let real = symbolica::symbol!("inspect_real_square"; Real).to_atom();
        let result = get_expression(ScalarIntegral::A0, &[real.pow(2), Atom::num(1)]).unwrap();
        assert!(count(&result.coefficients()[0], Symbol::IF) > 0);
        assert!(matches!(
            result.coefficients()[2].as_view(),
            AtomView::Num(_)
        ));
    });
}
