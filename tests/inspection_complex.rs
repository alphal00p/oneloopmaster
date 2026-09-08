//! Continuation checks for the inspection-only one-mass triangle identity.
use oneloop::{EvaluationBackend, PrecisionEvaluator, ScalarIntegral, get_expression};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::{Complex, Float, SingleFloat},
};

#[test]
fn complex_one_mass_triangle_matches_original_mapped_continuation() {
    std::thread::Builder::new().stack_size(128 * 1024 * 1024).spawn(|| {
        const BITS: u32 = 512;
        let s = symbolica::symbol!("inspect_complex_s"; Real).to_atom();
        // Intentionally no Real/Positive mass attribute: complex is the default.
        let m = symbolica::symbol!("inspect_complex_m").to_atom();
        let mu = symbolica::symbol!("inspect_complex_mu"; Positive).to_atom();
        let parameters = [s.clone(), m.clone(), mu.clone()];
        let zero = Atom::num(0);
        let bodies = get_expression(ScalarIntegral::C0, &[zero.clone(),zero.clone(),s,zero.clone(),m,zero,mu]).unwrap();
        let mut nodes = 0;
        for body in bodies.coefficients() {
            body.visitor(&mut |view| {
                nodes += 1;
                if let Some(symbol) = view.get_symbol() {
                    assert!(!symbol.get_name().contains("::__olo_"));
                    assert_ne!(symbol.get_name(), "oneloopmaster::C0");
                }
                true
            });
        }
        assert!(nodes < 2_000);
        eprintln!("complete complex one-mass C0: {nodes} output nodes");
        let exact = Atom::evaluator_multiple(&bodies.coefficients().iter().map(Atom::as_view).collect::<Vec<_>>(), &parameters)
            .direct_translation(true).build().unwrap();
        let converter = Float::new(BITS);
        let mut compact = exact.map_coeff_with_prec(&|value| Complex::new(
            converter.from_rational(&value.re), converter.from_rational(&value.im)), BITS);
        let mut reference = PrecisionEvaluator::with_binary_precision_and_backend(
            ScalarIntegral::C0, BITS, EvaluationBackend::Expression).unwrap();
        let number = |text: &str| Float::parse(text, Some(BITS)).unwrap();
        let real = |value: Float| Complex::new(value, Float::new(BITS));
        let blank = || std::array::from_fn::<_, 3, _>(|_| real(number("0")));
        let absolute = |value: Float| if value.is_negative() { -value } else { value };
        let mut rows = 0;
        let mut failures = 0;
        let mut failed_inputs = std::collections::BTreeMap::new();
        for mass_re in ["-32", "-1", "-0.125", "0", "0.125", "1", "32"] {
            for width in ["0", "-1e-30", "-0.01", "-1", "-32"] {
                let mass = Complex::new(number(mass_re), number(width));
                let momentum_scale = if mass.re.is_zero() { number("1") } else { absolute(mass.re.clone()) };
                for ratio in ["-100", "-4", "-1.00001", "-1", "-0.99999", "-1e-30", "0", "1e-30", "0.1", "0.99999", "1", "1.00001", "4", "100"] {
                    let momentum = number(ratio) * momentum_scale.clone();
                    for scale in ["0.03125", "1", "1024"] {
                        let point = [real(momentum.clone()), mass.clone(), real(number(scale))];
                        let mut actual = blank();
                        compact.evaluate(&point, &mut actual);
                        for mass_index in 0..3 {
                            let mut input = vec![real(number("0")); 7];
                            input[(mass_index + 1) % 3] = point[0].clone();
                            input[mass_index + 3] = mass.clone();
                            input[6] = point[2].clone();
                            let mut expected = blank();
                            reference.evaluate(&input, &mut expected).unwrap();
                            for (index, (a, b)) in actual.iter().zip(&expected).enumerate() {
                                for (part, a, b) in [("real", &a.re, &b.re), ("imaginary", &a.im, &b.im)] {
                                    assert!(a.is_finite() && b.is_finite(), "nonfinite s={momentum}, m={mass}, mu={scale}, mass_index={mass_index}, coefficient={index} {part}: {a} vs {b}");
                                    let error = absolute(a.clone() - b.clone());
                                    let tolerance = number("1e-60") * (absolute(b.clone()) + number("1e-40"));
                                    if error >= tolerance {
                                        if failures < 8 {
                                            eprintln!("s={momentum}, m={mass}, mu={scale}, mass_index={mass_index}, coefficient={index} {part}: {a} vs {b}; error={error:e}, tolerance={tolerance:e}");
                                        }
                                        failures += 1;
                                        *failed_inputs.entry((mass_re, width)).or_insert(0) += 1;
                                    }
                                }
                            }
                            rows += 1;
                        }
                    }
                }
            }
        }
        eprintln!("complex/default-untyped one-mass C0: {rows} original mapped-reference rows at {BITS} bits, all three coefficients");
        eprintln!("failed input mass/width groups: {failed_inputs:?}");
        assert_eq!(failures, 0, "failed coefficient components");
    }).unwrap().join().unwrap();
}

#[test]
fn complete_complex_triangle_preserves_a_thousand_digits() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            const BITS: u32 = 3456;
            const REFERENCE_BITS: u32 = 3840;
            let s = symbolica::symbol!("inspect_thousand_s"; Real).to_atom();
            let m = symbolica::symbol!("inspect_thousand_m").to_atom();
            let mu = symbolica::symbol!("inspect_thousand_mu"; Positive).to_atom();
            let parameters = [s.clone(), m.clone(), mu.clone()];
            let zero = Atom::num(0);
            let bodies = get_expression(
                ScalarIntegral::C0,
                &[zero.clone(), zero.clone(), s, zero.clone(), m, zero, mu],
            )
            .unwrap();
            let converter = Float::new(BITS);
            let mut compact = Atom::evaluator_multiple(
                &bodies
                    .coefficients()
                    .iter()
                    .map(Atom::as_view)
                    .collect::<Vec<_>>(),
                &parameters,
            )
            .direct_translation(true)
            .build()
            .unwrap()
            .map_coeff_with_prec(
                &|value| {
                    Complex::new(
                        converter.from_rational(&value.re),
                        converter.from_rational(&value.im),
                    )
                },
                BITS,
            );
            let mut reference = PrecisionEvaluator::with_binary_precision_and_backend(
                ScalarIntegral::C0,
                REFERENCE_BITS,
                EvaluationBackend::Expression,
            )
            .unwrap();
            let complex = |re: &str, im: &str, bits| {
                Complex::new(
                    Float::parse(re, Some(bits)).unwrap(),
                    Float::parse(im, Some(bits)).unwrap(),
                )
            };
            let absolute = |value: Float| if value.is_negative() { -value } else { value };
            let tolerance = Float::parse("1e-1000", Some(REFERENCE_BITS)).unwrap();
            for (re, im) in [
                ("0", "-1"),
                ("0", "-1e-30"),
                ("1.25", "-1e-30"),
                ("-1.25", "-0.5"),
            ] {
                for momentum in ["-100", "0", "1"] {
                    let point = [
                        complex(momentum, "0", BITS),
                        complex(re, im, BITS),
                        complex("4.75", "0", BITS),
                    ];
                    let mut actual = std::array::from_fn::<_, 3, _>(|_| complex("0", "0", BITS));
                    compact.evaluate(&point, &mut actual);
                    let input = [
                        complex("0", "0", REFERENCE_BITS),
                        complex("0", "0", REFERENCE_BITS),
                        complex(momentum, "0", REFERENCE_BITS),
                        complex("0", "0", REFERENCE_BITS),
                        complex(re, im, REFERENCE_BITS),
                        complex("0", "0", REFERENCE_BITS),
                        complex("4.75", "0", REFERENCE_BITS),
                    ];
                    let mut expected =
                        std::array::from_fn::<_, 3, _>(|_| complex("0", "0", REFERENCE_BITS));
                    reference.evaluate(&input, &mut expected).unwrap();
                    for (a, b) in actual.iter().zip(expected.iter()) {
                        for (a, b) in [(&a.re, &b.re), (&a.im, &b.im)] {
                            assert!(a.is_finite() && b.is_finite());
                            let error = absolute(a.clone() - b.clone());
                            let bound = tolerance.clone()
                                * (Float::with_val(REFERENCE_BITS, 1) + absolute(b.clone()));
                            assert!(
                                error < bound,
                                "s={momentum},m=({re},{im}): error={error:e}, bound={bound:e}"
                            );
                        }
                    }
                }
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
