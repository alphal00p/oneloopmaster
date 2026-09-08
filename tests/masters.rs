//! Public master Symbols: direct hooks, transparent maps, and original Fortran data.
use oneloop::{A0, B0, C0, D0, OneLoopExpressions, dB0};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    domains::{
        float::{Complex, Float, SingleFloat},
        rational::Rational,
    },
    evaluate::ExpressionEvaluator,
    parse,
    prelude::{FloatLike, Real, RealLike},
};

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

fn masters() -> [Symbol; 5] {
    [A0(), B0(), dB0(), C0(), D0()]
}

fn calls(symbol: Symbol, args: &[Atom]) -> [Atom; 3] {
    [0, -1, -2].map(|tag| {
        let mut values = vec![Atom::num(tag)];
        values.extend_from_slice(args);
        symbol.call(&values)
    })
}

fn compile_calls(expressions: &[Atom; 3], args: &[Atom]) -> ExpressionEvaluator<Complex<Rational>> {
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    Atom::evaluator_multiple(&views, args)
        .direct_translation(true)
        .build()
        .unwrap()
}

fn double(evaluator: ExpressionEvaluator<Complex<Rational>>) -> ExpressionEvaluator<Complex<f64>> {
    evaluator.map_coeff(&|value| Complex::new(value.re.to_f64(), value.im.to_f64()))
}

fn arbitrary(
    evaluator: ExpressionEvaluator<Complex<Rational>>,
    bits: u32,
) -> ExpressionEvaluator<Complex<Float>> {
    let converter = Float::new(bits);
    evaluator.map_coeff_with_prec(
        &|value| {
            Complex::new(
                converter.from_rational(&value.re),
                converter.from_rational(&value.im),
            )
        },
        bits,
    )
}

fn comparison_failure(
    actual: Complex<f64>,
    expected: Complex<f64>,
    normalization: f64,
    label: &str,
) -> Option<String> {
    let error = (actual.re - expected.re).hypot(actual.im - expected.im) * normalization;
    let tolerance = 2e-10 + 2e-8 * expected.re.hypot(expected.im) * normalization;
    if actual.re.is_finite() && actual.im.is_finite() && error.is_finite() && error <= tolerance {
        None
    } else {
        Some(format!(
            "{label}: actual={actual:?}, expected={expected:?}, normalized error={error}"
        ))
    }
}

fn assert_close(actual: Complex<f64>, expected: Complex<f64>, normalization: f64, label: &str) {
    if let Some(error) = comparison_failure(actual, expected, normalization, label) {
        panic!("{error}");
    }
}

#[test]
fn public_master_hooks_and_maps_match_every_acceptance_fixture() {
    on_stack(|| {
        let symbols = masters();
        let context = OneLoopExpressions::new();
        // Repeated initialization must be harmless; these raw full-name Symbols
        // must be the same registered objects, including their evaluation hooks.
        let _again = OneLoopExpressions::new();
        let raw = [
            symbolica::symbol!("oneloopmaster::A0"),
            symbolica::symbol!("oneloopmaster::B0"),
            symbolica::symbol!("oneloopmaster::dB0"),
            symbolica::symbol!("oneloopmaster::C0"),
            symbolica::symbol!("oneloopmaster::D0"),
        ];
        assert_eq!(symbols, raw);
        assert_eq!(masters(), raw);
        assert_eq!(parse!("oneloopmaster::A0(-1, 2, 1)"), A0().call((-1, 2, 1)));

        let mut direct = Vec::new();
        let mut mapped = Vec::new();
        for (family, (&symbol, arity)) in symbols.iter().zip([2, 4, 4, 7, 11]).enumerate() {
            let args = (0..arity)
                .map(|index| {
                    symbolica::symbol!(format!("master_fixture_{family}_{index}")).to_atom()
                })
                .collect::<Vec<_>>();
            let coefficients = calls(symbol, &args);
            direct.push(double(compile_calls(&coefficients, &args)));
            mapped.push(double(context.evaluator(&coefficients, &args).unwrap()));
        }

        let mut inventory = [0; 5];
        let mut pole_inventory = [0; 5];
        let mut failures = Vec::new();
        for (line_number, line) in include_str!("data/parity.txt").lines().enumerate() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let values = line
                .split_whitespace()
                .map(|value| value.parse::<f64>().unwrap())
                .collect::<Vec<_>>();
            let (family, masses, momenta, power) = match values[0] {
                1. => (0, 1, 0, -1),
                2. => (1, 2, 1, 0),
                -2. => (2, 2, 1, 1),
                3. => (3, 3, 3, 1),
                4. => (4, 4, 6, 2),
                other => panic!("invalid family {other}"),
            };
            assert_eq!(values.len(), 2 + momenta + 2 * masses + 6);
            let mass_start = 2 + momenta;
            let reference_start = mass_start + 2 * masses;
            let mut args = values[2..mass_start]
                .iter()
                .map(|&value| Complex::new(value, 0.))
                .collect::<Vec<_>>();
            args.extend(
                values[mass_start..reference_start]
                    .chunks_exact(2)
                    .map(|value| Complex::new(value[0], value[1])),
            );
            args.push(Complex::new(values[1], 0.));
            let normalization = values[1].powi(power);
            let mut hook_output = [Complex::new(0., 0.); 3];
            let mut native_output = [Complex::new(0., 0.); 3];
            direct[family].evaluate(&args, &mut hook_output);
            mapped[family].evaluate(&args, &mut native_output);
            for (tag, expected) in values[reference_start..].chunks_exact(2).enumerate() {
                let expected = Complex::new(expected[0], expected[1]);
                if let Some(error) = comparison_failure(
                    hook_output[tag],
                    expected,
                    normalization,
                    &format!("direct line {} coefficient {tag}", line_number + 1),
                ) {
                    failures.push(error);
                }
                if let Some(error) = comparison_failure(
                    native_output[tag],
                    expected,
                    normalization,
                    &format!("mapped line {} coefficient {tag}", line_number + 1),
                ) {
                    failures.push(error);
                }
                if let Some(error) = comparison_failure(
                    hook_output[tag],
                    native_output[tag],
                    normalization,
                    "hook versus native",
                ) {
                    failures.push(error);
                }
            }
            inventory[family] += 1;
            if values[reference_start + 2..]
                .iter()
                .any(|&value| value != 0.)
            {
                pole_inventory[family] += 1;
            }
        }
        assert_eq!(inventory, [4, 13, 12, 125, 139]);
        assert!(pole_inventory.into_iter().all(|count| count > 0));
        assert!(
            failures.is_empty(),
            "{} master comparisons failed:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}

#[test]
fn constant_and_variable_master_calls_preserve_arbitrary_precision() {
    on_stack(|| {
        let context = OneLoopExpressions::new();
        let symbols = masters();
        // Elementary controls include the regular derivative pseudothreshold
        // and exact massive vacuum limits. Every tag is checked independently.
        let arguments = [
            vec![-2, 1],
            vec![0, 3, 3, 5],
            vec![1, 1, 4, 1],
            vec![0, 0, 0, 1, 1, 1, 1],
            vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        ];
        // Compile each native graph once, then remap its exact constants at
        // each precision. Rebuilding C0/D0 for each bit count is unnecessary.
        let compiled = (0..5)
            .map(|family| {
                let constants = arguments[family]
                    .iter()
                    .map(|&v| Atom::num(v))
                    .collect::<Vec<_>>();
                let parameters = (0..constants.len())
                    .map(|index| {
                        symbolica::symbol!(format!("precision_master_{family}_{index}")).to_atom()
                    })
                    .collect::<Vec<_>>();
                let constant_calls = calls(symbols[family], &constants);
                let variable_calls = calls(symbols[family], &parameters);
                eprintln!("compiling master precision controls for family {family}");
                [
                    compile_calls(&constant_calls, &[]),
                    compile_calls(&variable_calls, &parameters),
                    context.evaluator(&variable_calls, &parameters).unwrap(),
                ]
            })
            .collect::<Vec<_>>();
        for bits in [128, 256] {
            let zero = || Float::new(bits);
            let number = |value| Float::with_val(bits, value);
            let complex = |re, im| Complex::new(re, im);
            let log_two = number(2).log();
            let pi = Float::with_val(bits, symbolica::domains::backend::float::Constant::Pi);
            let expected = [
                [
                    complex(number(-2) * (number(1) - log_two.clone()), number(-2) * pi),
                    complex(number(-2), zero()),
                    complex(zero(), zero()),
                ],
                [
                    complex(-(number(3) / number(5)).log(), zero()),
                    complex(number(1), zero()),
                    complex(zero(), zero()),
                ],
                [
                    complex(number(3) * log_two - number(2), zero()),
                    complex(zero(), zero()),
                    complex(zero(), zero()),
                ],
                [
                    complex(-number(1) / number(2), zero()),
                    complex(zero(), zero()),
                    complex(zero(), zero()),
                ],
                [
                    complex(number(1) / number(6), zero()),
                    complex(zero(), zero()),
                    complex(zero(), zero()),
                ],
            ];
            for family in 0..5 {
                eprintln!("checking master precision controls for family {family} at {bits} bits");
                let [mut constant_hook, mut variable_hook, mut variable_native] = compiled[family]
                    .clone()
                    .map(|evaluator| arbitrary(evaluator, bits));
                let inputs = arguments[family]
                    .iter()
                    .map(|&value| complex(number(value), zero()))
                    .collect::<Vec<_>>();
                let mut from_constant =
                    core::array::from_fn::<_, 3, _>(|_| complex(zero(), zero()));
                let mut from_variable = from_constant.clone();
                let mut from_native = from_constant.clone();
                constant_hook.evaluate(&[], &mut from_constant);
                variable_hook.evaluate(&inputs, &mut from_variable);
                variable_native.evaluate(&inputs, &mut from_native);
                let tolerance = if bits == 128 { 1e-32 } else { 1e-69 };
                for tag in 0..3 {
                    for (route, value) in [
                        ("constant hook", &from_constant[tag]),
                        ("variable hook", &from_variable[tag]),
                        ("native map", &from_native[tag]),
                    ] {
                        let difference = value.clone() - expected[family][tag].clone();
                        assert!(
                            difference.re.to_f64().hypot(difference.im.to_f64()) < tolerance,
                            "{route} family={family} tag={tag} bits={bits}: {value} != {}",
                            expected[family][tag]
                        );
                    }
                }
            }

            // Numeric arguments below f64 resolution must survive unchanged.
            // Give the imaginary component more precision than the real part
            // to ensure the callback considers both components and every arg.
            let x = parse!("master_tiny_width");
            let expression = A0().call((-1, &x, 1));
            let mut evaluator = arbitrary(
                compile_calls(&[expression.clone(), expression.clone(), expression], &[x]),
                bits,
            );
            let tiny = number(1) / number(2).pow(100);
            let mut output = core::array::from_fn::<_, 3, _>(|_| complex(zero(), zero()));
            evaluator.evaluate(
                &[complex(Float::with_val(64, -2), -tiny.clone())],
                &mut output,
            );
            assert_eq!(output[0].im, -tiny);
        }

        // Repeated and mixed calls exercise external-function indexing when
        // constant calls and runtime calls occur in the same evaluator.
        let x = parse!("mixed_master_x");
        let mixed = A0().call((0, 1, 1))
            + B0().call((-1, &x, 1, 1, 1))
            + C0().call((-2, 0, 0, -1, 0, 0, 0, 1))
            + D0().call((-2, 0, 0, 0, 0, -1, -2, 0, 0, 0, 0, 1));
        let mut evaluator = double(compile_calls(&[mixed.clone(), mixed.clone(), mixed], &[x]));
        let mut actual = [Complex::new(0., 0.); 3];
        evaluator.evaluate(&[Complex::new(-1., 0.)], &mut actual);
        // C0(-1)/epsilon²=-1 and the massless box has 4/(s*t)=2.
        assert_close(
            actual[0],
            Complex::new(3., 0.),
            1.,
            "mixed master expression",
        );
    });
}

#[test]
fn malformed_master_tags_and_arities_are_deliberate_errors() {
    on_stack(|| {
        for (symbol, arity) in masters().into_iter().zip([2, 4, 4, 7, 11]) {
            let info = symbol.get_evaluation_info().expect("master EvaluationInfo");
            assert!(
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(
                    || info.get_evaluator::<Complex<f64>>(&[])
                ))
                .is_err(),
                "missing Laurent tag must fail"
            );
            for tag in [
                parse!("1"),
                parse!("-3"),
                parse!("1/2"),
                parse!("invalid_master_tag"),
            ] {
                let failure = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    info.get_evaluator::<Complex<f64>>(&[tag.as_view()])
                }));
                let panic = failure.err().expect("invalid Laurent tag must fail");
                let message = panic.downcast_ref::<String>().unwrap();
                assert!(
                    message.contains("Laurent tag") && message.contains("expected 0, -1, or -2")
                );
            }
            let tag = Atom::num(0);
            let evaluator = info
                .get_evaluator::<Complex<f64>>(&[tag.as_view()])
                .unwrap();
            for count in [0, arity - 1, arity + 1] {
                let args = vec![Complex::new(1., 0.); count];
                let failure =
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| evaluator(&args)));
                let panic = failure.expect_err("invalid numeric arity must fail");
                let message = panic.downcast_ref::<String>().unwrap();
                assert!(
                    message.contains("numeric arguments after its Laurent tag"),
                    "{message}"
                );
                // Also exercise Symbolica's actual constant-call dispatch:
                // a tag-only call uses its precision-aware constant route.
                let constants = vec![Atom::num(1); count];
                let expressions = calls(symbol, &constants);
                let failure = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    let mut evaluator = double(compile_calls(&expressions, &[]));
                    let mut output = [Complex::new(0., 0.); 3];
                    evaluator.evaluate(&[], &mut output);
                }));
                let panic = failure.expect_err("malformed constant master must fail");
                let message = panic.downcast_ref::<String>().unwrap();
                assert!(
                    message.contains("numeric arguments after its Laurent tag"),
                    "{message}"
                );
            }
        }
    });
}
