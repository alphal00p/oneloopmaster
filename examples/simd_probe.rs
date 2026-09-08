//! Bounded correctness probe for SymJIT's scalar-layout complex batch API.
//!
//! Run only while holding the single Symbolica runtime slot. The optional
//! first argument selects one case; the default runs all cases. No timing
//! results from this diagnostic should be treated as performance claims.
use symbolica::evaluate::{FunctionRegistrationOptions, InliningPolicy};
use symbolica::prelude::*;

const CASES: &[&str] = &[
    "arithmetic",
    "truth",
    "polylog",
    "top_if",
    "never_if",
    "nested_if",
    "never_polylog",
    "never_if_polylog",
];

fn main() {
    let selected = std::env::args().nth(1).unwrap_or_else(|| "all".into());
    assert!(selected == "all" || CASES.contains(&selected.as_str()));
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(move || run(&selected))
        .unwrap()
        .join()
        .unwrap();
}

fn run(selected: &str) {
    let mut failures = Vec::new();
    for &case in CASES {
        if selected != "all" && selected != case {
            continue;
        }
        check_case(case, &mut failures);
    }
    assert!(
        failures.is_empty(),
        "{} SIMD probe failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn check_case(case: &str, failures: &mut Vec<String>) {
    let condition_symbol = symbol!(format!("simd_probe_{case}_condition"));
    let z_symbol = symbol!(format!("simd_probe_{case}_z"));
    let condition = condition_symbol.to_atom();
    let z = z_symbol.to_atom();
    let params = [condition.clone(), z.clone()];
    let conditional = Symbol::IF.call((&condition, &z / &condition, &z + 1));
    let mut map = FunctionMap::new();
    let never = || FunctionRegistrationOptions::new().inlining(InliningPolicy::Never);
    let inner = symbol!(format!("simd_probe_{case}_inner"));
    let body = match case {
        "arithmetic" => &z * &z + &condition,
        "truth" => Symbol::IF.call((&condition, Atom::num(7), Atom::num(-3))),
        "polylog" => z.polylog(2),
        "top_if" => conditional,
        "never_if" | "nested_if" | "never_polylog" | "never_if_polylog" => {
            let inner_body = match case {
                "never_polylog" => z.polylog(2),
                "never_if_polylog" => Symbol::IF.call((&condition, z.polylog(2), &z + 1)),
                _ => conditional,
            };
            map.add_function_with_options(
                inner,
                vec![condition_symbol, z_symbol],
                inner_body,
                never(),
            )
            .unwrap();
            let inner_call = inner.call((&condition, &z));
            if case == "nested_if" {
                let outer = symbol!(format!("simd_probe_{case}_outer"));
                map.add_function_with_options(
                    outer,
                    vec![condition_symbol, z_symbol],
                    inner_call + &z,
                    never(),
                )
                .unwrap();
                outer.call((&condition, &z))
            } else {
                inner_call
            }
        }
        _ => unreachable!(),
    };
    // Three outputs exercise the layout required by Laurent coefficients.
    let output_atoms = [
        body.clone(),
        Atom::num(2) * &body + &z,
        body.conj() - &condition,
    ];
    let outputs = output_atoms.each_ref().map(Atom::as_view);
    let exact = Atom::evaluator_multiple(&outputs, &params)
        .function_map(map)
        .direct_translation(true)
        .build()
        .unwrap();
    let mut native = exact
        .clone()
        .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    // Cover every component-zero pattern and nonzero quadrant. In particular,
    // purely imaginary conditions are true even though their real part is zero.
    let mut rows = [
        ((0., 0.), (0.2, 0.7)),
        ((0., 1.), (-0.4, 0.3)),
        ((0., -1.), (0.6, -0.2)),
        ((1., 0.), (1.4, -0.6)),
        ((-2., 0.), (-1.1, -0.8)),
        ((2., 3.), (0.05, 0.9)),
        ((2., -3.), (2.2, 0.4)),
        ((-2., 3.), (-0.7, 1.3)),
        ((-2., -3.), (0.3, -0.8)),
        ((-0., 0.), (-0.5, 0.6)),
        ((0., -0.), (1.7, -0.4)),
        ((-0., -0.), (-0.9, 0.2)),
    ]
    .map(|((cre, cim), (re, im))| [Complex::new(cre, cim), Complex::new(re, im)])
    .to_vec();
    if case == "truth" {
        // No division by the condition in this case: exact truth must also
        // survive subnormal inputs, without squaring or an epsilon threshold.
        let tiny = f64::from_bits(1);
        rows.extend([
            [Complex::new(0., tiny), Complex::new(0.2, 0.7)],
            [Complex::new(0., -tiny), Complex::new(-0.4, 0.3)],
            [Complex::new(tiny, 0.), Complex::new(0.6, -0.2)],
            [Complex::new(-tiny, tiny), Complex::new(1.4, -0.6)],
        ]);
    }
    let mut expected = vec![[Complex::new(0., 0.); 3]; rows.len()];
    for (input, out) in rows.iter().zip(&mut expected) {
        native.evaluate(input, out);
        if case == "truth" {
            let truth = input[0].re != 0. || input[0].im != 0.;
            assert_eq!(out[0], Complex::new(if truth { 7. } else { -3. }, 0.));
        }
    }

    for (mode, simd, simd_branch) in [
        ("scalar", false, false),
        ("simd_fallback", true, false),
        ("simd_branches", true, true),
    ] {
        let settings = JITCompilationSettings::new()
            .optimization_level(2)
            .direct_translation(true)
            .with_option("use_simd", simd.to_string())
            .with_option("simd_branch", simd_branch.to_string())
            .with_option("enable_simd512", "false")
            .with_option("use_threads", "false")
            .with_option("fastmath", "false")
            .with_option("fast_complex", "false");
        let mut compiled = match exact.jit_compile::<Complex<f64>>(settings) {
            Ok(compiled) => compiled,
            Err(error) => {
                failures.push(format!("{case}/{mode}: compile error: {error}"));
                continue;
            }
        };
        let before = failures.len();
        for (index, input) in rows.iter().enumerate() {
            let mut actual = [Complex::new(f64::NAN, f64::NAN); 3];
            compiled.evaluate(input, &mut actual);
            compare(
                case,
                mode,
                "point",
                index,
                &actual,
                &expected[index],
                failures,
            );
        }
        for size in [3, 4, 5, 8, 12, 16] {
            if size > rows.len() {
                continue;
            }
            let input = rows[..size].iter().flatten().copied().collect::<Vec<_>>();
            let mut actual = vec![Complex::new(f64::NAN, f64::NAN); 3 * size];
            compiled.batch_evaluate(&input, &mut actual, size);
            for (index, result) in actual.chunks_exact(3).enumerate() {
                compare(
                    case,
                    mode,
                    &format!("batch{size}"),
                    index,
                    result,
                    &expected[index],
                    failures,
                );
            }
        }
        eprintln!("{case}/{mode}: {} mismatches", failures.len() - before);
    }
}

fn compare(
    case: &str,
    mode: &str,
    batch: &str,
    row: usize,
    actual: &[Complex<f64>],
    expected: &[Complex<f64>],
    failures: &mut Vec<String>,
) {
    for (column, (a, b)) in actual.iter().zip(expected).enumerate() {
        if case == "truth" && column == 0 {
            if a != b {
                failures.push(format!(
                    "{case}/{mode}/{batch} row{row}: exact truth {a:?}, expected {b:?}"
                ));
            }
            continue;
        }
        let error = (a.re - b.re).hypot(a.im - b.im);
        let tolerance = 2e-12 * (1. + b.re.hypot(b.im));
        if !error.is_finite() || error > tolerance {
            failures.push(format!(
                "{case}/{mode}/{batch} row{row} col{column}: {a:?}, expected {b:?}"
            ));
        }
    }
}
