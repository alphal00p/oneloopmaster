//! Finite boxes with a pair of massive lines at its pseudo-threshold, in all
//! 24 labellings of their propagators.
//!
//! At the pseudo-threshold `p = (√mi - √mj)²` of the lines `i` and `j` that
//! the invariant `p` joins, their Källén function `(mi + mj - p)² - 4 mi mj`
//! vanishes, and their r-function `r + 1/r = k = (mi + mj - p) / (√mi √mj)`
//! has the double root `r = 1`. With masses without exact square roots `k`
//! rounds to either side of 2, depending on the labelling. That split the
//! root into `1 ± sqrt(rounding)`, which cost `Float` half its digits, and
//! the four-mass chart took the pair or not by the rounding of `k² - 4`: the
//! expression routes (the binary64 expression backend, SymJIT and the
//! interpreted `Float`) returned values up to 126 times their size off in some
//! labellings (alphal00p/oneloopmaster#30, review). `box_integral::chart_pair`
//! takes `k = 2` exactly and keeps the pair out of the four-mass chart's `p23`
//! and `p12` slots, and SymJIT takes these four-mass points from the native
//! kernels.
//!
//! Cases: dyadic squared masses `a² m` with `m` not a square, so that the
//! pseudo-thresholds are exact in binary64; four massive lines with a leg, a
//! channel, or both at their pseudo-threshold (the review's box among them),
//! in the Euclidean and the physical region.
//! References: the mean of the invariants at their pseudo-thresholds moved by
//! `±2^-300` at 1024 bits (`±2^-200` agrees to 1e-122 of the value). The
//! native `Float` route evaluates the generated code, which is not
//! regenerated and keeps about half its digits here; it is not checked. Nor
//! are three- and two-mass boxes and the contour box, which keep the rounded
//! `k` and in `Float` about half the digits.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, DoubleFloat, Float, Real};

const BITS: u32 = 256;

struct Case {
    name: &'static str,
    /// [p1², p2², p3², p4², s12, s23] of the labelling 0123.
    invariants: [f64; 6],
    masses: [f64; 4],
    mu_squared: f64,
    /// The finite part as [re, im] (the poles vanish).
    finite: [&'static str; 2],
}

const CASES: [Case; 6] = [
    Case {
        // Lines 0 and 1 (p1) at their pseudo-threshold, Euclidean: the
        // expression backend and the interpreted Float were 3.8 off in three
        // labellings.
        name: "four-masses-leg",
        invariants: [15.72265625, -23.9375, -19.375, -54.5, -58.0625, -58.625],
        masses: [22.640625, 0.62890625, 1.4150390625, 15.72265625],
        mu_squared: 4.0,
        finite: ["6.86268763297036602047657901858474590E-04", "0"],
    },
    Case {
        // The same, where the interpreted Float kept 2e-17.
        name: "four-masses-leg-half-digits",
        invariants: [4.69140625, -41.375, -0.625, -29.9375, -34.75, -55.9375],
        masses: [29.3212890625, 10.5556640625, 4.69140625, 10.5556640625],
        mu_squared: 4.0,
        finite: ["4.17938763643721600473979200984578736E-04", "0"],
    },
    Case {
        // Lines 1 and 3 (s23): the interpreted Float was 15 off in four
        // labellings.
        name: "four-masses-channel",
        invariants: [-29.6875, -20.3125, -19.6875, 13.25, -11.0, 24.85546875],
        masses: [99.421875, 44.1875, 44.1875, 2.76171875],
        mu_squared: 4.0,
        finite: ["1.17520749245973525720956286463999677E-04", "0"],
    },
    Case {
        // Lines 0 and 1 (p1) and lines 0 and 2 (s12): the expression backend
        // and SymJIT were 32 off in four labellings.
        name: "four-masses-leg-and-channel",
        invariants: [0.76171875, 9.625, -15.4375, -8.4375, 0.76171875, 2.125],
        masses: [19.04296875, 12.1875, 12.1875, 4.7607421875],
        mu_squared: 4.0,
        finite: ["1.15650708262821876229351189729039608E-03", "0"],
    },
    Case {
        // The same, Euclidean: the interpreted Float was 126 off in four
        // labellings.
        name: "four-masses-leg-and-channel-euclidean",
        invariants: [
            4.66796875, -34.8125, -41.5625, -28.875, 4.66796875, -51.3125,
        ],
        masses: [168.046875, 116.69921875, 116.69921875, 42.01171875],
        mu_squared: 4.0,
        finite: ["1.29830553447542804732251921412684490E-05", "0"],
    },
    Case {
        // The review's box: lines 0 and 1 (p1), the interpreted Float 1.4
        // off in four labellings.
        name: "four-masses-review",
        invariants: [
            1.0048828125,
            -52.609375,
            -2.984375,
            -35.5,
            -47.265625,
            3.65625,
        ],
        masses: [0.5126953125, 2.953125, 0.5126953125, 2.05078125],
        mu_squared: 4.0,
        finite: ["4.17487846559143050477294811528511928E-03", "0"],
    },
];

fn labellings() -> Vec<[usize; 4]> {
    let mut result = Vec::new();
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    if a != b && a != c && a != d && b != c && b != d && c != d {
                        result.push([a, b, c, d]);
                    }
                }
            }
        }
    }
    assert_eq!(result.len(), 24);
    result
}

fn relabel(case: &Case, [a, b, c, d]: [usize; 4]) -> Vec<f64> {
    let slot = |i: usize, k: usize| match (i.min(k), i.max(k)) {
        (0, 1) => 0,
        (1, 2) => 1,
        (2, 3) => 2,
        (0, 3) => 3,
        (0, 2) => 4,
        (1, 3) => 5,
        _ => unreachable!(),
    };
    let p = |i, k| case.invariants[slot(i, k)];
    let m = case.masses;
    vec![
        p(a, b),
        p(b, c),
        p(c, d),
        p(d, a),
        p(a, c),
        p(b, d),
        m[a],
        m[b],
        m[c],
        m[d],
        case.mu_squared,
    ]
}

/// Relative error of the finite part, and the poles against its size.
fn check(
    values: [[Float; 2]; 3],
    case: &Case,
    tolerance: f64,
    label: &str,
    failures: &mut Vec<String>,
) {
    let parse = |s: &str| Float::parse(s, Some(BITS)).unwrap();
    let expected = [parse(case.finite[0]), parse(case.finite[1])];
    let norm = |z: &[Float; 2]| (z[0].clone() * &z[0] + z[1].clone() * &z[1]).sqrt();
    let size = norm(&expected);
    let zero = Float::new(BITS);
    for (coefficient, value) in values.iter().enumerate() {
        let target = if coefficient == 0 {
            &expected
        } else {
            &[zero.clone(), zero.clone()]
        };
        let error = norm(&[value[0].clone() - &target[0], value[1].clone() - &target[1]]);
        if !error.is_finite() || error > size.clone() * Float::with_val(BITS, tolerance) {
            failures.push(format!(
                "{label}, coefficient {coefficient}: {:e} {:e} (relative error {:e})",
                value[0],
                value[1],
                error / &size
            ));
        }
    }
}

fn widen(values: &[Complex<f64>; 3]) -> [[Float; 2]; 3] {
    values.map(|z| [Float::with_val(BITS, z.re), Float::with_val(BITS, z.im)])
}

#[test]
fn binary64_routes_agree_in_every_labelling() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let mut handwritten = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
            let mut generated = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
            generated.set_handwritten(false);
            let mut failures = Vec::new();
            for case in &CASES {
                for order in labellings() {
                    let input: Vec<Complex<f64>> = relabel(case, order)
                        .into_iter()
                        .map(|v| Complex::new(v, 0.0))
                        .collect();
                    for route in ["handwritten", "generated", "expression", "symjit"] {
                        let mut output = [Complex::new(f64::NAN, f64::NAN); 3];
                        match route {
                            "handwritten" => handwritten.evaluate(&input, &mut output).unwrap(),
                            "generated" => generated.evaluate(&input, &mut output).unwrap(),
                            "expression" => evaluate_with_backend(
                                ScalarIntegral::D0,
                                &input,
                                &mut output,
                                EvaluationBackend::Expression,
                            )
                            .unwrap(),
                            _ => evaluate_with_backend(
                                ScalarIntegral::D0,
                                &input,
                                &mut output,
                                EvaluationBackend::SymJit,
                            )
                            .unwrap(),
                        }
                        let label = format!("{}, {route}, labelling {order:?}", case.name);
                        check(widen(&output), case, 1e-12, &label, &mut failures);
                    }
                }
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn double_double_agrees_in_every_labelling() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let mut native =
                NativeEvaluator::<DoubleFloat>::with_binary_precision(ScalarIntegral::D0, 106)
                    .unwrap();
            let mut failures = Vec::new();
            for case in &CASES {
                for order in labellings() {
                    let input: Vec<Complex<DoubleFloat>> = relabel(case, order)
                        .into_iter()
                        .map(|v| Complex::new(DoubleFloat::from(v), DoubleFloat::from(0.0)))
                        .collect();
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(DoubleFloat::from(0.0), DoubleFloat::from(0.0))
                    });
                    native.evaluate(&input, &mut output).unwrap();
                    let wide = |x: &DoubleFloat| {
                        let parts = x.into_inner();
                        Float::with_val(BITS, parts.hi()) + &Float::with_val(BITS, parts.lo())
                    };
                    let values = output.map(|z| [wide(&z.re), wide(&z.im)]);
                    let label = format!("{}, DoubleFloat, labelling {order:?}", case.name);
                    check(values, case, 1e-26, &label, &mut failures);
                }
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn interpreted_float_agrees_in_every_labelling() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let bits = 128;
            let mut evaluator = PrecisionEvaluator::with_binary_precision_and_backend(
                ScalarIntegral::D0,
                bits,
                EvaluationBackend::Expression,
            )
            .unwrap();
            let mut failures = Vec::new();
            for case in &CASES {
                for order in labellings() {
                    let input: Vec<Complex<Float>> = relabel(case, order)
                        .into_iter()
                        .map(|v| Complex::new(Float::with_val(bits, v), Float::new(bits)))
                        .collect();
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(Float::new(bits), Float::new(bits))
                    });
                    evaluator.evaluate(&input, &mut output).unwrap();
                    let label = format!(
                        "{}, interpreted Float {bits}, labelling {order:?}",
                        case.name
                    );
                    check(
                        output.map(|z| [z.re, z.im]),
                        case,
                        1e-30,
                        &label,
                        &mut failures,
                    );
                }
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}
