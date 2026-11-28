//! Boxes with three massive propagators pairwise at their pseudo-thresholds
//! (a pseudo-threshold triangle), in all 24 labellings of their propagators
//! and on every route.
//!
//! With `p_jk = (√m_j - √m_k)²` for the three pairs of the lines `j, k` of the
//! triangle (equal masses at `p_jk = 0` among them: two collinear massless legs
//! between equal-mass lines), the modified Cayley matrix has a rank-one block
//! and `det Y = 0`. Every finite form of the sector code divided by a power of
//! it, and returned values up to `1e77` off, also at 1024 bits (MadGraph's
//! stress points of `u u~ > z g g` and `u u~ > z z g g`: two exactly collinear
//! gluons on a b-quark loop). `native::pseudo_triangle` evaluates them in
//! closed form, `D0 = κ / (v_1 v_2 v_3) · Δ²[L B0(L; m_o, κ²)](L_1, L_2, L_3)`,
//! on every route.
//!
//! Cases: line 0 outside the triangle (massive, massless or complex), its
//! three invariants to the triangle distinct, two or three of their `L`
//! equal, in the Euclidean and the physical region; and all six pairs at their
//! pseudo-thresholds (`D0 = 1/(6 √(m_0 m_1 m_2 m_3))`).
//! References: the closed form at 130 digits (`mpmath`, the bubbles by
//! quadrature, split at the real parts of the roots of their denominators;
//! equal `L` spread by `±1e-35`). Where the box has real masses or Euclidean
//! invariants they agree with the direct integral over the simplex (to its 15
//! digits), and with the limit of the 1024-bit sector code at the invariants
//! moved by `±2^-300` (to 20 digits or more) where that converges.
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
    /// The squared masses as [re, im].
    masses: [[f64; 2]; 4],
    /// The finite part as [re, im] (the poles vanish).
    finite: [&'static str; 2],
}

const B: [f64; 2] = [22.090000000000003, 0.0];

const CASES: [Case; 12] = [
    Case {
        // MadGraph's stress point of u u~ > z g g (call 7050): two exactly
        // collinear gluons on a b-quark loop. Every route was 3e10 or more
        // off on 4499b48 (at 1024 bits, 1e77).
        name: "collinear-gluons",
        invariants: [
            8315.251344000004,
            0.0,
            0.0,
            8315.253007050354,
            8315.252165165437,
            0.0,
        ],
        masses: [B; 4],
        finite: [
            "-2.73622268832694949369560741595724296e-6",
            "-2.45306321453653302399423297211483412e-10",
        ],
    },
    Case {
        // The same with p1² = p4²: two of the L coincide.
        name: "collinear-equal-legs",
        invariants: [8315.251344, 0.0, 0.0, 8315.251344, 8315.252165165437, 0.0],
        masses: [B; 4],
        finite: [
            "-2.73622287166984472639326260897536104e-6",
            "-2.45306370778388725296289234287827326e-10",
        ],
    },
    Case {
        name: "equal-masses-physical",
        invariants: [13.5, 0.0, 0.0, 17.0, 11.25, 0.0],
        masses: [[2.0, 0.0], [1.5, 0.0], [1.5, 0.0], [1.5, 0.0]],
        finite: [
            "-2.87910101735170706877135848640733283e-2",
            "-1.46767510946563685335343697898402937e-2",
        ],
    },
    Case {
        name: "equal-masses-mixed",
        invariants: [-3.5, 0.0, 0.0, 17.0, 11.25, 0.0],
        masses: [[2.0, 0.0], [1.5, 0.0], [1.5, 0.0], [1.5, 0.0]],
        finite: [
            "-5.37618282237617437448883440956121396e-2",
            "6.06365054265331967071226889961335908e-2",
        ],
    },
    Case {
        // A massless line 0: a three-mass box.
        name: "three-masses",
        invariants: [-3.5, 0.0, 0.0, -7.0, -1.25, 0.0],
        masses: [[0.0, 0.0], [1.5, 0.0], [1.5, 0.0], [1.5, 0.0]],
        finite: ["6.46932106819954838218238676939865701e-2", "0"],
    },
    Case {
        name: "three-masses-physical",
        invariants: [13.5, 0.0, 0.0, 17.0, 11.25, 0.0],
        masses: [[0.0, 0.0], [1.5, 0.0], [1.5, 0.0], [1.5, 0.0]],
        finite: ["-2.70896982743376596120998108128684047e-2", "0"],
    },
    Case {
        // p1² on shell: L at the branch point of the bubble (the box is
        // finite with one on-shell leg next to the massless line).
        name: "three-masses-leg-on-shell",
        invariants: [1.5, 0.0, 0.0, -7.0, -1.25, 0.0],
        masses: [[0.0, 0.0], [1.5, 0.0], [1.5, 0.0], [1.5, 0.0]],
        finite: ["1.3083655093539604007012442450622305e-1", "0"],
    },
    Case {
        // Unequal masses at p_jk = (√m_j - √m_k)².
        name: "unequal-masses",
        invariants: [-3.5, 1.0, 1.0, -7.0, -1.25, 4.0],
        masses: [[2.0, 0.0], [1.0, 0.0], [4.0, 0.0], [9.0, 0.0]],
        finite: ["1.27641909943134113055606706647040695e-2", "0"],
    },
    Case {
        name: "unequal-masses-physical",
        invariants: [23.5, 1.0, 1.0, 27.0, 21.25, 4.0],
        masses: [[2.0, 0.0], [1.0, 0.0], [4.0, 0.0], [9.0, 0.0]],
        finite: [
            "-7.37912640766910510761457250514854624e-3",
            "-3.10110675079130542652453451136217279e-3",
        ],
    },
    Case {
        // All three L equal (L = -3 for κ = 4).
        name: "unequal-masses-coincident",
        invariants: [-2.25, 1.0, 4.0, -3.0, -4.5, 9.0],
        masses: [[2.0, 0.0], [1.0, 0.0], [4.0, 0.0], [16.0, 0.0]],
        finite: ["9.66754557217061803866647526251602181e-3", "0"],
    },
    Case {
        // A complex line 0 (a top quark with its width) on the b-quark triangle.
        name: "complex-mass",
        invariants: [36000.5, 0.0, 0.0, 33000.25, 35000.75, 0.0],
        masses: [[30032.890000000003, -259.375], B, B, B],
        finite: [
            "-5.02819545839569506834295770354649413e-6",
            "-2.49076603017632943047448445700643417e-6",
        ],
    },
    Case {
        // All six pairs at their pseudo-thresholds: 1/(6 · 1 · 2 · 3 · 5).
        name: "all-six-pairs",
        invariants: [1.0, 1.0, 4.0, 16.0, 4.0, 9.0],
        masses: [[1.0, 0.0], [4.0, 0.0], [9.0, 0.0], [25.0, 0.0]],
        finite: ["5.55555555555555555555555555555555556e-3", "0"],
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

/// The arguments `[re, im]` of the labelling `[a, b, c, d]`.
fn relabel(case: &Case, [a, b, c, d]: [usize; 4]) -> Vec<[f64; 2]> {
    let slot = |i: usize, k: usize| match (i.min(k), i.max(k)) {
        (0, 1) => 0,
        (1, 2) => 1,
        (2, 3) => 2,
        (0, 3) => 3,
        (0, 2) => 4,
        (1, 3) => 5,
        _ => unreachable!(),
    };
    let p = |i, k| [case.invariants[slot(i, k)], 0.0];
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
        [1.0, 0.0],
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

fn spawn(body: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(body)
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn binary64_routes_agree_in_every_labelling() {
    spawn(|| {
        let mut handwritten = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
        let mut generated = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
        generated.set_handwritten(false);
        let mut failures = Vec::new();
        for case in &CASES {
            for order in labellings() {
                let input: Vec<Complex<f64>> = relabel(case, order)
                    .into_iter()
                    .map(|[re, im]| Complex::new(re, im))
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
                    check(widen(&output), case, 1e-15, &label, &mut failures);
                }
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    });
}

#[test]
fn double_double_agrees_in_every_labelling() {
    spawn(|| {
        let mut native =
            NativeEvaluator::<DoubleFloat>::with_binary_precision(ScalarIntegral::D0, 106).unwrap();
        let mut failures = Vec::new();
        for case in &CASES {
            for order in labellings() {
                let input: Vec<Complex<DoubleFloat>> = relabel(case, order)
                    .into_iter()
                    .map(|[re, im]| Complex::new(DoubleFloat::from(re), DoubleFloat::from(im)))
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
                check(values, case, 1e-30, &label, &mut failures);
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    });
}

#[test]
fn float_routes_agree_in_every_labelling() {
    spawn(|| {
        let bits = 128;
        let mut native =
            NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
        let mut interpreted = PrecisionEvaluator::with_binary_precision_and_backend(
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
                    .map(|[re, im]| {
                        Complex::new(Float::with_val(bits, re), Float::with_val(bits, im))
                    })
                    .collect();
                for route in ["native", "interpreted"] {
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(Float::new(bits), Float::new(bits))
                    });
                    if route == "native" {
                        native.evaluate(&input, &mut output).unwrap();
                    } else {
                        interpreted.evaluate(&input, &mut output).unwrap();
                    }
                    let label = format!("{}, {route} Float {bits}, labelling {order:?}", case.name);
                    check(
                        output.map(|z| [z.re, z.im]),
                        case,
                        1e-35,
                        &label,
                        &mut failures,
                    );
                }
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    });
}
