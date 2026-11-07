//! The massless triangle (OneLOop's trif0) next to light-like legs, in every
//! order of its arguments.
//!
//! A massless triangle `q0 = 0, q1 = k1, q2 = k1 + k2` with binary64 massless
//! momenta `k1, k2` has `p1 = k1²` and `p2 = k2²` of about `1e-16 s`, not zero.
//! The integral is the same for all six orders of its propagators, but the
//! direct forms lost the small root of the quadratic in `b ± d` and the
//! complement `1 - x` of a dilogarithm argument next to one: in
//! double-double the value was off by up to the whole value (relative error
//! 1.0) in some orders and by 1e-2 in all of them; binary64 likewise
//! (`triangle::triangle_finite_massless`; binary64 and double-double run its
//! hand-written transcription). With one invariant small against two nearly
//! equal ones (gap `g`), the orders that put the small one first lost up to
//! `30 eps/g²` relative: in `λ = b² - 4a`, in the sums `b`, `u`, `v`, and in
//! the winding term `i pi n/(1 - X1/X2)` of the quotient of the roots
//! (`triangle::massless_kallen`, `triangle::ordered_sum`).
//!
//! Every route is checked: native binary64 and double-double (the
//! hand-written kernel), and the interpreted expressions
//! (`EvaluationBackend::Expression`, the default without
//! `generated-evaluators`) in binary64 and at 256 bits. The references are
//! native `Float` at 512 bits (the generated code, which keeps more than 400
//! bits here), which must agree in every order, and for Euclidean points the
//! Usyukina-Davydychev closed form (60 digits, mpmath), which also runs
//! without `generated-evaluators` (expression routes only). The prebuilt
//! SymJIT C0 asset is not regenerated and keeps the direct forms.
#[cfg(feature = "generated-evaluators")]
use oneloop::NativeEvaluator;
use oneloop::{
    EvaluationBackend, NativeFloat, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
#[cfg(feature = "generated-evaluators")]
use symbolica::domains::float::DoubleFloat;
use symbolica::domains::float::{Complex, Float, Real, RealLike};

/// Precision of the exact invariants and of the native `Float` references.
const BITS: u32 = 512;
/// Working precision of the interpreted expressions' arbitrary-precision route.
const EXPRESSION_BITS: u32 = 256;

#[cfg(feature = "generated-evaluators")]
struct Rng(u64);

#[cfg(feature = "generated-evaluators")]
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A binary64 massless momentum: `E` is the rounded length of the spatial part.
#[cfg(feature = "generated-evaluators")]
fn massless(rng: &mut Rng) -> [f64; 4] {
    let p = [0, 1, 2].map(|_| 2.0 * rng.uniform() - 1.0);
    [
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt(),
        p[0],
        p[1],
        p[2],
    ]
}

fn exact(value: f64) -> Float {
    Float::with_val(BITS, value)
}

/// The exact Minkowski square of a sum of binary64 momenta.
#[cfg(feature = "generated-evaluators")]
fn square(vectors: &[[f64; 4]]) -> Float {
    let component = |i: usize| vectors.iter().fold(exact(0.0), |sum, v| sum + &exact(v[i]));
    let [e, x, y, z] = [0, 1, 2, 3].map(component);
    e.clone() * &e - x.clone() * &x - y.clone() * &y - z.clone() * &z
}

#[cfg(feature = "generated-evaluators")]
fn wide(value: &DoubleFloat) -> Float {
    let pair = value.into_inner();
    exact(pair.hi()) + &exact(pair.lo())
}

/// The finite part of a massless C0 at `mu² = 1`, native code.
#[cfg(feature = "generated-evaluators")]
fn finite<T: NativeFloat>(p: &[T; 3], bits: u32) -> Complex<T> {
    let mut evaluator =
        NativeEvaluator::<T>::with_binary_precision(ScalarIntegral::C0, bits).unwrap();
    let zero = p[0].zero();
    let mut input: Vec<_> = p
        .iter()
        .map(|x| Complex::new(x.clone(), zero.clone()))
        .collect();
    input.extend((0..3).map(|_| Complex::new(zero.clone(), zero.clone())));
    input.push(Complex::new(zero.one(), zero.clone()));
    let mut output = vec![Complex::new(zero.clone(), zero.clone()); 3];
    evaluator.evaluate(&input, &mut output).unwrap();
    output.swap_remove(0)
}

/// The finite part through the interpreted expressions in binary64.
fn expression_binary64(p: &[f64; 3]) -> [Float; 2] {
    let mut input: Vec<_> = p.iter().map(|&x| Complex::new(x, 0.0)).collect();
    input.extend([Complex::new(0.0, 0.0); 3]);
    input.push(Complex::new(1.0, 0.0));
    let mut output = [Complex::new(0.0, 0.0); 3];
    evaluate_with_backend(
        ScalarIntegral::C0,
        &input,
        &mut output,
        EvaluationBackend::Expression,
    )
    .unwrap();
    [exact(output[0].re), exact(output[0].im)]
}

/// The interpreted expressions at [EXPRESSION_BITS].
fn expression_evaluator() -> PrecisionEvaluator {
    PrecisionEvaluator::with_binary_precision_and_backend(
        ScalarIntegral::C0,
        EXPRESSION_BITS,
        EvaluationBackend::Expression,
    )
    .unwrap()
}

/// The finite part through `evaluator` at binary64 arguments.
fn expression_precise(evaluator: &mut PrecisionEvaluator, p: &[f64; 3]) -> [Float; 2] {
    let value = |x: f64| {
        Complex::new(
            Float::with_val(EXPRESSION_BITS, x),
            Float::new(EXPRESSION_BITS),
        )
    };
    let mut input: Vec<_> = p.iter().map(|&x| value(x)).collect();
    input.extend([0.0, 0.0, 0.0, 1.0].map(value));
    let mut output: Vec<_> = (0..3).map(|_| value(0.0)).collect();
    evaluator.evaluate(&input, &mut output).unwrap();
    let finite = output.swap_remove(0);
    [finite.re, finite.im].map(|x| x.at_precision(BITS))
}

fn relative(value: &[Float; 2], reference: &[Float; 2]) -> f64 {
    let (re, im) = (
        value[0].clone() - &reference[0],
        value[1].clone() - &reference[1],
    );
    let error = (re.clone() * &re + im.clone() * &im).sqrt();
    let scale =
        (reference[0].clone() * &reference[0] + reference[1].clone() * &reference[1]).sqrt();
    let error = (error / scale).to_f64();
    if error.is_nan() { f64::INFINITY } else { error }
}

/// The six orders of `[p1, p2, p3] = [s01, s12, s20]`: the propagator orders
/// `(a, b, c)` give `[s_ab, s_bc, s_ca]`.
fn orders<V: Clone>(p: &[V; 3]) -> Vec<[V; 3]> {
    let pair = |i: usize, j: usize| match (i.min(j), i.max(j)) {
        (0, 1) => p[0].clone(),
        (1, 2) => p[1].clone(),
        _ => p[2].clone(),
    };
    [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
        [0, 2, 1],
        [2, 1, 0],
        [1, 0, 2],
    ]
    .map(|[a, b, c]| [pair(a, b), pair(b, c), pair(c, a)])
    .to_vec()
}

/// Worst relative errors over the six orders of a point, per route.
#[derive(Default)]
struct Worst {
    double_double: f64,
    binary64: f64,
    expression: f64,
    expression_precise: f64,
}

impl Worst {
    #[cfg(feature = "generated-evaluators")]
    fn max(&mut self, other: &Self) {
        self.double_double = self.double_double.max(other.double_double);
        self.binary64 = self.binary64.max(other.binary64);
        self.expression = self.expression.max(other.expression);
        self.expression_precise = self.expression_precise.max(other.expression_precise);
    }

    fn within(&self, bounds: &Self) -> bool {
        self.double_double <= bounds.double_double
            && self.binary64 <= bounds.binary64
            && self.expression <= bounds.expression
            && self.expression_precise <= bounds.expression_precise
    }
}

impl std::fmt::Display for Worst {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if cfg!(feature = "generated-evaluators") {
            write!(
                f,
                "double-double {:.2e}, binary64 {:.2e}, ",
                self.double_double, self.binary64
            )?;
        }
        write!(
            f,
            "expression binary64 {:.2e}, expression {EXPRESSION_BITS} bits {:.2e}",
            self.expression, self.expression_precise
        )
    }
}

/// The 512-bit native value at the arguments `p`.
#[cfg(feature = "generated-evaluators")]
fn native_reference(p: &[Float; 3]) -> [Float; 2] {
    let value = finite(&p.clone().map(|x| x.at_precision(BITS)), BITS);
    [value.re, value.im]
}

/// Worst errors over the six orders of `p`, each against the 512-bit native
/// value at its own (rounded) arguments, or against `reference` where given:
/// double-double at the double-double arguments, the binary64 routes (native,
/// expression in binary64 and at [EXPRESSION_BITS]) at the binary64 ones.
/// The 512-bit values must agree across the orders.
#[cfg(feature = "generated-evaluators")]
fn worst_over_orders(
    p: &[Float; 3],
    reference: Option<&[Float; 2]>,
    expression: &mut PrecisionEvaluator,
) -> Worst {
    let mut worst = Worst::default();
    let mut first: Option<[Float; 2]> = None;
    for order in orders(p) {
        let dd = order.clone().map(|x| x.to_double_float());
        let reference_value = native_reference(&dd.map(|x| wide(&x)));
        match &first {
            None => first = Some(reference_value.clone()),
            Some(first) => assert!(
                relative(&reference_value, first) < 1e-60,
                "{BITS} bits disagree across orders at {:?}",
                p.clone().map(|x| x.to_f64())
            ),
        }
        let value = finite(&dd, 106);
        let error = relative(
            &[wide(&value.re), wide(&value.im)],
            reference.unwrap_or(&reference_value),
        );
        worst.double_double = worst.double_double.max(error);

        let binary64 = order.map(|x| x.to_f64());
        let reference_value = native_reference(&binary64.map(exact));
        let reference_value = reference.unwrap_or(&reference_value);
        let value = finite(&binary64, 53);
        let error = relative(&[exact(value.re), exact(value.im)], reference_value);
        worst.binary64 = worst.binary64.max(error);
        let error = relative(&expression_binary64(&binary64), reference_value);
        worst.expression = worst.expression.max(error);
        let error = relative(&expression_precise(expression, &binary64), reference_value);
        worst.expression_precise = worst.expression_precise.max(error);
    }
    worst
}

/// The expression routes alone over the six orders of binary64 `p`, against
/// `reference` (without `generated-evaluators`).
#[cfg(not(feature = "generated-evaluators"))]
fn expression_over_orders(
    p: &[f64; 3],
    reference: &[Float; 2],
    expression: &mut PrecisionEvaluator,
) -> Worst {
    let mut worst = Worst::default();
    for order in orders(p) {
        let error = relative(&expression_binary64(&order), reference);
        worst.expression = worst.expression.max(error);
        let error = relative(&expression_precise(expression, &order), reference);
        worst.expression_precise = worst.expression_precise.max(error);
    }
    worst
}

fn on_stack(action: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(action)
        .unwrap()
        .join()
        .unwrap();
}

/// Bounds far below the defect (1e-2 to 1.0 in double-double and binary64)
/// and above the conditioning of the formula: with one light-like leg the
/// `log² p1` terms cancel between the dilogarithms and the logarithms (three
/// digits at `p1 ~ 1e-17 s`). Before the fix the 256-bit expressions kept
/// only 42 to 58 digits at these points (now 69 to 77).
#[cfg(feature = "generated-evaluators")]
const BOUNDS: Worst = Worst {
    double_double: 1e-27,
    binary64: 1e-12,
    expression: 1e-12,
    expression_precise: 1e-60,
};

/// [BOUNDS] scaled by `max(1, 1e-2/g)` where the two invariants of largest
/// modulus are a relative gap `g` apart: the trif0 representation's terms
/// grow like `1/g` there and cancel (a few `1e-15/g` in binary64 and
/// `3e-31/g` in double-double, in every order; the generated code of the
/// orders without the small invariant first has the same).
#[cfg(feature = "generated-evaluators")]
fn gap_scaled(p: &[Float; 3]) -> Worst {
    let mut p = p.clone().map(|x| x.to_f64());
    p.sort_by(|a, b| b.abs().total_cmp(&a.abs()));
    let scale = (1e-2 * p[0].abs() / (p[0] - p[1]).abs()).max(1.0);
    Worst {
        double_double: BOUNDS.double_double * scale,
        binary64: BOUNDS.binary64 * scale,
        expression: BOUNDS.expression * scale,
        expression_precise: BOUNDS.expression_precise * scale,
    }
}

#[cfg(feature = "generated-evaluators")]
fn check(label: &str, points: Vec<[Float; 3]>, bounds: impl Fn(&[Float; 3]) -> Worst) {
    let mut expression = expression_evaluator();
    let mut worst = Worst::default();
    let mut failures = vec![];
    for p in &points {
        let errors = worst_over_orders(p, None, &mut expression);
        worst.max(&errors);
        if !errors.within(&bounds(p)) {
            failures.push(format!("{:?}: {errors}", p.clone().map(|x| x.to_f64())));
        }
    }
    eprintln!("{label}: {} points x 6 orders: worst {worst}", points.len());
    assert!(
        failures.is_empty(),
        "{label}: {} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// `q0 = 0, q1 = k1, q2 = k1 + k2` with binary64 massless `k1, k2`: two
/// light-like legs (the TIRust review's point).
#[test]
#[cfg(feature = "generated-evaluators")]
fn two_light_like_legs_in_every_order() {
    on_stack(|| {
        let mut rng = Rng(0x9e3779b97f4a7c15);
        let points = (0..40)
            .map(|_| {
                let (k1, k2) = (massless(&mut rng), massless(&mut rng));
                [square(&[k1]), square(&[k2]), square(&[k1, k2])]
            })
            .collect();
        check("two light-like legs", points, |_| BOUNDS);
    });
}

/// One binary64 massless leg next to two massive ones.
#[test]
#[cfg(feature = "generated-evaluators")]
fn one_light_like_leg_in_every_order() {
    on_stack(|| {
        let mut rng = Rng(0x2545f4914f6cdd1d);
        let points = (0..40)
            .map(|_| {
                let k1 = massless(&mut rng);
                let k2 = [
                    2.0 + rng.uniform(),
                    rng.uniform(),
                    rng.uniform(),
                    rng.uniform(),
                ];
                [square(&[k1]), square(&[k2]), square(&[k1, k2])]
            })
            .collect();
        check("one light-like leg", points, |_| BOUNDS);
    });
}

/// Invariants 1e-3 to 1e-14 times the largest, of either sign, in every order.
#[test]
#[cfg(feature = "generated-evaluators")]
fn hierarchical_invariants_in_every_order() {
    on_stack(|| {
        let mut rng = Rng(0x5851f42d4c957f2d);
        let sign = |rng: &mut Rng| if rng.next() & 1 == 0 { -1.0 } else { 1.0 };
        let points = (0..40)
            .map(|_| {
                let s = sign(&mut rng) * (0.5 + rng.uniform());
                let small =
                    |rng: &mut Rng| sign(rng) * 10f64.powf(-3.0 - 11.0 * rng.uniform()) * s.abs();
                let p1 = small(&mut rng);
                let p2 = if rng.next() & 1 == 0 {
                    small(&mut rng)
                } else {
                    6.0 * rng.uniform() - 3.0
                };
                [exact(p1), exact(p2), exact(s)]
            })
            .collect();
        check("hierarchical invariants", points, |_| BOUNDS);
    });
}

/// One invariant small against two nearly equal ones, in every order:
/// relative gap `g` from 1e-1 to 1e-4 and the small invariant `1e-2 g²` to
/// `1e-10 g²` times the others (either sign; `λ` stays away from zero). In
/// the orders with the small invariant first, `b² - 4a`, the sums `b`, `u`,
/// `v` and the winding term of `log(X1/X2)/(1 - X1/X2)` lost up to
/// `30 eps/g²` relative (the generated code: up to 6 relative at g = 1e-4).
/// Bounds: [gap_scaled].
#[test]
#[cfg(feature = "generated-evaluators")]
fn small_invariant_next_to_nearly_equal_ones_in_every_order() {
    on_stack(|| {
        let mut rng = Rng(0x94d049bb133111eb);
        let sign = |rng: &mut Rng| if rng.next() & 1 == 0 { -1.0 } else { 1.0 };
        let points = (0..60)
            .map(|_| {
                let s = sign(&mut rng) * (0.5 + rng.uniform());
                let gap = sign(&mut rng) * 10f64.powf(-1.0 - 3.0 * rng.uniform());
                let small =
                    sign(&mut rng) * gap * gap * 10f64.powf(-2.0 - 8.0 * rng.uniform()) * s.abs();
                [exact(small), exact(s), exact(s * (1.0 + gap))]
            })
            .collect();
        check(
            "small invariant next to nearly equal ones",
            points,
            gap_scaled,
        );
    });
}

/// Euclidean points against the Usyukina-Davydychev form, 40 digits of a
/// 60-digit mpmath evaluation (`C0 = Phi(p1/p3, p2/p3)/p3`).
#[test]
fn euclidean_points_match_the_closed_form_in_every_order() {
    on_stack(|| {
        let rational = |n: i64, log2_d: i32| exact(n as f64) * &exact(2f64.powi(-log2_d));
        let cases = [
            (
                [rational(-1, 55), rational(-3, 56), rational(-2, 0)],
                "-747.1259439885710731333279418644793640615",
            ),
            (
                [rational(-5, 60), rational(-7, 53), rational(-13, 2)],
                "-456.5311206079765291305840964389150849504",
            ),
            (
                [rational(-1, 0), rational(-2, 0), rational(-3, 0)],
                "-1.254138667411437103276142783391067434113",
            ),
        ];
        let mut expression = expression_evaluator();
        for (p, value) in cases {
            let reference = [Float::parse(value, Some(BITS)).unwrap(), exact(0.0)];
            #[cfg(feature = "generated-evaluators")]
            let worst = worst_over_orders(&p, Some(&reference), &mut expression);
            #[cfg(not(feature = "generated-evaluators"))]
            let worst = expression_over_orders(&p.map(|x| x.to_f64()), &reference, &mut expression);
            eprintln!("{value}: {worst}");
            assert!(
                worst.within(&Worst {
                    double_double: 1e-30,
                    binary64: 1e-14,
                    expression: 1e-14,
                    expression_precise: 1e-38,
                }),
                "{value}: {worst}"
            );
        }
    });
}
