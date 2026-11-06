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
//! hand-written transcription). The references are native `Float` at 512
//! bits (the generated code, which keeps more than 400 bits here), which must
//! agree in every order, and for Euclidean points the Usyukina-Davydychev
//! closed form (60 digits, mpmath).
#![cfg(feature = "generated-evaluators")]
use oneloop::{NativeEvaluator, NativeFloat, ScalarIntegral};
use symbolica::domains::float::{Complex, DoubleFloat, Float, Real, RealLike};

const BITS: u32 = 512;
/// Working precision of the native `Float` references.
const REFERENCE_BITS: u32 = 512;

struct Rng(u64);

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
fn square(vectors: &[[f64; 4]]) -> Float {
    let component = |i: usize| vectors.iter().fold(exact(0.0), |sum, v| sum + &exact(v[i]));
    let [e, x, y, z] = [0, 1, 2, 3].map(component);
    e.clone() * &e - x.clone() * &x - y.clone() * &y - z.clone() * &z
}

fn wide(value: &DoubleFloat) -> Float {
    let pair = value.into_inner();
    exact(pair.hi()) + &exact(pair.lo())
}

/// The finite part of a massless C0 at `mu² = 1`.
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

/// Worst double-double and binary64 errors over the six orders of `p`, each
/// against the 512-bit value at its own (rounded) arguments, or against
/// `reference` where given. The 512-bit values must agree across the orders.
fn worst_over_orders(p: &[Float; 3], reference: Option<&[Float; 2]>) -> (f64, f64) {
    let (mut worst_dd, mut worst_binary64) = (0.0f64, 0.0f64);
    let mut first: Option<[Float; 2]> = None;
    for order in orders(p) {
        let dd = order.clone().map(|x| x.to_double_float());
        let dd_exact = dd.clone().map(|x| wide(&x));
        let reference_value = finite(
            &dd_exact.map(|x| x.at_precision(REFERENCE_BITS)),
            REFERENCE_BITS,
        );
        let reference_value =
            [reference_value.re, reference_value.im].map(|x| x.at_precision(BITS));
        match &first {
            None => first = Some(reference_value.clone()),
            Some(first) => assert!(
                relative(&reference_value, first) < 1e-60,
                "{REFERENCE_BITS} bits disagree across orders at {:?}",
                p.clone().map(|x| x.to_f64())
            ),
        }
        let value = finite(&dd, 106);
        let error = relative(
            &[wide(&value.re), wide(&value.im)],
            reference.unwrap_or(&reference_value),
        );
        worst_dd = worst_dd.max(error);

        let binary64 = order.map(|x| x.to_f64());
        let reference_value = finite(
            &binary64.map(|x| Float::with_val(REFERENCE_BITS, x)),
            REFERENCE_BITS,
        );
        let reference_value =
            [reference_value.re, reference_value.im].map(|x| x.at_precision(BITS));
        let value = finite(&binary64, 53);
        let error = relative(
            &[exact(value.re), exact(value.im)],
            reference.unwrap_or(&reference_value),
        );
        worst_binary64 = worst_binary64.max(error);
    }
    (worst_dd, worst_binary64)
}

fn on_stack(action: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(action)
        .unwrap()
        .join()
        .unwrap();
}

/// Double-double and binary64 bounds, far below the defect (1e-2 to 1.0 in
/// double-double) and above the conditioning of the formula: with one
/// light-like leg the `log² p1` terms cancel between the dilogarithms and
/// the logarithms (three digits at `p1 ~ 1e-17 s`).
const DOUBLE_DOUBLE: f64 = 1e-27;
const BINARY64: f64 = 1e-12;

fn check(label: &str, points: Vec<[Float; 3]>) {
    let (mut worst_dd, mut worst_binary64) = (0.0f64, 0.0f64);
    let mut failures = vec![];
    for p in &points {
        let (dd, binary64) = worst_over_orders(p, None);
        worst_dd = worst_dd.max(dd);
        worst_binary64 = worst_binary64.max(binary64);
        if !(dd <= DOUBLE_DOUBLE && binary64 <= BINARY64) {
            failures.push(format!(
                "{:?}: double-double {dd:.2e}, binary64 {binary64:.2e}",
                p.clone().map(|x| x.to_f64())
            ));
        }
    }
    eprintln!(
        "{label}: {} points x 6 orders: worst double-double {worst_dd:.2e}, binary64 {worst_binary64:.2e}",
        points.len()
    );
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
fn two_light_like_legs_in_every_order() {
    on_stack(|| {
        let mut rng = Rng(0x9e3779b97f4a7c15);
        let points = (0..40)
            .map(|_| {
                let (k1, k2) = (massless(&mut rng), massless(&mut rng));
                [square(&[k1]), square(&[k2]), square(&[k1, k2])]
            })
            .collect();
        check("two light-like legs", points);
    });
}

/// One binary64 massless leg next to two massive ones.
#[test]
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
        check("one light-like leg", points);
    });
}

/// Invariants 1e-3 to 1e-14 times the largest, of either sign, in every order.
#[test]
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
        check("hierarchical invariants", points);
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
        for (p, value) in cases {
            let reference = [Float::parse(value, Some(BITS)).unwrap(), exact(0.0)];
            let (dd, binary64) = worst_over_orders(&p, Some(&reference));
            eprintln!("{value}: double-double {dd:.2e}, binary64 {binary64:.2e}");
            assert!(
                dd <= 1e-30 && binary64 <= 1e-14,
                "{value}: double-double {dd:e}, binary64 {binary64:e}"
            );
        }
    });
}
