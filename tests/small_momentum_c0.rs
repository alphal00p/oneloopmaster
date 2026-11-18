#![cfg(feature = "generated-evaluators")]
//! Small-momentum C0: binary64 and double-double finite parts against the
//! general expression at 1024 bits, for external invariants from zero and
//! 1e-14 to 1/80 of the largest squared mass (one or none massless, equal,
//! separated, nearby, hierarchical and complex masses; one, two or three
//! nonzero invariants of either sign, in every propagator order).
use oneloop::{NativeEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, DoubleFloat, Float, RealLike};

const BITS: u32 = 1024;

type Point = [(f64, f64); 7];

fn reference(input: &Point) -> [Complex<Float>; 3] {
    let mut native =
        NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::C0, BITS).unwrap();
    native.set_small_momentum(false);
    let point: Vec<_> = input
        .iter()
        .map(|&(re, im)| Complex::new(Float::with_val(BITS, re), Float::with_val(BITS, im)))
        .collect();
    let mut out = vec![Complex::new(Float::new(BITS), Float::new(BITS)); 3];
    native.evaluate(&point, &mut out).unwrap();
    [out[0].clone(), out[1].clone(), out[2].clone()]
}

/// Largest slot difference over the finite part's modulus.
fn relative<T: RealLike>(value: &[Complex<T>], reference: &[Complex<Float>; 3]) -> f64 {
    let scale = reference[0].re.to_f64().hypot(reference[0].im.to_f64());
    let mut worst = 0.0f64;
    for slot in 0..3 {
        let dr = Float::with_val(BITS, value[slot].re.to_f64()) - reference[slot].re.clone();
        let di = Float::with_val(BITS, value[slot].im.to_f64()) - reference[slot].im.clone();
        worst = worst.max(dr.to_f64().hypot(di.to_f64()));
    }
    worst / scale
}

/// Like `relative`, with both double-double limbs.
fn relative_dd(value: &[Complex<DoubleFloat>], reference: &[Complex<Float>; 3]) -> f64 {
    let wide = |x: &DoubleFloat| {
        let d = x.clone().into_inner();
        Float::with_val(BITS, d.hi()) + Float::with_val(BITS, d.lo())
    };
    let scale = reference[0].re.to_f64().hypot(reference[0].im.to_f64());
    let mut worst = 0.0f64;
    for slot in 0..3 {
        let dr = wide(&value[slot].re) - reference[slot].re.clone();
        let di = wide(&value[slot].im) - reference[slot].im.clone();
        worst = worst.max(dr.to_f64().hypot(di.to_f64()));
    }
    worst / scale
}

/// The series parameter ρ = Σ |s_ij| / (Re m_i² + Re m_j²) of `small_momentum_c0`.
fn rho(point: &Point) -> f64 {
    let re = |k: usize| point[k].0;
    [(0, 3, 4), (1, 4, 5), (2, 3, 5)]
        .iter()
        .map(|&(s, i, j)| re(s).abs() / (re(i) + re(j)))
        .sum()
}

fn cases() -> Vec<Point> {
    let mu = (1.0e6, 0.0);
    let (m, w) = (6467.2, 0.0);
    let top = (29929.0, -233.6);
    let masses: Vec<[(f64, f64); 3]> = vec![
        [(0.0, 0.0), (m, w), (m, w)],
        [(m, w), (m, w), (m, w)],
        [(m, w), (2.5 * m, w), (7.0 * m, w)],
        [(m, w), (m, w), (4.6 * m, w)],
        [(0.0, 0.0), (m, w), (3.0 * m, w)],
        [(0.0, 0.0), (m, w), (1.286 * m, w)],
        [(m, w), (1.001 * m, w), (3.0 * m, w)],
        [(m, w), (1.1 * m, w), (1.25 * m, w)],
        [(1.0e-6 * m, w), (m, w), (2.0 * m, w)],
        [(0.0, 0.0), top, top],
        [top, top, top],
        [(m, -130.0), (m, -130.0), top],
    ];
    let mut out = Vec::new();
    for masses in masses {
        let big = masses.iter().map(|z| z.0.hypot(z.1)).fold(0.0, f64::max);
        let mut invariants: Vec<[f64; 3]> = vec![[0.0; 3]];
        for exponent in [-14, -11, -8, -6, -4, -3, -2] {
            let t = big * 10f64.powi(exponent);
            for shape in [
                [0.0, 0.0, -1.2345],
                [0.0, 0.0, 0.987],
                [0.42, 0.0, -1.1],
                [0.31, -0.77, 0.52],
            ] {
                if exponent == -2 && shape[1] != 0.0 {
                    continue; // ρ above the series' 1/16 for three invariants.
                }
                invariants.push([shape[0] * t, shape[1] * t, shape[2] * t]);
            }
        }
        for s in invariants {
            // Every relabelling of the propagators (legs follow the masses).
            for order in [
                [0, 1, 2],
                [1, 2, 0],
                [2, 0, 1],
                [0, 2, 1],
                [2, 1, 0],
                [1, 0, 2],
            ] {
                // Original invariants: s[0] between propagators (0, 1), s[1] (1, 2), s[2] (0, 2).
                let original = |x: usize, y: usize| match (x.min(y), x.max(y)) {
                    (0, 1) => s[0],
                    (1, 2) => s[1],
                    _ => s[2],
                };
                let inv = |a: usize, b: usize| original(order[a], order[b]);
                let m3 = order.map(|l| masses[l]);
                out.push([
                    (inv(0, 1), 0.0),
                    (inv(1, 2), 0.0),
                    (inv(0, 2), 0.0),
                    m3[0],
                    m3[1],
                    m3[2],
                    mu,
                ]);
            }
        }
    }
    out
}

#[test]
fn small_momentum_c0_is_accurate_in_binary64_and_double_double() {
    let mut double = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
    let mut dd = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::C0).unwrap();
    let mut general = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
    general.set_small_momentum(false);
    let (mut worst64, mut worst_dd, mut worst_general, mut worst_dd_general) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let cases = cases();
    for case in &cases {
        let exact = reference(case);
        let point: Vec<_> = case.iter().map(|&(re, im)| Complex::new(re, im)).collect();
        let mut out = vec![Complex::new(0.0, 0.0); 3];
        double.evaluate(&point, &mut out).unwrap();
        let e64 = relative(&out, &exact);
        general.evaluate(&point, &mut out).unwrap();
        let eg = relative(&out, &exact);
        let point: Vec<_> = case
            .iter()
            .map(|&(re, im)| Complex::new(DoubleFloat::from(re), DoubleFloat::from(im)))
            .collect();
        let mut out = vec![Complex::new(DoubleFloat::from(0.0), DoubleFloat::from(0.0)); 3];
        dd.evaluate(&point, &mut out).unwrap();
        let edd = relative_dd(&out, &exact);
        if std::env::var("VERBOSE").is_ok() {
            println!("{case:?}: binary64 {e64:.2e} general {eg:.2e} double-double {edd:.2e}");
        }
        worst64 = worst64.max(e64);
        // Double-double takes the series up to ρ ≈ 1/415 (12 orders), binary64 up
        // to ρ ≈ 1/24; beyond, the general expression, which loses a few bits.
        if rho(case) <= 1.0 / 500.0 {
            worst_dd = worst_dd.max(edd);
        } else {
            worst_dd_general = worst_dd_general.max(edd);
        }
        worst_general = worst_general.max(eg);
    }
    println!(
        "{} points: worst binary64 {worst64:e}, double-double {worst_dd:e} (beyond its series {worst_dd_general:e}); general expression {worst_general:e}",
        cases.len()
    );
    assert!(worst64 < 4e-15, "binary64 {worst64:e}");
    assert!(worst_dd < 1e-29, "double-double {worst_dd:e}");
    // The general expression is what the small-momentum path repairs.
    assert!(worst_general > 1e-6);
}

/// Outside the series' domain the evaluator is the general expression, bit
/// for bit: ρ above 1/16, two massless propagators (collinear or soft
/// singularities), and physical invariants.
#[test]
fn outside_the_series_domain_c0_is_the_general_expression() {
    let mut double = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
    let mut general = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
    general.set_small_momentum(false);
    let (m, mu) = (6467.2, 1.0e6);
    let points: [[f64; 7]; 5] = [
        [-900.0, 0.0, 300.0, m, m, m, mu],
        [0.0, 0.0, -1.0, 0.0, 0.0, m, mu],
        [-1.0e-3, 0.0, 0.0, m, 0.0, 0.0, mu],
        [29929.0, 0.0, 1.0e5, 29929.0, 29929.0, 0.0, mu],
        [-2.0e4, 1.0e3, -5.0e3, m, 2.0 * m, 3.0 * m, mu],
    ];
    for point in points {
        let input: Vec<_> = point.iter().map(|&x| Complex::new(x, 0.0)).collect();
        let (mut a, mut b) = (
            vec![Complex::new(0.0, 0.0); 3],
            vec![Complex::new(0.0, 0.0); 3],
        );
        double.evaluate(&input, &mut a).unwrap();
        general.evaluate(&input, &mut b).unwrap();
        for slot in 0..3 {
            assert_eq!(
                a[slot].re.to_bits(),
                b[slot].re.to_bits(),
                "{point:?} slot {slot}"
            );
            assert_eq!(
                a[slot].im.to_bits(),
                b[slot].im.to_bits(),
                "{point:?} slot {slot}"
            );
        }
    }
}
