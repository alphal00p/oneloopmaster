//! Binary64 D0 next to the two-particle threshold of two massive propagators.
//!
//! The sector formulas take the masses' square roots: `k = (m_i + m_j - p) /
//! (√m_i √m_j)` (squared masses) and `√(k² - 4)`, or the contour's edge
//! discriminants. Next to a threshold `p = (√m_i + √m_j)²` the binary64 kernels
//! lost about `2^-53 |a| / √|λ|` (`a = m_i + m_j - p`, `λ = a² - 4 m_i m_j`),
//! the three-mass sectors about `2^-53 a² / |λ|`, and at a threshold of masses
//! without exact square roots `√(k² - 4)` came out near `√(2^-53)` instead of
//! zero. Such points are now evaluated by the double-double kernels and rounded
//! once (`d0_next_to_threshold` in `src/native/mod.rs`).
//!
//! References: native code in `Float` at 256 bits on the same binary64 inputs.
#![cfg(feature = "generated-evaluators")]
use oneloop::{EvaluationBackend, NativeEvaluator, ScalarIntegral, evaluate_with_backend};
use symbolica::domains::float::{Complex, Float, RealLike};

type C = Complex<f64>;

const BITS: u32 = 256;

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

/// The worst relative error of the finite part over every binary64 route
/// (native hand-written sectors and generated code, expressions, SymJIT), per
/// point.
fn binary64_errors(points: &[[f64; 11]]) -> Vec<(f64, [f64; 11])> {
    let mut precise =
        NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, BITS).unwrap();
    let mut native = [true, false].map(|handwritten| {
        let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
        evaluator.set_handwritten(handwritten);
        evaluator
    });
    points
        .iter()
        .map(|point| {
            let high: Vec<Complex<Float>> = point
                .iter()
                .map(|&x| Complex::new(Float::with_val(BITS, x), Float::new(BITS)))
                .collect();
            let mut out = vec![Complex::new(Float::new(BITS), Float::new(BITS)); 3];
            precise.evaluate(&high, &mut out).unwrap();
            let exact = C::new(out[0].re.to_f64(), out[0].im.to_f64());
            let input: Vec<C> = point.iter().map(|&x| C::new(x, 0.0)).collect();
            let mut values: Vec<C> = native
                .iter_mut()
                .map(|evaluator| {
                    let mut output = vec![C::new(f64::NAN, f64::NAN); 3];
                    evaluator.evaluate(&input, &mut output).unwrap();
                    output[0]
                })
                .collect();
            for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluate_with_backend(ScalarIntegral::D0, &input, &mut output, backend).unwrap();
                values.push(output[0]);
            }
            let worst = values
                .iter()
                .map(|value| {
                    (value.re - exact.re).hypot(value.im - exact.im) / exact.re.hypot(exact.im)
                })
                .fold(0.0, |worst: f64, error| {
                    if error.is_nan() { f64::INFINITY } else { worst.max(error) }
                });
            (worst, *point)
        })
        .collect()
}

fn check(points: &[[f64; 11]], tolerance: f64) {
    let failures: Vec<_> = binary64_errors(points)
        .into_iter()
        .filter(|(error, _)| !(*error <= tolerance))
        .collect();
    assert!(failures.is_empty(), "{failures:#?}");
}

/// `point` with the invariant `slot` at `(1 + r)` times its value, for r = 0 and
/// r = ±1e-15, ±1e-12, ±1e-9, ±1e-6.
fn neighbourhood(point: [f64; 11], slot: usize) -> Vec<[f64; 11]> {
    let mut points = vec![point];
    for r in [1e-15, 1e-12, 1e-9, 1e-6] {
        for sign in [1.0, -1.0] {
            let mut moved = point;
            moved[slot] = point[slot] * (1.0 + sign * r);
            points.push(moved);
        }
    }
    points
}

/// Squared mass of the heavy quarks of two onia at threshold (two quark pairs
/// at P/2 each, every momentum exact in binary64); it has no exact square root.
const M: f64 = 2.78515625;

/// Two boxes of a seven-point integral with two onia at threshold (TIRust's
/// two-onia heptagon, its second onium moved off threshold by 47%):
/// - masses `(M, M, 0, 0)` with `p1² = 4 M`, the generated code's sector;
/// - four masses `(1, M, M, 1.47 M)` with `p2² = 4 M`, the hand-written
///   Denner-Nierste-Scharf sector.
///
/// Binary64 was 5.3e-9 and 4.0e-9 off at the exact thresholds, 5e-10 and
/// 7e-10 one part in 1e15 from them, 2e-11 at 1e-12, 7e-13 at 1e-9 and 3e-14
/// at 1e-6 (double-double exact throughout). Now within 1e-15.
#[test]
fn onia_boxes_at_their_thresholds_in_binary64() {
    on_stack(|| {
        let two_mass =
            [11.140625, M, -43.02734375, -32.34375, 29.50390625, -57.265625, M, M, 0.0, 0.0, 4.0];
        let four_mass =
            [-3.203125, 4.0 * M, 4.0 * M, 15.859375, -12.875, 53.4375, 1.0, M, M, 1.47 * M, 4.0];
        let mut points = neighbourhood(two_mass, 0);
        points.extend(neighbourhood(four_mass, 1));
        check(&points, 1e-15);
    });
}

/// The four-mass box with both of its quark pairs at threshold (`p2² = p3² =
/// 4 M`, the second onium on threshold too): binary64 was 1.7e-8 off.
#[test]
fn box_with_two_thresholds_in_binary64() {
    on_stack(|| {
        let point = [-3.203125, 4.0 * M, 4.0 * M, 15.859375, -12.875, 53.4375, 1.0, M, M, M, 4.0];
        check(&[point], 1e-15);
    });
}

/// Three-mass boxes next to a threshold of a channel (masses with exact square
/// roots): `s23` 3.2e-9 above `(2 + 3)²` with masses `(4, 4, 0, 9)`, and `s12`
/// 4.2e-9 above `(0.75 + 0.5)²` with masses `(0.5625, 0, 0.25, 1)`. Binary64
/// was 1.5e-11 and 2.2e-9 off (about `2^-53 a² / |λ|`). Now within 1e-15.
#[test]
fn three_mass_boxes_next_to_a_channel_threshold_in_binary64() {
    on_stack(|| {
        let points = [
            [
                3.890625, 39.890625, 24.3125, 28.625, 11.984375, 25.000000079118351,
                4.0, 4.0, 0.0, 9.0, 4.0,
            ],
            [
                18.8125, -11.90625, -8.328125, 27.234375, 1.5625000065335182, 8.921875,
                0.5625, 0.0, 0.25, 1.0, 4.0,
            ],
        ];
        check(&points, 1e-15);
    });
}

/// Pseudo-thresholds are not singular on the physical sheet: equal masses
/// next to a vanishing leg (a massless leg between two heavy-quark
/// propagators) keep the binary64 kernels and stay within 1e-13.
#[test]
fn equal_masses_at_a_vanishing_leg_in_binary64() {
    on_stack(|| {
        let point = [0.0, 3.890625, 24.3125, 28.625, 11.984375, -7.25, M, M, 0.0, 9.0, 4.0];
        let mut points = neighbourhood(point, 1);
        points[0][0] = 0.0;
        check(&points, 1e-13);
    });
}
