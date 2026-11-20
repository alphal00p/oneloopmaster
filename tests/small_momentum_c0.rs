#![cfg(feature = "generated-evaluators")]
//! Small-momentum C0: binary64 and double-double finite parts against the
//! general expression at 1024 bits, for external invariants from zero and
//! 1e-14 to 1/16 of the squared masses (one or none massless, equal,
//! separated, nearby, hierarchical and complex masses; one, two or three
//! nonzero invariants of either sign, in every propagator order), the cases
//! of the three reviews of #21 with points on both sides of the far-mass
//! separation, of binary64's closeness and of the rescaling for squared
//! masses far below the largest, and random three-mass points (down to 1e-34
//! of the largest).
//!
//! Where the series serves a point (`c0_small_momentum` returns a value), it
//! must be accurate to a few units of the precision. Where it declines, the
//! evaluator must return the general expression, and the general expression
//! must still be accurate there.
use oneloop::{NativeEvaluator, ScalarIntegral, c0_small_momentum};
use symbolica::domains::float::{Complex, DoubleFloat, Float, RealLike};

const BITS: u32 = 1024;
/// Served points.
const SERVED_BINARY64: f64 = 4e-15;
const SERVED_DOUBLE_DOUBLE: f64 = 1e-30;
/// Declined points inside the series' range of orders (binary64 ρ ≤ 1/24,
/// double-double ρ ≤ 1/415): the general expression there.
const DECLINED_BINARY64: f64 = 1e-10;
const DECLINED_DOUBLE_DOUBLE: f64 = 1e-22;
/// The series' range of orders: the fixed mass sets are served throughout.
const ORDERS_BINARY64: f64 = 1.0 / 24.0;
const ORDERS_DOUBLE_DOUBLE: f64 = 1.0 / 415.0;

type Point = [(f64, f64); 7];

fn binary64_tolerance() -> f64 {
    2f64.powi(-57)
}

fn double_double_tolerance() -> DoubleFloat {
    DoubleFloat::from(2f64.powi(-55)) * DoubleFloat::from(2f64.powi(-55))
}

/// The series parameter ρ = Σ |s_ij| / (Re m_i² + Re m_j²) of `c0_small_momentum`.
fn rho(point: &Point) -> f64 {
    [(0, 3, 4), (1, 4, 5), (2, 3, 5)]
        .iter()
        .map(|&(s, i, j)| {
            if point[s].0 == 0.0 {
                0.0
            } else {
                point[s].0.abs() / (point[i].0 + point[j].0)
            }
        })
        .sum()
}

/// Relative errors against 1024 bits and the routes of one point.
#[derive(Debug)]
struct Outcome {
    binary64: f64,
    binary64_general: f64,
    double_double: f64,
    double_double_general: f64,
    served_binary64: bool,
    served_double_double: bool,
}

struct Evaluators {
    binary64: NativeEvaluator<f64>,
    binary64_general: NativeEvaluator<f64>,
    double_double: NativeEvaluator<DoubleFloat>,
    double_double_general: NativeEvaluator<DoubleFloat>,
    reference: NativeEvaluator<Float>,
}

fn wide(x: &DoubleFloat) -> Float {
    let d = (*x).into_inner();
    Float::with_val(BITS, d.hi()) + Float::with_val(BITS, d.lo())
}

/// Largest slot difference over the finite part's modulus.
fn relative(value: &[Complex<Float>], reference: &[Complex<Float>]) -> f64 {
    let scale = reference[0].re.to_f64().hypot(reference[0].im.to_f64());
    let mut worst = 0.0f64;
    for slot in 0..3 {
        let dr = value[slot].re.clone() - reference[slot].re.clone();
        let di = value[slot].im.clone() - reference[slot].im.clone();
        let error = dr.to_f64().hypot(di.to_f64());
        worst = if error.is_nan() {
            f64::INFINITY
        } else {
            worst.max(error)
        };
    }
    worst / scale
}

impl Evaluators {
    fn new() -> Self {
        let mut binary64_general = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
        binary64_general.set_small_momentum(false);
        let mut double_double_general =
            NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::C0).unwrap();
        double_double_general.set_small_momentum(false);
        let mut reference =
            NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::C0, BITS).unwrap();
        reference.set_small_momentum(false);
        Self {
            binary64: NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap(),
            binary64_general,
            double_double: NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::C0).unwrap(),
            double_double_general,
            reference,
        }
    }

    fn binary64(evaluator: &mut NativeEvaluator<f64>, point: &[Complex<f64>]) -> Vec<Complex<f64>> {
        let mut out = vec![Complex::new(0.0, 0.0); 3];
        evaluator.evaluate(point, &mut out).unwrap();
        out
    }

    fn double_double(
        evaluator: &mut NativeEvaluator<DoubleFloat>,
        point: &[Complex<DoubleFloat>],
    ) -> Vec<Complex<DoubleFloat>> {
        let mut out = vec![Complex::new(DoubleFloat::from(0.0), DoubleFloat::from(0.0)); 3];
        evaluator.evaluate(point, &mut out).unwrap();
        out
    }

    /// Evaluates `point` every way and checks that the evaluator returns
    /// what `c0_small_momentum` serves (finite part, zero poles) or else the
    /// general expression, bit for bit.
    fn check(&mut self, point: &Point) -> Outcome {
        let exact = {
            let input: Vec<_> = point
                .iter()
                .map(|&(re, im)| Complex::new(Float::with_val(BITS, re), Float::with_val(BITS, im)))
                .collect();
            let mut out = vec![Complex::new(Float::new(BITS), Float::new(BITS)); 3];
            self.reference.evaluate(&input, &mut out).unwrap();
            out
        };
        let input: Vec<_> = point.iter().map(|&(re, im)| Complex::new(re, im)).collect();
        let value = Self::binary64(&mut self.binary64, &input);
        let general = Self::binary64(&mut self.binary64_general, &input);
        let served = c0_small_momentum(&input, &binary64_tolerance());
        let expected = match &served {
            Some(finite) => vec![*finite, Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
            None => general.clone(),
        };
        for slot in 0..3 {
            assert!(
                value[slot].re.to_bits() == expected[slot].re.to_bits()
                    && value[slot].im.to_bits() == expected[slot].im.to_bits(),
                "{point:?}: binary64 evaluator slot {slot} is neither the series nor the general expression"
            );
        }
        let lift = |v: &[Complex<f64>]| -> Vec<Complex<Float>> {
            v.iter()
                .map(|z| Complex::new(Float::with_val(BITS, z.re), Float::with_val(BITS, z.im)))
                .collect()
        };
        let binary64 = relative(&lift(&value), &exact);
        let binary64_general = relative(&lift(&general), &exact);

        let input: Vec<_> = point
            .iter()
            .map(|&(re, im)| Complex::new(DoubleFloat::from(re), DoubleFloat::from(im)))
            .collect();
        let value = Self::double_double(&mut self.double_double, &input);
        let general = Self::double_double(&mut self.double_double_general, &input);
        let served_dd = c0_small_momentum(&input, &double_double_tolerance());
        let zero = Complex::new(DoubleFloat::from(0.0), DoubleFloat::from(0.0));
        let expected = match &served_dd {
            Some(finite) => vec![*finite, zero, zero],
            None => general.clone(),
        };
        for slot in 0..3 {
            assert!(
                value[slot] == expected[slot],
                "{point:?}: double-double evaluator slot {slot} is neither the series nor the general expression"
            );
        }
        let lift = |v: &[Complex<DoubleFloat>]| -> Vec<Complex<Float>> {
            v.iter()
                .map(|z| Complex::new(wide(&z.re), wide(&z.im)))
                .collect()
        };
        Outcome {
            binary64,
            binary64_general,
            double_double: relative(&lift(&value), &exact),
            double_double_general: relative(&lift(&general), &exact),
            served_binary64: served.is_some(),
            served_double_double: served_dd.is_some(),
        }
    }
}

/// Worst errors of a set of points, split by route.
#[derive(Default, Debug)]
struct Summary {
    points: usize,
    served_binary64: usize,
    served_double_double: usize,
    worst_served_binary64: f64,
    worst_served_double_double: f64,
    worst_declined_binary64: f64,
    worst_declined_double_double: f64,
    worst_general_binary64: f64,
    smallest_declined_rho_binary64: f64,
    smallest_declined_rho_double_double: f64,
}

impl Summary {
    fn add(&mut self, point: &Point, outcome: &Outcome) {
        if std::env::var("VERBOSE").is_ok() {
            println!("{point:?}: {outcome:?}");
        }
        if self.points == 0 {
            self.smallest_declined_rho_binary64 = f64::INFINITY;
            self.smallest_declined_rho_double_double = f64::INFINITY;
        }
        self.points += 1;
        let rho = rho(point);
        if !outcome.served_binary64 {
            self.smallest_declined_rho_binary64 = self.smallest_declined_rho_binary64.min(rho);
        }
        if !outcome.served_double_double {
            self.smallest_declined_rho_double_double =
                self.smallest_declined_rho_double_double.min(rho);
        }
        if outcome.served_binary64 {
            self.served_binary64 += 1;
            self.worst_served_binary64 = self.worst_served_binary64.max(outcome.binary64);
        } else if rho <= 1.0 / 24.0 {
            self.worst_declined_binary64 = self.worst_declined_binary64.max(outcome.binary64);
        }
        if outcome.served_double_double {
            self.served_double_double += 1;
            self.worst_served_double_double =
                self.worst_served_double_double.max(outcome.double_double);
        } else if rho <= 1.0 / 415.0 {
            self.worst_declined_double_double =
                self.worst_declined_double_double.max(outcome.double_double);
        }
        self.worst_general_binary64 = self.worst_general_binary64.max(outcome.binary64_general);
        let _ = outcome.double_double_general;
    }

    fn assert_accurate(&self, name: &str) {
        println!("{name}: {self:?}");
        assert!(
            self.worst_served_binary64 < SERVED_BINARY64,
            "{name}: served binary64 {:e}",
            self.worst_served_binary64
        );
        assert!(
            self.worst_served_double_double < SERVED_DOUBLE_DOUBLE,
            "{name}: served double-double {:e}",
            self.worst_served_double_double
        );
        assert!(
            self.worst_declined_binary64 < DECLINED_BINARY64,
            "{name}: declined binary64 {:e}",
            self.worst_declined_binary64
        );
        assert!(
            self.worst_declined_double_double < DECLINED_DOUBLE_DOUBLE,
            "{name}: declined double-double {:e}",
            self.worst_declined_double_double
        );
    }
}

/// Every relabelling of the propagators of `masses` and invariants `s`
/// (`s[0]` between propagators 0 and 1, `s[1]` between 1 and 2, `s[2]`
/// between 0 and 2; legs follow the masses).
fn relabellings(masses: [(f64, f64); 3], s: [f64; 3]) -> Vec<Point> {
    let mu = (1.0e6, 0.0);
    [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
        [0, 2, 1],
        [2, 1, 0],
        [1, 0, 2],
    ]
    .iter()
    .map(|order| {
        let original = |x: usize, y: usize| match (x.min(y), x.max(y)) {
            (0, 1) => s[0],
            (1, 2) => s[1],
            _ => s[2],
        };
        let inv = |a: usize, b: usize| original(order[a], order[b]);
        let m = order.map(|l| masses[l]);
        [
            (inv(0, 1), 0.0),
            (inv(1, 2), 0.0),
            (inv(0, 2), 0.0),
            m[0],
            m[1],
            m[2],
            mu,
        ]
    })
    .collect()
}

fn cases() -> Vec<Point> {
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
            out.extend(relabellings(masses, s));
        }
    }
    out
}

#[test]
fn small_momentum_c0_is_accurate_in_binary64_and_double_double() {
    let mut evaluators = Evaluators::new();
    let mut summary = Summary::default();
    for case in cases() {
        let outcome = evaluators.check(&case);
        summary.add(&case, &outcome);
    }
    summary.assert_accurate("fixed mass sets");
    // Within its orders the series serves every fixed point.
    assert!(summary.smallest_declined_rho_binary64 > ORDERS_BINARY64);
    assert!(summary.smallest_declined_rho_double_double > ORDERS_DOUBLE_DOUBLE);
    // The general expression is what the small-momentum path repairs.
    assert!(summary.worst_general_binary64 > 1e-6);
}

/// Squared masses `m` with invariants `(0.4, −0.35, 0.25) (Re m_i² + Re m_j²) ρ`.
fn scaled(m: [(f64, f64); 3], target: f64) -> Point {
    let den = [m[0].0 + m[1].0, m[1].0 + m[2].0, m[0].0 + m[2].0];
    [
        (0.4 * den[0] * target, 0.0),
        (-0.35 * den[1] * target, 0.0),
        (0.25 * den[2] * target, 0.0),
        m[0],
        m[1],
        m[2],
        (7.3, 0.0),
    ]
}

fn real(x: f64) -> (f64, f64) {
    (x, 0.0)
}

/// The cases of the review of #21: a cluster of two masses with a nearby,
/// nonzero third mass (A to E), two tiny masses whose expansion overflowed
/// (F) and wide widths beyond the cluster expansion's reach (G).
#[test]
fn review_of_21_cases() {
    let mut evaluators = Evaluators::new();
    let a = [real(1.0), real(0.66655), real(0.45)];
    let sm = [real(8315.0), real(6467.2), real(15625.0)];
    let d = [
        real(0.0),
        real(-11.506190046397727),
        real(-14.539770758302355),
        real(675.0403603627077),
        real(445.21451896598046),
        real(322.5358104582457),
        real(7.3),
    ];
    let f = [
        real(0.0),
        real(-5.107859805229796e-24),
        real(0.0),
        real(0.037034422757829706),
        real(6.334868353316832e-9),
        real(5.15208355328473e-9),
        real(7.3),
    ];
    let g = [
        real(2.4365910182777774e-7),
        real(-2.232458914626797e-7),
        real(-3.042667720433323e-7),
        (4084365.994061994, -234743.79550801977),
        (2835925.136725336, -1027301.8922336465),
        real(2774783.2466520304),
        real(7.3),
    ];
    // At ρ = 1/30 and 1/500 the series may decline (the general expression
    // keeps 13 digits); below, it serves.
    for (name, point) in [
        ("A", scaled(a, 1.0 / 30.0)),
        ("C", scaled([real(1.0), real(0.7), real(0.45)], 1.0 / 30.0)),
        ("D", d),
        ("E", scaled(sm, 1.0 / 30.0)),
        ("B", scaled(a, 1.0 / 500.0)),
    ] {
        let outcome = evaluators.check(&point);
        println!("{name}: {outcome:?}");
        assert!(
            outcome.binary64 < 1e-12,
            "{name}: binary64 {:e}",
            outcome.binary64
        );
        assert!(
            outcome.double_double < 1e-27,
            "{name}: double-double {:e}",
            outcome.double_double
        );
        assert!(
            !outcome.served_binary64 || outcome.binary64 < SERVED_BINARY64,
            "{name}: served binary64 {:e}",
            outcome.binary64
        );
    }
    for (name, point) in [
        ("A at 1e-6", scaled(a, 1e-6)),
        ("E at 1e-6", scaled(sm, 1e-6)),
        ("F", f),
        ("G", g),
    ] {
        let outcome = evaluators.check(&point);
        println!("{name}: {outcome:?}");
        assert!(
            outcome.binary64 < SERVED_BINARY64,
            "{name}: binary64 {:e}",
            outcome.binary64
        );
        assert!(outcome.served_binary64, "{name}: binary64 declined");
        // Double-double may decline the comparable masses at ρ = 1e-6 (its
        // expansions are longer); F and G it serves.
        let bound = if outcome.served_double_double {
            SERVED_DOUBLE_DOUBLE
        } else {
            1e-24
        };
        assert!(
            outcome.double_double < bound,
            "{name}: double-double {:e}",
            outcome.double_double
        );
        if name == "F" || name == "G" {
            assert!(
                outcome.served_double_double,
                "{name}: double-double declined"
            );
        }
    }
}

/// The cases of the second review of #21, which the series must serve:
/// two light masses where the general expression keeps 1e-7 (A), a far mass
/// 4 times lighter than a close pair (B) and masses about 1.5 times apart at
/// binary64's highest orders (C), and points on both sides of the far-mass
/// separation (8) and of binary64's closeness (a factor 2).
#[test]
fn second_review_of_21_cases() {
    let mut evaluators = Evaluators::new();
    let pair = |s: f64, close: f64, far: f64| {
        [
            real(0.0),
            real(s),
            real(0.0),
            real(1.0),
            real(far),
            real(close),
            real(7.3),
        ]
    };
    // (1, x, 0) with p1² only at ρ = 0.037.
    let apart = |x: f64| {
        [
            real(-0.037 * (1.0 + x)),
            real(0.0),
            real(0.0),
            real(1.0),
            real(x),
            real(0.0),
            real(7.3),
        ]
    };
    let mut cases: Vec<(String, Point, bool, bool)> = vec![
        (
            "A1".into(),
            [
                real(0.0),
                real(2.7154560641709367e-12),
                real(0.0),
                real(0.002833936460072468),
                real(4.960376656625851e-11),
                real(6.391638675976065e-11),
                real(7.3),
            ],
            true,
            false,
        ),
        (
            "A2".into(),
            [
                real(0.0),
                real(7.725726775100916e-12),
                real(0.0),
                real(0.0009112614124386656),
                real(1.8375402578626648e-10),
                real(1.5832519539193448e-11),
                real(7.3),
            ],
            true,
            false,
        ),
        (
            "B1".into(),
            pair(-0.0025784667908368836, 0.9999732, 0.2498547556760716),
            true,
            true,
        ),
        (
            "B3".into(),
            pair(-0.046237225960014644, 0.9998, 0.2498547556760716),
            true,
            false,
        ),
        (
            "C1".into(),
            [
                real(-0.09274666666666667),
                real(0.0),
                real(0.0),
                real(1.0),
                real(1.5066666666666666),
                real(0.0),
                real(7.3),
            ],
            true,
            false,
        ),
        (
            "C2".into(),
            [
                real(-0.009335604456199737),
                real(0.0),
                real(0.0),
                (0.1171875, -0.048828125),
                real(0.126953125),
                (0.1171875, -0.048828125),
                real(7.3),
            ],
            true,
            false,
        ),
        (
            "C3".into(),
            [
                real(-3.5837784979204557),
                real(0.0),
                real(0.0),
                real(48.748828181845894),
                (48.748828181845894, -17.609712998081104),
                real(48.748828181845894),
                real(7.3),
            ],
            true,
            false,
        ),
    ];
    for q in [7.9f64, 8.1] {
        for far in [1.0 / q, q] {
            cases.push((
                format!("far mass {far:.4} (double-double)"),
                pair(-0.0025784667908368836 * far.max(1.0), 0.9999732, far),
                true,
                true,
            ));
            cases.push((
                format!("far mass {far:.4} (binary64)"),
                pair(-0.025 * far.max(1.0), 0.9998, far),
                true,
                false,
            ));
        }
    }
    for x in [1.49, 1.51, 1.99, 2.01] {
        cases.push((format!("(1, {x}, 0)"), apart(x), true, false));
    }
    assert_cases(&mut evaluators, cases);
}

/// Named points and the carriers that must serve them; every point must be
/// accurate on the route it takes (a declined point: the general expression,
/// for double-double within the series' orders).
fn assert_cases(evaluators: &mut Evaluators, cases: Vec<(String, Point, bool, bool)>) {
    for (name, point, binary64, double_double) in cases {
        let outcome = evaluators.check(&point);
        println!("{name}: rho {:.3e} {outcome:?}", rho(&point));
        assert!(
            !binary64 || outcome.served_binary64,
            "{name}: binary64 declined"
        );
        assert!(
            !double_double || outcome.served_double_double,
            "{name}: double-double declined"
        );
        let bound = if outcome.served_binary64 {
            SERVED_BINARY64
        } else {
            1e-12
        };
        assert!(
            outcome.binary64 < bound,
            "{name}: binary64 {:e}",
            outcome.binary64
        );
        if rho(&point) <= ORDERS_DOUBLE_DOUBLE {
            let bound = if outcome.served_double_double {
                SERVED_DOUBLE_DOUBLE
            } else {
                1e-27
            };
            assert!(
                outcome.double_double < bound,
                "{name}: double-double {:e}",
                outcome.double_double
            );
        }
    }
}

/// The third review of #21: a nonzero squared mass far below the largest (D1 to
/// D7, down to 1e-34 of it) must be served at the bounds, with onset sweeps
/// and points on both sides of the rescaling threshold `2^-40`.
#[test]
fn third_review_of_21_cases() {
    let mut evaluators = Evaluators::new();
    let point = |s: [f64; 3], m: [f64; 3]| {
        [
            real(s[0]),
            real(s[1]),
            real(s[2]),
            real(m[0]),
            real(m[1]),
            real(m[2]),
            real(7.3),
        ]
    };
    let mut cases: Vec<(String, Point, bool, bool)> = vec![
        (
            "D1".into(),
            point([-6e-34, 0.0, 0.0], [1e-32, 2e-32, 1.0]),
            true,
            false,
        ),
        (
            "D2".into(),
            point([-1.167e-25, 0.0, 0.0], [1e-24, 2e-24, 1.0]),
            true,
            false,
        ),
        (
            "D3".into(),
            point(
                [-1.0423410757830895e-29, 0.0, 0.0],
                [
                    8.348713658487341e-29,
                    1.8483749301695216e-28,
                    0.04278863316278276,
                ],
            ),
            true,
            false,
        ),
        (
            "D4".into(),
            point([-0.0023, 0.0, 0.0], [1.0, 0.0, 1e-30]),
            true,
            true,
        ),
        (
            "D5".into(),
            point([-0.0023, 0.0, 0.0], [1.0, 0.0, 1e-34]),
            true,
            true,
        ),
        (
            "D6".into(),
            point(
                [-0.0, -0.0, 4.4797187706966715e-22],
                [
                    8.7186308281056e-20,
                    11957381.75392906,
                    1.4348713837390516e-19,
                ],
            ),
            true,
            true,
        ),
        (
            "D7".into(),
            point(
                [-0.0, 1.6683316122991433e-36, 0.0],
                [
                    1.9700872362686774e-6,
                    2.5629803088749135e-34,
                    5.6057353225659384e-34,
                ],
            ),
            true,
            true,
        ),
    ];
    // The onsets: (1, 0, r) at ρ = 2.3e-3 (double-double) and (r, 2r, 1) at
    // ρ = 0.039 (binary64), both with p1² only, and at ρ = 1e-6.
    for r in [1e-12, 1e-20, 1e-24, 1e-26, 1e-28, 1e-30, 1e-34] {
        cases.push((
            format!("(1, 0, {r:e})"),
            point([-0.0023, 0.0, 0.0], [1.0, 0.0, r]),
            true,
            true,
        ));
        cases.push((
            format!("(1, 0, {r:e}) at 1e-6"),
            point([-1e-6, 0.0, 0.0], [1.0, 0.0, r]),
            true,
            true,
        ));
    }
    let threshold = 2f64.powi(-40);
    for r in [
        1e-12,
        1e-22,
        1e-24,
        1e-32,
        0.99 * threshold,
        1.01 * threshold,
    ] {
        cases.push((
            format!("({r:e}, 2r, 1)"),
            point([-0.039 * 3.0 * r, 0.0, 0.0], [r, 2.0 * r, 1.0]),
            true,
            false,
        ));
        cases.push((
            format!("({r:e}, 2r, 1) at 1e-6"),
            point([-1e-6 * 3.0 * r, 0.0, 0.0], [r, 2.0 * r, 1.0]),
            true,
            true,
        ));
    }
    // Across the threshold with a tiny mass on its own and a heavy far node of
    // two propagators.
    for r in [0.99 * threshold, 1.01 * threshold] {
        cases.push((
            format!("(1, 1, {r:e})"),
            point([-0.001, 0.0, -0.001], [1.0, 1.0, r]),
            true,
            true,
        ));
    }
    assert_cases(&mut evaluators, cases);
}

/// xorshift64, so that the points are reproducible without an RNG crate.
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

    fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.uniform()
    }

    fn log_range(&mut self, a: f64, b: f64) -> f64 {
        (a.ln() + (b.ln() - a.ln()) * self.uniform()).exp()
    }

    fn pick(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Squared masses of a class (times a random scale), with a nonzero third
/// mass wherever the class allows it.
fn random_masses(rng: &mut Rng, class: usize) -> [(f64, f64); 3] {
    let mut m = match class {
        // Two masses about a third apart, the third nearby.
        0 => {
            let d = rng.range(0.25, 0.4);
            [
                real(1.0),
                real(1.0 - d),
                real((1.0 - d) * rng.range(0.4, 0.8)),
            ]
        }
        // A close pair and a nonzero far mass anywhere.
        1 => {
            let d = rng.log_range(1e-6, 0.33);
            [real(1.0), real(1.0 - d), real(rng.log_range(1e-3, 30.0))]
        }
        // Two light masses and a heavy one.
        2 => [
            real(rng.log_range(1e-8, 0.1)),
            real(rng.log_range(1e-8, 0.1)),
            real(1.0),
        ],
        // A massless propagator and two masses up to a factor 5 apart.
        3 => [real(0.0), real(1.0), real(rng.log_range(0.2, 5.0))],
        // Two tiny masses, 1e-34 to 1e-12 of the third, up to a factor 3 apart.
        5 => {
            let a = 10f64.powf(rng.range(-34.0, -12.0));
            [real(1.0), real(a), real(a * rng.log_range(0.3, 3.0))]
        }
        // Uniform masses.
        _ => [
            real(rng.uniform()),
            real(rng.uniform()),
            real(rng.uniform()),
        ],
    };
    // Widths up to 0.6 m on some points.
    if rng.pick(3) == 0 {
        for z in m.iter_mut() {
            if z.0 > 0.0 && rng.uniform() < 0.8 {
                z.1 = -z.0 * rng.range(0.0, 0.6);
            }
        }
    }
    let scale = 10f64.powf(rng.range(-6.0, 8.0));
    m.map(|z| (z.0 * scale, z.1 * scale))
}

/// Invariants with ρ = `target`: one, two or three nonzero, random signs.
fn random_invariants(rng: &mut Rng, m: &[(f64, f64); 3], target: f64) -> [f64; 3] {
    let mut w = [rng.uniform(), rng.uniform(), rng.uniform()];
    match rng.pick(5) {
        0 => w = [w[0], 0.0, 0.0],
        1 => w = [0.0, w[1], 0.0],
        2 => w[rng.pick(3)] = 0.0,
        _ => {}
    }
    for x in w.iter_mut() {
        if rng.uniform() < 0.5 {
            *x = -*x;
        }
    }
    let den = [m[0].0 + m[1].0, m[1].0 + m[2].0, m[0].0 + m[2].0];
    let raw: f64 = w.iter().map(|x| x.abs()).sum();
    if raw == 0.0 {
        return [0.0; 3];
    }
    let t = target / raw;
    [w[0] * den[0] * t, w[1] * den[1] * t, w[2] * den[2] * t]
}

/// Random three-mass points in the domain (fixed seed), ρ log-uniform from
/// 1e-16 to 1/16, in a random propagator order.
#[test]
fn random_three_mass_points() {
    let classes = [
        "pair a third apart",
        "close pair, far mass",
        "two light",
        "(0, m, M)",
        "uniform",
        "two tiny masses",
    ];
    let mut rng = Rng(0x5EED_C0C0_2026_1119);
    let mut evaluators = Evaluators::new();
    for (class, name) in classes.iter().enumerate() {
        let mut summary = Summary::default();
        while summary.points < 300 {
            let m = random_masses(&mut rng, class);
            let target = rng.log_range(1e-16, 1.0 / 16.0);
            let s = random_invariants(&mut rng, &m, target);
            let all = relabellings(m, s);
            let point = all[rng.pick(6)];
            if rho(&point) <= 1.0 / 16.0 {
                summary.add(&point, &evaluators.check(&point));
            }
        }
        summary.assert_accurate(name);
    }
}

/// The route and the value do not depend on the propagator order.
#[test]
fn propagator_order_does_not_change_the_route_or_the_value() {
    let mut rng = Rng(0xC0FF_EE00_2026_1119);
    let tolerance = binary64_tolerance();
    let tolerance_dd = double_double_tolerance();
    for _ in 0..400 {
        let class = rng.pick(5);
        let m = random_masses(&mut rng, class);
        let target = rng.log_range(1e-14, 1.0 / 24.0);
        let s = random_invariants(&mut rng, &m, target);
        let points = relabellings(m, s);
        let values: Vec<Option<Complex<f64>>> = points
            .iter()
            .map(|p| {
                let input: Vec<_> = p.iter().map(|&(re, im)| Complex::new(re, im)).collect();
                c0_small_momentum(&input, &tolerance)
            })
            .collect();
        let values_dd: Vec<Option<Complex<DoubleFloat>>> = points
            .iter()
            .map(|p| {
                let input: Vec<_> = p
                    .iter()
                    .map(|&(re, im)| Complex::new(DoubleFloat::from(re), DoubleFloat::from(im)))
                    .collect();
                c0_small_momentum(&input, &tolerance_dd)
            })
            .collect();
        let first = values[0];
        for value in &values {
            assert_eq!(
                value.is_some(),
                first.is_some(),
                "{:?}: binary64 route",
                points[0]
            );
            if let (Some(a), Some(b)) = (value, first) {
                let spread = (a.re - b.re).hypot(a.im - b.im) / b.re.hypot(b.im);
                assert!(
                    spread < 4e-15,
                    "{:?}: binary64 spread {spread:e}",
                    points[0]
                );
            }
        }
        let first = &values_dd[0];
        for value in &values_dd {
            assert_eq!(
                value.is_some(),
                first.is_some(),
                "{:?}: double-double route",
                points[0]
            );
            if let (Some(a), Some(b)) = (value, first) {
                let dr = wide(&a.re) - wide(&b.re);
                let di = wide(&a.im) - wide(&b.im);
                let size = wide(&b.re).to_f64().hypot(wide(&b.im).to_f64());
                let spread = dr.to_f64().hypot(di.to_f64()) / size;
                assert!(
                    spread < 1e-29,
                    "{:?}: double-double spread {spread:e}",
                    points[0]
                );
            }
        }
    }
}

/// Outside the series' domain the evaluator is the general expression, bit
/// for bit: ρ above 1/16, two massless propagators (collinear or soft
/// singularities), physical invariants, and an invariant with an imaginary
/// part.
#[test]
fn outside_the_series_domain_c0_is_the_general_expression() {
    let mut double = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
    let mut general = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
    general.set_small_momentum(false);
    let (m, mu) = (6467.2, 1.0e6);
    let points: [[(f64, f64); 7]; 6] = [
        [-900.0, 0.0, 300.0, m, m, m, mu].map(real),
        [0.0, 0.0, -1.0, 0.0, 0.0, m, mu].map(real),
        [-1.0e-3, 0.0, 0.0, m, 0.0, 0.0, mu].map(real),
        [29929.0, 0.0, 1.0e5, 29929.0, 29929.0, 0.0, mu].map(real),
        [-2.0e4, 1.0e3, -5.0e3, m, 2.0 * m, 3.0 * m, mu].map(real),
        [
            (-1.0e-3, 1.0e-4),
            real(0.0),
            real(0.0),
            real(m),
            real(m),
            real(m),
            real(mu),
        ],
    ];
    for point in points {
        let input: Vec<_> = point.iter().map(|&(re, im)| Complex::new(re, im)).collect();
        assert!(c0_small_momentum(&input, &binary64_tolerance()).is_none());
        if point[0].1 != 0.0 {
            continue; // The evaluator requires real invariants.
        }
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
