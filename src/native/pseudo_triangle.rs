//! D0 with three massive propagators pairwise at their pseudo-thresholds (a
//! *pseudo-threshold triangle*), through an exact closed form in bubbles.
//!
//! Three real positive masses `m_j = v_j²` (`j` in the triangle `T`) whose
//! invariants all sit at their pseudo-thresholds, `p_jk = (v_j - v_k)²` (equal
//! masses at `p_jk = 0` among them: two collinear massless legs between
//! equal-mass lines), give the modified Cayley matrix the rank-one block
//! `Y_jk = m_j + m_k - p_jk = 2 v_j v_k`. The Feynman-parameter denominator
//! `F = x^T Y x / 2` is then `m_o x_o² + x_o Σ_j Y_oj x_j + (Σ_j v_j x_j)²`,
//! with `o` the fourth propagator (mass `m_o`, real or complex, possibly zero).
//! `det Y = 0`, and every finite form of the box divides by a power of it:
//! the chart and DNS forms of `box_finite_four_masses` take double roots at
//! `r = 1` and dilogarithms at their branch points, the contour box does not
//! support a vanishing channel, and the sector code returned values up to
//! `1e77` off at 1024 bits (MadGraph's stress points of `u u~ > z g g` and
//! `u u~ > z z g g`: two exactly collinear gluons on a b-quark loop).
//!
//! With the projective weights `x_o + Σ_j (v_j/κ) x_j = 1` (Cheng-Wu, any
//! `κ > 0`), `z_j = v_j x_j / κ`, `x_o = t`, `z = (1 - t) y` (`y` on the
//! triangle), the denominator is `A(t) - t (1 - t) Σ_j L_j y_j`, linear on
//! the triangle, with `A(t) = t m_o + (1 - t) κ²` and
//! `L_j = (v_j - κ) (m_o / v_j - κ) + (κ / v_j) p_oj`. Hermite-Genocchi turns
//! the triangle integral into a second divided difference, and the `t`
//! integral of each term is the bubble `B0(L; m_o, κ²)` (finite part; its
//! constant and linear terms in `L` drop out of the difference):
//!
//! `D0 = κ / (v_1 v_2 v_3) · Δ²[L B0(L; m_o, κ²)](L_1, L_2, L_3)`.
//!
//! For equal masses `κ = v` makes `L_j = p_oj`: `D0 = Δ²[p B0(p; m_o, m)] / m`.
//! With all six pairs at their pseudo-thresholds `F = (Σ_i v_i x_i)²` and
//! `D0 = 1 / (6 v_0 v_1 v_2 v_3)`.
//!
//! The case is decided by exact predicates on the represented inputs
//! ([NativeFloat::exact_rational]), after a rounded prefilter. The divided
//! difference cancels by `(scale / gap)²` where two `L_j` are close (MadGraph's
//! three `p_oj` agree to 1e-7), so the bubbles are evaluated in `Float` with
//! three times those bits on top of the carrier's (numerica tracks the
//! cancellation of each gap and that of the difference as lost precision), and
//! the result is rounded once.
//! `L_j` closer than `2^-(bits + 16)` of the scale are spread symmetrically
//! about their mean to that distance: the divided difference is symmetric, so
//! the change is second order, `2^-(2 bits + 32)`. A cluster next to a branch
//! point of the bubble (`(√m_o ± κ)²`, where the difference is not smooth),
//! complex `L_j` (a complex `m_o` with unequal triangle masses), and an
//! infrared-divergent box (massless `m_o` between two on-shell legs) are left
//! to the sector code, as is the equal-mass vacuum box, which it evaluates
//! exactly.
#![cfg(feature = "generated-evaluators")]
use super::NativeEvaluator;
use super::primitives::{C, NativeFloat, add, inverse, mul};
use crate::ScalarIntegral;
use std::cell::RefCell;
use symbolica::domains::float::{Float, Numeric, Real, RealLike, SingleFloat};
use symbolica::domains::rational::Rational;

/// `(invariant, propagator i, propagator j)`: the legs `p1..p4` and the
/// channels `s12 = p_02`, `s23 = p_13`.
const PAIRS: [(usize, usize, usize); 6] = [
    (0, 0, 1),
    (1, 1, 2),
    (2, 2, 3),
    (3, 0, 3),
    (4, 0, 2),
    (5, 1, 3),
];

/// The three propagators other than `o`, for each `o`.
const TRIANGLES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// The three pairs of a triangle.
fn pairs([a, b, c]: [usize; 3]) -> [(usize, usize); 3] {
    [(a, b), (b, c), (a, c)]
}

fn invariant(i: usize, j: usize) -> usize {
    PAIRS
        .iter()
        .find(|&&(_, a, b)| (a, b) == (i.min(j), i.max(j)))
        .map(|&(k, _, _)| k)
        .expect("distinct propagators")
}

/// `(a - b)² - 2 p (a + b) + p²`, the Källén function of `p` between the
/// masses `a` and `b`, exactly.
fn kallen(p: &Rational, a: &Rational, b: &Rational) -> Rational {
    let difference = a - b;
    let twice = p + p;
    &(&(&difference * &difference) - &(&twice * &(a + b))) + &(p * p)
}

/// [kallen] in the carrier, for the prefilter.
fn rounded_kallen<T: NativeFloat>(p: &T, a: &T, b: &T) -> T {
    let difference = a.clone() - b;
    let twice = p.clone() + p;
    difference.clone() * &difference - &(twice * &(a.clone() + b)) + &(p.clone() * p)
}

/// The represented inputs as exact rationals, `(re, im)` per argument.
struct Exact {
    re: Vec<Rational>,
    im: Vec<Rational>,
}

impl Exact {
    fn of<T: NativeFloat>(input: &[C<T>]) -> Option<Self> {
        let mut re = Vec::with_capacity(10);
        let mut im = Vec::with_capacity(10);
        for z in &input[..10] {
            re.push(z.re.exact_rational()?);
            im.push(z.im.exact_rational()?);
        }
        Some(Self { re, im })
    }

    fn mass(&self, line: usize) -> &Rational {
        &self.re[6 + line]
    }

    fn real_mass(&self, line: usize) -> bool {
        self.im[6 + line].is_zero()
    }

    fn positive_mass(&self, line: usize) -> bool {
        self.real_mass(line) && self.mass(line) > &Rational::zero()
    }

    /// Whether the pair `(i, j)` (real positive masses) sits at its
    /// pseudo-threshold: a vanishing Källén function and `m_i + m_j >= p_ij`.
    fn pseudo(&self, i: usize, j: usize) -> bool {
        let p = &self.re[invariant(i, j)];
        let (a, b) = (self.mass(i), self.mass(j));
        kallen(p, a, b).is_zero() && &(a + b) >= p
    }
}

/// Whether the pair `(i, j)` of real positive masses may sit at its
/// pseudo-threshold, for the exact test: `0 <= p <= max(a, b)` (since
/// `p = (√a - √b)²`), and a rounded Källén function within
/// `2^-40 (a + b + p)²`.
fn near_pseudo<T: NativeFloat>(input: &[C<T>], i: usize, j: usize, two40: &T) -> bool {
    let (p, a, b) = (
        &input[invariant(i, j)].re,
        &input[6 + i].re,
        &input[6 + j].re,
    );
    if p.negative_sign() && !p.is_zero() || (p > a && p > b) {
        return false;
    }
    let lambda = rounded_kallen(p, a, b);
    let scale = a.clone() + b + p;
    lambda.norm() * two40 <= scale.clone() * &scale
}

/// The propagator `o` outside a pseudo-threshold triangle, if the box has one.
fn outside<T: NativeFloat>(input: &[C<T>]) -> Option<(usize, Exact)> {
    let positive: [bool; 4] = core::array::from_fn(|line| {
        let z = &input[6 + line];
        z.im.is_zero() && z.re > z.re.zero()
    });
    if positive.iter().filter(|&&p| p).count() < 3 {
        return None;
    }
    // The prefilter of each pair of real positive masses, once.
    let two40 = input[0].re.from_i64(1 << 40);
    let near: [bool; 6] = core::array::from_fn(|k| {
        let (_, i, j) = PAIRS[k];
        positive[i] && positive[j] && near_pseudo(input, i, j, &two40)
    });
    let passes = |o: usize| {
        pairs(TRIANGLES[o])
            .iter()
            .all(|&(i, j)| near[invariant(i, j)])
    };
    if !(0..4).any(passes) {
        return None;
    }
    let exact = Exact::of(input)?;
    let o = (0..4).filter(|&o| passes(o)).find(|&o| {
        pairs(TRIANGLES[o])
            .iter()
            .all(|&(i, j)| exact.positive_mass(i) && exact.positive_mass(j) && exact.pseudo(i, j))
    })?;
    Some((o, exact))
}

/// A `Float` at `bits` from an exact rational.
fn wide(bits: u32, value: &Rational) -> Float {
    Float::new(bits).from_rational(value)
}

/// `B0(p; a, b)` (finite part, `mu² = b`) at `bits`, from a cached evaluator.
fn bubble(bits: u32, p: &Float, a: &C<Float>, b: &Float) -> Option<C<Float>> {
    thread_local! {
        static B0: RefCell<Option<NativeEvaluator<Float>>> = const { RefCell::new(None) };
    }
    B0.with(|cell| {
        let mut cached = cell.borrow_mut();
        if cached
            .as_ref()
            .is_none_or(|evaluator| evaluator.binary_precision() != bits)
        {
            *cached = Some(NativeEvaluator::with_binary_precision(ScalarIntegral::B0, bits).ok()?);
        }
        let evaluator = cached.as_mut()?;
        let zero = Float::new(bits);
        let input = [
            C::new(p.clone(), zero.clone()),
            a.clone(),
            C::new(b.clone(), zero.clone()),
            C::new(b.clone(), zero.clone()),
        ];
        let mut output: [C<Float>; 3] =
            core::array::from_fn(|_| C::new(zero.clone(), zero.clone()));
        evaluator.evaluate(&input, &mut output).ok()?;
        let value = output[0].clone();
        (value.re.is_finite() && value.im.is_finite()).then_some(value)
    })
}

/// The triangle's data at `bits`: the roots `v_j`, `κ = max v_j`, and the
/// `L_j` (real), from the exact inputs.
struct Triangle {
    lines: [usize; 3],
    equal: bool,
}

impl Triangle {
    fn at(&self, exact: &Exact, o: usize, bits: u32) -> ([Float; 3], Float, [Float; 3]) {
        let v = self.lines.map(|l| wide(bits, exact.mass(l)).sqrt());
        let kappa = v
            .iter()
            .cloned()
            .reduce(|a, b| if b > a { b } else { a })
            .expect("three lines");
        let m_o = wide(bits, exact.mass(o));
        let l = core::array::from_fn(|n| {
            let p = wide(bits, &exact.re[invariant(o, self.lines[n])]);
            if self.equal {
                p
            } else {
                // (v - κ)(m_o / v - κ) + (κ / v) p
                let ratio = kappa.clone() / &v[n];
                (v[n].clone() - &kappa) * (m_o.clone() / &v[n] - &kappa) + ratio * &p
            }
        });
        (v, kappa, l)
    }
}

/// The three `L` after spreading clusters closer than `delta` about their
/// mean (module documentation), and the clusters' centres.
fn spread(l: [Float; 3], delta: &Float) -> ([Float; 3], Vec<Float>) {
    let mut order = [0usize, 1, 2];
    order.sort_by(|&a, &b| {
        l[a].partial_cmp(&l[b])
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    let [a, b, c] = order;
    let three = delta.clone() + delta + delta;
    let mut out = l.clone();
    let mut centres = Vec::new();
    if l[c].clone() - &l[a] < three {
        let centre = (l[a].clone() + &l[b] + &l[c]) / &l[a].from_i64(3);
        out[a] = centre.clone() - delta;
        out[b] = centre.clone();
        out[c] = centre.clone() + delta;
        centres.push(centre);
    } else {
        for (x, y) in [(a, b), (b, c)] {
            if l[y].clone() - &l[x] < *delta {
                let centre = (l[x].clone() + &l[y]) / &l[x].from_i64(2);
                out[x] = centre.clone() - delta;
                out[y] = centre.clone() + delta;
                centres.push(centre);
            }
        }
    }
    (out, centres)
}

/// The D0 `[finite, pole, double pole]` of a box with a pseudo-threshold
/// triangle (module documentation) for a carrier of `bits` bits (the
/// carrier's own for binary64 and double-double), or `None`.
#[inline]
pub(crate) fn d0<T: NativeFloat>(input: &[C<T>], bits: u32) -> Option<[C<T>; 3]> {
    // At least three real positive masses, inline; the rest out of line.
    let positive = |line: usize| {
        let z = &input[6 + line];
        z.im.is_zero() && z.re > z.re.zero()
    };
    if (positive(0) as u8) + (positive(1) as u8) + (positive(2) as u8) + (positive(3) as u8) < 3 {
        return None;
    }
    triangle_d0(input, bits)
}

#[inline(never)]
fn triangle_d0<T: NativeFloat>(input: &[C<T>], bits: u32) -> Option<[C<T>; 3]> {
    let (o, exact) = outside(input)?;
    let lines = TRIANGLES[o];
    let target = T::FIXED_BITS.unwrap_or(bits);
    let zero = input[0].re.zero();
    let to_carrier = |x: &Float| -> Option<T> { Some(zero.from_rational(&x.try_to_rational()?)) };
    let finite = |value: C<Float>| -> Option<[C<T>; 3]> {
        let z = C::new(zero.clone(), zero.clone());
        Some([
            C::new(to_carrier(&value.re)?, to_carrier(&value.im)?),
            z.clone(),
            z,
        ])
    };
    let invariants_zero = exact.re[..6].iter().all(Rational::is_zero);
    let masses_equal =
        (1..4).all(|l| exact.re[6 + l] == exact.re[6] && exact.im[6 + l] == exact.im[6]);
    if invariants_zero && masses_equal {
        // The vacuum box, which the sector code evaluates exactly.
        return None;
    }
    let equal = lines.iter().all(|&l| exact.mass(l) == exact.mass(lines[0]));
    let o_real = exact.real_mass(o);
    if !o_real && !equal {
        return None;
    }
    let o_zero = o_real && exact.mass(o).is_zero();
    if o_zero {
        // A massless line between two on-shell legs: infrared divergent.
        let on_shell = |l: usize| exact.re[invariant(o, l)] == *exact.mass(l);
        if on_shell((o + 1) % 4) && on_shell((o + 3) % 4) {
            return None;
        }
    }
    let triangle = Triangle { lines, equal };
    // All six pairs at their pseudo-thresholds: F = (Σ v_i x_i)².
    if exact.positive_mass(o) && lines.iter().all(|&l| exact.pseudo(o, l)) {
        let bits = target + 64;
        let v = (0..4)
            .map(|l| wide(bits, exact.mass(l)).sqrt())
            .reduce(|a, b| a * b)?;
        let value = (v * &Float::with_val(bits, 6.0)).inv();
        return finite(C::new(value, Float::new(bits)));
    }
    // The scale and the gaps at a provisional precision.
    let provisional = target + 64;
    let (_, kappa, l) = triangle.at(&exact, o, provisional);
    let m_o_re = wide(provisional, exact.mass(o));
    let m_o_im = wide(provisional, &exact.im[6 + o]);
    let mut scale = kappa.clone() * &kappa;
    for x in l.iter().chain([&m_o_re, &m_o_im]) {
        if x.norm() > scale {
            scale = x.norm();
        }
    }
    let unit =
        |bits: u32, k: u32| Float::with_val(bits, 2).powf(&Float::with_val(bits, -(k as f64)));
    let delta = scale.clone() * &unit(provisional, target + 16);
    let branches: Vec<Float> = if o_real {
        let root = m_o_re.clone().sqrt();
        vec![
            (root.clone() + &kappa) * &(root.clone() + &kappa),
            (root.clone() - &kappa) * &(root - &kappa),
        ]
    } else {
        Vec::new()
    };
    let (spread_l, centres) = spread(l, &delta);
    // A cluster next to a branch point of the bubble: not smooth there.
    let near = scale.clone() * &unit(provisional, target / 2 + 8);
    if centres
        .iter()
        .any(|c| branches.iter().any(|b| (c.clone() - b).norm() < near))
    {
        return None;
    }
    let mut gap = scale.clone();
    for (i, j) in [(0, 1), (1, 2), (0, 2)] {
        let g = (spread_l[i].clone() - &spread_l[j]).norm();
        if g < gap {
            gap = g;
        }
    }
    // log2(scale / gap) through the logarithm: the ratio reaches 2^(bits + 16).
    let cancellation = ((scale / &gap).log().to_f64() / core::f64::consts::LN_2)
        .ceil()
        .max(0.0) as u32;
    // Numerica tracks the cancellation of each gap `L_j - L_k` (`cancellation`
    // bits) and that of the divided difference (twice as many) as lost
    // precision: three times the bits, on top of the carrier's.
    let bits = (target + 3 * cancellation + 64).div_ceil(64) * 64;
    // At the working precision: the clusters spread by the same rule.
    let (v, kappa, l) = triangle.at(&exact, o, bits);
    let delta = {
        let mut scale = kappa.clone() * &kappa;
        for x in l
            .iter()
            .chain([&wide(bits, exact.mass(o)), &wide(bits, &exact.im[6 + o])])
        {
            if x.norm() > scale {
                scale = x.norm();
            }
        }
        scale * &unit(bits, target + 16)
    };
    let (l, _) = spread(l, &delta);
    let m_o = C::new(wide(bits, exact.mass(o)), wide(bits, &exact.im[6 + o]));
    let kappa2 = kappa.clone() * &kappa;
    let mut sum = C::new(Float::new(bits), Float::new(bits));
    for n in 0..3 {
        let b0 = bubble(bits, &l[n], &m_o, &kappa2)?;
        let g = mul(&C::new(l[n].clone(), Float::new(bits)), &b0);
        let mut denominator = Float::with_val(bits, 1);
        for k in (0..3).filter(|&k| k != n) {
            denominator *= &(l[n].clone() - &l[k]);
        }
        sum = add(
            &sum,
            &mul(&g, &inverse(&C::new(denominator, Float::new(bits)))),
        );
    }
    let prefactor = kappa / &(v[0].clone() * &v[1] * &v[2]);
    finite(mul(&sum, &C::new(prefactor, Float::new(bits))))
}
