//! The soft three-mass box (OneLOop's sector 16, `box_three_masses_16`) for
//! every carrier, where the generated code gets it wrong.
//!
//! One massless line with both of its legs on their massive partners' mass
//! shells is soft divergent. `box_series` keeps such a box in its labelling
//! even when a channel vanishes (the channel through the soft line and the
//! opposite line is zero for an accidental light-like leg of the
//! process): either channel exchange would move an on-shell pair of the soft
//! line onto a channel, where the finite three-mass formula applies and
//! misses the pole. At the pseudo-threshold of one of the sector's three
//! massive pairs (their Källén function vanishes: the soft line's two
//! neighbours, for example equal masses and a light-like channel between
//! them, or one neighbour and the opposite line), `box_three_masses_16` takes
//! the pair's double root `r = 1` exactly: rounding would split it into
//! `1 ± sqrt(rounding)`. For the neighbours, `r24`, the divided differences of
//! merging arguments then lose half their digits (all of them in `Float`); for
//! the other two pairs, `r23` and `r34`, `Float` loses half its digits through
//! its significance tracking. The generated code and the SymJIT asset predate
//! both, so the numerical entry points evaluate these points here instead: a
//! transcription of `box_three_masses_16`, with the sheets of [`Sheet`] and
//! the divided difference [`dd`], branch for branch.
use super::super::primitives::{C, NativeFloat, mul, sqrt};
use super::sheet::{Sheet, dd};
use super::{Constants, add, div, is_zero, neg, negative, physical_sqrt, r_function, real, sub};

/// Legs `(invariant, partner mass)` of each line (as `expressions::SOFT_LEGS`).
const LEGS: [[(usize, usize); 2]; 4] = [
    [(0, 1), (3, 3)],
    [(0, 0), (1, 2)],
    [(1, 1), (2, 3)],
    [(2, 2), (3, 0)],
];

/// The canonical orders of `box_three_masses_ir` (massless line first) for a
/// massless line at position 0 to 3: invariants and masses.
const ORDER: [[usize; 6]; 4] = [
    [0, 1, 2, 3, 4, 5],
    [1, 2, 3, 0, 5, 4],
    [2, 3, 0, 1, 4, 5],
    [3, 0, 1, 2, 5, 4],
];

fn on_shell<T: NativeFloat>(value: &C<T>, mass: &C<T>) -> bool {
    is_zero(&sub(value, mass))
}

/// Whether the Källén function `sum² - 4 ma mb` of two lines with squared
/// masses `ma`, `mb` vanishes exactly, given `sum = ma + mb - s` for their
/// channel `s`: their threshold or pseudo-threshold, where `r + 1/r = ±2` has
/// a double root. (The factor 4 by two exact doublings.)
fn kallen_vanishes<T: NativeFloat>(sum: &C<T>, ma: &C<T>, mb: &C<T>) -> bool {
    let product = mul(ma, mb);
    let double = add(&product, &product);
    is_zero(&sub(&mul(sum, sum), &add(&double, &double)))
}

/// The massless line of a soft three-mass box in the labelling of `point`
/// (`[p1², p2², p3², p4², s12, s23, m1², m2², m3², m4², mu²]`) that the
/// generated code gets wrong, exact tests: a vanishing channel, which
/// `box_series` keeps in this labelling, or one of the sector's three massive
/// pairs at its pseudo-threshold (`sum > 0`, not its threshold): the soft
/// line's neighbours, or one of them and the opposite line.
pub(crate) fn soft_line<T: NativeFloat>(point: &[C<T>]) -> Option<usize> {
    if point.len() != 11 {
        return None;
    }
    let (p, m) = (&point[..6], &point[6..10]);
    let vanishing = is_zero(&p[4]) || is_zero(&p[5]);
    (0..4).find(|&i| {
        let [(a, pa), (b, pb)] = LEGS[i];
        let pseudo_threshold = || {
            // (invariant, line, line): the leg `p[k]` joins lines `k` and
            // `k + 1`; the neighbours' channel is s23 for lines 0 and 2, s12
            // for 1 and 3.
            let (next, opposite, previous) = ((i + 1) % 4, (i + 2) % 4, (i + 3) % 4);
            [
                (next, next, opposite),
                (opposite, opposite, previous),
                (5 - i % 2, pa, pb),
            ]
            .into_iter()
            .any(|(k, ka, kb)| {
                let sum = sub(&add(&m[ka], &m[kb]), &p[k]);
                !negative(&sum.re) && kallen_vanishes(&sum, &m[ka], &m[kb])
            })
        };
        is_zero(&m[i])
            && (1..4).all(|k| !is_zero(&m[(i + k) % 4]))
            && on_shell(&p[a], &m[pa])
            && on_shell(&p[b], &m[pb])
            && (vanishing || pseudo_threshold())
    })
}

/// Whether `point` has exactly one massless line and an invariant on the
/// shell of a massive line, which every [soft_line] point of any labelling
/// has.
#[inline]
pub(crate) fn candidate<T: NativeFloat>(point: &[C<T>]) -> bool {
    if point.len() != 11 {
        return false;
    }
    let masses = &point[6..10];
    masses.iter().filter(|m| is_zero(m)).count() == 1
        && point[..6]
            .iter()
            .any(|p| masses.iter().any(|m| !is_zero(m) && on_shell(p, m)))
}

/// [d0] of `point` in the labelling the numerical entry points evaluate
/// (`crate::labelling::soft_order`). Points without exactly one massless line
/// and an invariant on a massive line's shell return `None` at once.
#[inline]
pub(crate) fn d0_any_labelling<T: NativeFloat>(
    point: &[C<T>],
    c: &Constants<T>,
) -> Option<[C<T>; 3]> {
    if !candidate(point) {
        return None;
    }
    match crate::labelling::soft_order(point) {
        Some(order) => {
            let relabelled: [C<T>; 11] = core::array::from_fn(|k| point[order[k]].clone());
            d0(&relabelled, c)
        }
        None => d0(point, c),
    }
}

/// D0 `[finite, pole, double pole]` of a [soft_line] point.
pub(crate) fn d0<T: NativeFloat>(point: &[C<T>], c: &Constants<T>) -> Option<[C<T>; 3]> {
    let shift = soft_line(point)?;
    let p: [&C<T>; 6] = core::array::from_fn(|k| &point[ORDER[shift][k]]);
    let m: [&C<T>; 4] = core::array::from_fn(|k| &point[6 + (k + shift) % 4]);
    let [finite, pole] = sector_16(p[1], p[2], p[4], p[5], m[1], m[2], m[3], &point[10], c);
    let zero = real(c.one.zero());
    Some([finite, pole, zero])
}

/// `box_three_masses_16` (with `with_box_normalization`, whose double-pole
/// term vanishes): `[finite, pole]`.
#[allow(clippy::too_many_arguments)] // OneLOop sector 16's named invariants.
fn sector_16<T: NativeFloat>(
    p2: &C<T>,
    p3: &C<T>,
    p12: &C<T>,
    p23: &C<T>,
    mass_2: &C<T>,
    mass_3: &C<T>,
    mass_4: &C<T>,
    mu_squared: &C<T>,
    c: &Constants<T>,
) -> [C<T>; 2] {
    // The massive endpoints by magnitude, exchanging the adjacent invariants.
    let size = |z: &C<T>| z.re.clone() * &z.re + &(z.im.clone() * &z.im);
    let swap = negative(&(size(mass_2) - &size(mass_4)));
    let (p2, p3, mass_2, mass_4) = if swap {
        (p3, p2, mass_4, mass_2)
    } else {
        (p2, p3, mass_2, mass_4)
    };
    let mu = sqrt(mu_squared);
    let sqrt_2 = physical_sqrt(mass_2);
    let sqrt_3 = physical_sqrt(mass_3);
    let sqrt_4 = physical_sqrt(mass_4);
    let r13 = div(&sub(mass_3, p12), &mul(&mu, &sqrt_3));
    // The double root 1 at a pair's pseudo-threshold (as the expression's
    // `pair_r`).
    let pair_r = |p: &C<T>, ma: &C<T>, mb: &C<T>, sqrt_a: &C<T>, sqrt_b: &C<T>| {
        let sum = sub(&add(ma, mb), p);
        let x = div(&sum, &mul(sqrt_a, sqrt_b));
        if kallen_vanishes(&sum, ma, mb) && !negative(&x.re) {
            c.one_c.clone()
        } else {
            r_function(&x, c).0
        }
    };
    let r23 = pair_r(p2, mass_2, mass_3, &sqrt_2, &sqrt_3);
    let r24 = pair_r(p23, mass_2, mass_4, &sqrt_2, &sqrt_4);
    let r34 = pair_r(p3, mass_3, mass_4, &sqrt_3, &sqrt_4);
    let [q13, q23, q24, q34] = [&r13, &r23, &r24, &r34].map(|r| Sheet::lower(r.clone()));
    let qproduct = q13.mul(&q23);
    let qsquare = qproduct.mul(&qproduct);
    let qss = qsquare.div(&q24);
    let y = q23.mul(&q34);
    let (y1, y2) = (y.mul(&q24), y.div(&q24));
    let z = q23.div(&q34);
    let (z1, z2) = (z.mul(&q24), z.div(&q24));
    let factor = div(&c.one_c, &mul(&mul(&sqrt_2, &sqrt_4), &sub(p12, mass_3)));
    let log24 = div(&mul(&q24.log_over_one_minus(c), &r24), &add(&c.one_c, &r24));
    let finite = sub(
        &sub(
            &add(
                &mul(&log24, &qss.log(c)),
                &mul(&dd(&q24.mul(&q24), &Sheet::lower(c.one_c.clone()), c), &r24),
            ),
            &mul(&mul(&dd(&y1, &y2, c), &r23), &r34),
        ),
        &div(&mul(&dd(&z1, &z2, c), &r23), &r34),
    );
    [mul(&factor, &finite), neg(&mul(&factor, &log24))]
}
