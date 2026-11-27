//! The soft three-mass box (OneLOop's sector 16, `box_three_masses_16`) for
//! every carrier, where the generated code takes another labelling.
//!
//! One massless line with both of its legs on their massive partners' mass
//! shells is soft divergent. `box_series` keeps such a box in its labelling
//! even when a channel vanishes (the channel through the soft line and the
//! opposite line is zero for an accidental light-like leg of the
//! process): either channel exchange would move an on-shell pair of the soft
//! line onto a channel, where the finite three-mass formula applies and
//! misses the pole. The generated code and the SymJIT asset predate that
//! test, so the numerical entry points evaluate these points here instead: a
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

/// The massless line of a soft three-mass box with a vanishing channel in
/// the labelling of `point` (`[p1², p2², p3², p4², s12, s23, m1², m2², m3²,
/// m4², mu²]`), exact tests: the points that `box_series` keeps in this
/// labelling and the generated code would not.
pub(crate) fn soft_line<T: NativeFloat>(point: &[C<T>]) -> Option<usize> {
    if point.len() != 11 || !(is_zero(&point[4]) || is_zero(&point[5])) {
        return None;
    }
    let (p, m) = (&point[..6], &point[6..10]);
    (0..4).find(|&i| {
        let [(a, pa), (b, pb)] = LEGS[i];
        is_zero(&m[i])
            && (1..4).all(|k| !is_zero(&m[(i + k) % 4]))
            && on_shell(&p[a], &m[pa])
            && on_shell(&p[b], &m[pb])
    })
}

/// Whether `point` has a massless line and a vanishing invariant, which every
/// [soft_line] point of any labelling has.
#[inline]
pub(crate) fn candidate<T: NativeFloat>(point: &[C<T>]) -> bool {
    point.len() == 11 && point[6..10].iter().any(is_zero) && point[..6].iter().any(is_zero)
}

/// [d0] of `point` in the labelling the numerical entry points evaluate
/// (`crate::labelling::soft_order`). Points without a massless line or a
/// vanishing invariant return `None` at once.
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
    let (r23, _) = r_function(
        &div(&sub(&add(mass_2, mass_3), p2), &mul(&sqrt_2, &sqrt_3)),
        c,
    );
    let (r24, _) = r_function(
        &div(&sub(&add(mass_2, mass_4), p23), &mul(&sqrt_2, &sqrt_4)),
        c,
    );
    let (r34, _) = r_function(
        &div(&sub(&add(mass_3, mass_4), p3), &mul(&sqrt_3, &sqrt_4)),
        c,
    );
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
