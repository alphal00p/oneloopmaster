//! Compact exact sectors, sharing all Laurent work in one native call.
//!
//! These are algebraic simplifications of `triangle_massless` and
//! `box_massless_low_case` (including `with_box_normalization`). Selection uses
//! exact zero tests only, after the public input-domain checks. No approximate
//! degeneracy test, finite-width regulator or precision fallback is used.
use super::primitives::{C, NativeFloat, mul};

#[derive(Clone)]
pub(super) struct Constants<T> {
    one: T,
    two: T,
    four: T,
    half: T,
    pi: T,
    pi_squared: T,
}

impl<T: NativeFloat> Constants<T> {
    pub(super) fn new(prototype: &T) -> Self {
        let one = prototype.one();
        let two = prototype.from_i64(2);
        let pi = prototype.pi();
        let pi_squared = pi.clone() * &pi;
        Self {
            half: one.clone() / &two,
            one,
            two,
            four: prototype.from_i64(4),
            pi,
            pi_squared,
        }
    }
}

fn zero<T: NativeFloat>(prototype: &T) -> C<T> {
    C::new(prototype.zero(), prototype.zero())
}

pub(super) fn zeros<T: NativeFloat>(prototype: &T) -> [C<T>; 3] {
    std::array::from_fn(|_| zero(prototype))
}

#[inline]
fn scale<T: NativeFloat>(value: &C<T>, factor: &T) -> C<T> {
    C::new(value.re.clone() * factor, value.im.clone() * factor)
}

// For a nonzero real invariant s, log((-s-i0)/mu²). Separate real
// logarithms avoid overflowing/underflowing the ratio. The sheet is explicit,
// independent of an accidental signed zero from intermediate arithmetic.
#[inline]
fn momentum_log<T: NativeFloat>(s: &T, log_mu: &T, c: &Constants<T>) -> C<T> {
    C::new(
        s.norm().log() - log_mu,
        if *s > s.zero() {
            -c.pi.clone()
        } else {
            s.zero()
        },
    )
}

pub(super) fn massless_triangle<T: NativeFloat>(a: &[C<T>], c: &Constants<T>) -> Option<[C<T>; 3]> {
    if a[3..6]
        .iter()
        .any(|z| !z.re.is_fully_zero() || !z.im.is_fully_zero())
    {
        return None;
    }
    let mut nonzero = a[..3].iter().filter(|z| !z.re.is_fully_zero());
    let Some(s) = nonzero.next().map(|z| &z.re) else {
        return Some(zeros(&c.one));
    };
    let t = nonzero.next().map(|z| &z.re);
    if nonzero.next().is_some() {
        return None;
    }
    let log_mu = a[6].re.log();
    let ls = momentum_log(s, &log_mu, c);
    Some(if let Some(t) = t {
        let lt = momentum_log(t, &log_mu, c);
        let difference = s.clone() - t;
        let pole = if difference.is_fully_zero() {
            C::new(-s.inv(), s.zero())
        } else {
            scale(
                &C::new(lt.re.clone() - &ls.re, lt.im.clone() - &ls.im),
                &difference.inv(),
            )
        };
        let sum = C::new(ls.re + lt.re, ls.im + lt.im);
        let finite = scale(&mul(&pole, &sum), &(-c.half.clone()));
        [finite, pole, zero(&c.one)]
    } else {
        let reciprocal = s.inv();
        // The raw triangle sector's -pi²/(12s) cancels against the
        // +pi²/12 times double-pole normalization of the public master.
        let finite = scale(&mul(&ls, &ls), &c.half);
        [
            scale(&finite, &reciprocal),
            scale(&ls, &(-reciprocal.clone())),
            C::new(reciprocal, s.zero()),
        ]
    })
}

pub(super) fn massless_onshell_box<T: NativeFloat>(
    a: &[C<T>],
    c: &Constants<T>,
) -> Option<[C<T>; 3]> {
    if a[..4].iter().any(|z| !z.re.is_fully_zero())
        || a[6..10]
            .iter()
            .any(|z| !z.re.is_fully_zero() || !z.im.is_fully_zero())
        || a[4].re.is_fully_zero()
        || a[5].re.is_fully_zero()
    {
        return None;
    }
    let log_mu = a[10].re.log();
    let ls = momentum_log(&a[4].re, &log_mu, c);
    let lt = momentum_log(&a[5].re, &log_mu, c);
    let factor = (a[4].re.clone() * &a[5].re).inv();
    // Ls²+Lt²−(Ls−Lt)²−4π²/3 + (4π²/12) = 2 Ls Lt−π².
    let mut finite = scale(&mul(&ls, &lt), &c.two);
    finite.re -= &c.pi_squared;
    let pole = scale(
        &C::new(ls.re + lt.re, ls.im + lt.im),
        &(-c.two.clone() * &factor),
    );
    Some([
        scale(&finite, &factor),
        pole,
        C::new(c.four.clone() * factor, c.one.zero()),
    ])
}
