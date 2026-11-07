//! C0 at a vanishing Källén function (Gram determinant) through its exact
//! reduction to bubbles.
//!
//! With the modified Cayley matrix `Y_ij = m_i + m_j - p_ij` (`p_ii = 0`, `p_ij`
//! the invariant between propagators `i` and `j`), the Feynman-parameter
//! denominator is `x^T Y x / 2` on the simplex, and any `b` with `Y b = e`
//! (`e = (1, 1, 1)`) gives the exact identity
//!
//! `C0 = -Σ_i b_i B0^(i) + (N - d - 1) (Σ_i b_i) C0^(d+2)`
//!
//! (`B0^(i)` is the bubble without propagator `i`; OneLOop's normalisation of
//! both masters, so the bubbles' unit simple poles cancel against the last
//! term's). `Σ_i b_i = e^T Y^-1 e` is proportional to the Gram determinant,
//! `-λ(p1, p2, p3)/4`, so where the Källén function vanishes the triangle is
//! `C0 = -Σ_i b_i B0^(i)` exactly, in every mass sector, at every order of its
//! arguments, and the bubbles' poles cancel (`Σ_i b_i = 0`).
//!
//! The direct forms are not built for that point: the roots of their
//! quadratics coincide, and where a sector forms the discriminant `b² - 4ac`
//! from rounded coefficients (the one-mass sector, whose `b² - 4ac` is
//! `λ/(|m|² m)`) the roots split by rounding noise, every divided difference
//! is formed from that noise, and the arbitrary-precision result was wrong by
//! up to the whole value, at times with every bit tracked (numerica keeps the
//! tracked precision of a factor of exactly `±1` and of `log` next to one).
//!
//! Two cases are reduced, decided by exact predicates on the represented
//! input values ([NativeFloat::exact_rational]): the Källén function and
//! `det Y` in rational arithmetic, so that the decision does not depend on
//! the working precision (at 128 bits a rounded `det Y` of a zero-invariant
//! triangle is a rounding residue):
//!
//! - `det Y ≠ 0`: `b = adj(Y) e / det Y`, unique. An infrared-divergent
//!   triangle (a massless line between two on-shell legs, or a light-like leg
//!   between two massless lines) has `det Y = 0`, so this case is infrared
//!   finite and never takes a scaleless bubble.
//! - all invariants zero (`Y = m e^T + e m^T`, rank two): any `b` with
//!   `e^T b = 0` and `m^T b = 1` solves `Y b = e`; `b = (e_a - e_c)/(m_a - m_c)`
//!   for the pair of masses that differ most. Two massless lines (an infrared
//!   divergence) or three equal masses are left to the sector code.
//!
//! Otherwise (`det Y = 0` with a nonzero invariant: `Y` has a null vector
//! orthogonal to `e`, and no `b` with `Σ_i b_i = 0` exists) the sector code is
//! used.
use super::primitives::{C, NativeFloat, add, inverse, mul, truth};
use symbolica::domains::rational::Rational;

/// Exact complex rationals for the predicates.
#[derive(Clone)]
struct Exact {
    re: Rational,
    im: Rational,
}

impl Exact {
    fn of<T: NativeFloat>(z: &C<T>) -> Option<Self> {
        Some(Self {
            re: z.re.exact_rational()?,
            im: z.im.exact_rational()?,
        })
    }

    fn add(&self, other: &Self) -> Self {
        Self {
            re: &self.re + &other.re,
            im: &self.im + &other.im,
        }
    }

    fn sub(&self, other: &Self) -> Self {
        Self {
            re: &self.re - &other.re,
            im: &self.im - &other.im,
        }
    }

    fn mul(&self, other: &Self) -> Self {
        Self {
            re: &(&self.re * &other.re) - &(&self.im * &other.im),
            im: &(&self.re * &other.im) + &(&self.im * &other.re),
        }
    }

    fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }
}

/// `det Y` of `Y_ij = m_i + m_j - p_ij` (input order p12, p23, p13, m1, m2, m3).
fn exact_determinant(z: &[Exact]) -> Exact {
    let y = |i: usize, j: usize| {
        let invariant = match (i.min(j), i.max(j)) {
            (0, 1) => &z[0],
            (1, 2) => &z[1],
            _ => &z[2],
        };
        if i == j {
            z[3 + i].add(&z[3 + i])
        } else {
            z[3 + i].add(&z[3 + j]).sub(invariant)
        }
    };
    let [a, b, c] = [y(0, 0), y(1, 1), y(2, 2)];
    let [d, e, f] = [y(0, 1), y(1, 2), y(0, 2)];
    let a11 = b.mul(&c).sub(&e.mul(&e));
    let a12 = f.mul(&e).sub(&d.mul(&c));
    let a13 = d.mul(&e).sub(&b.mul(&f));
    a.mul(&a11).add(&d.mul(&a12)).add(&f.mul(&a13))
}

fn sub<T: NativeFloat>(a: &C<T>, b: &C<T>) -> C<T> {
    C::new(a.re.clone() - &b.re, a.im.clone() - &b.im)
}

fn neg<T: NativeFloat>(a: &C<T>) -> C<T> {
    C::new(-a.re.clone(), -a.im.clone())
}

fn modulus_squared<T: NativeFloat>(a: &C<T>) -> T {
    a.re.clone() * &a.re + a.im.clone() * &a.im
}

/// The C0 `[finite, pole, double pole]` at a vanishing Källén function of
/// the invariants, from the bubbles `b0([p², m_i², m_j², mu²])`, or `None`
/// where the reduction does not apply (module documentation).
pub(super) fn c0<T: NativeFloat>(
    input: &[C<T>],
    mut b0: impl FnMut(&[C<T>]) -> [C<T>; 3],
) -> Option<[C<T>; 3]> {
    // The invariants are real (validated); the masses are converted only at
    // a vanishing Källén function.
    let [p1, p2, p3] = [0, 1, 2].map(|i| input[i].re.exact_rational());
    let (p1, p2, p3) = (p1?, p2?, p3?);
    let sum = &(&p3 - &p1) - &p2;
    let lambda = &(&sum * &sum) - &(&(&p1 * &p2) * &Rational::from(4));
    if !lambda.is_zero() {
        return None;
    }
    let exact = input[..6].iter().map(Exact::of).collect::<Option<Vec<_>>>()?;
    let all_zero = exact[..3].iter().all(Exact::is_zero);
    let m = [&input[3], &input[4], &input[5]];
    // The invariant between propagators i and j (input order p12, p23, p13).
    let invariant = |i: usize, j: usize| match (i.min(j), i.max(j)) {
        (0, 1) => &input[0],
        (1, 2) => &input[1],
        _ => &input[2],
    };
    let coefficients: [Option<C<T>>; 3] = if !all_zero && !exact_determinant(&exact).is_zero() {
        let y = |i: usize, j: usize| {
            if i == j {
                add(m[i], m[i])
            } else {
                sub(&add(m[i], m[j]), invariant(i, j))
            }
        };
        let [a, b, c] = [y(0, 0), y(1, 1), y(2, 2)];
        let [d, e, f] = [y(0, 1), y(1, 2), y(0, 2)];
        // b = adj(Y) e / det Y with the cofactors of the symmetric
        // Y = [[a, d, f], [d, b, e], [f, e, c]].
        let a11 = sub(&mul(&b, &c), &mul(&e, &e));
        let a22 = sub(&mul(&a, &c), &mul(&f, &f));
        let a33 = sub(&mul(&a, &b), &mul(&d, &d));
        let a12 = sub(&mul(&f, &e), &mul(&d, &c));
        let a13 = sub(&mul(&d, &e), &mul(&b, &f));
        let a23 = sub(&mul(&d, &f), &mul(&a, &e));
        let determinant = add(&add(&mul(&a, &a11), &mul(&d, &a12)), &mul(&f, &a13));
        let inverse_determinant = inverse(&determinant);
        [
            add(&add(&a11, &a12), &a13),
            add(&add(&a12, &a22), &a23),
            add(&add(&a13, &a23), &a33),
        ]
        .map(|cofactor_sum| Some(mul(&cofactor_sum, &inverse_determinant)))
    } else if all_zero {
        let masses = &exact[3..6];
        if masses.iter().filter(|z| z.is_zero()).count() >= 2 {
            return None;
        }
        // The pair of masses that differ most (b = (e_a - e_c)/(m_a - m_c)).
        let (first, second) = [(0, 1), (1, 2), (0, 2)]
            .into_iter()
            .filter(|&(i, j)| !masses[i].sub(&masses[j]).is_zero())
            .map(|(i, j)| (i, j, modulus_squared(&sub(m[i], m[j]))))
            .reduce(|best, next| if next.2 > best.2 { next } else { best })
            .map(|(i, j, _)| (i, j))?;
        let weight = inverse(&sub(m[first], m[second]));
        let mut coefficients = [None, None, None];
        coefficients[second] = Some(neg(&weight));
        coefficients[first] = Some(weight);
        coefficients
    } else {
        return None;
    };
    let mu = &input[6];
    let zero = C::new(input[0].re.zero(), input[0].re.zero());
    let mut finite = zero.clone();
    let mut terms = 0;
    for (i, coefficient) in coefficients.iter().enumerate() {
        let Some(coefficient) = coefficient else {
            continue;
        };
        if !truth(coefficient) {
            continue;
        }
        terms += 1;
        // The bubble without propagator i: its two other propagators and
        // the invariant between them.
        let (j, k) = match i {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        };
        let arguments = [invariant(j, k).clone(), m[j].clone(), m[k].clone(), mu.clone()];
        let bubble = b0(&arguments);
        finite = sub(&finite, &mul(coefficient, &bubble[0]));
    }
    // Σ b_i = 0 with Y b = e needs two nonzero coefficients.
    (terms >= 2).then(|| [finite, zero.clone(), zero])
}
