//! C0 at a vanishing Källén function (Gram determinant) through its exact
//! reductions to bubbles.
//!
//! With the modified Cayley matrix `Y_ij = m_i + m_j - p_ij` ([crate::cayley])
//! the Feynman-parameter denominator is `x^T Y x / 2` on the simplex, and
//! `e^T adj(Y) e` (`e = (1, 1, 1)`) is proportional to the Källén function.
//! `B0^(i)` is the bubble without propagator `i` (OneLOop's normalisation of
//! both masters).
//!
//! - `det Y ≠ 0`: any `b` with `Y b = e` gives the exact identity
//!   `C0 = -Σ_i b_i B0^(i) + (N - d - 1) (Σ_i b_i) C0^(d+2)`, and
//!   `b = adj(Y) e / det Y` has `Σ_i b_i = e^T adj(Y) e / det Y = 0`, so
//!   `C0 = -Σ_i b_i B0^(i)` (the bubbles' poles cancel). An infrared-divergent
//!   triangle (a massless line between two on-shell legs, or a light-like leg
//!   between two massless lines) has `det Y = 0`, so this case is infrared
//!   finite and never takes a scaleless bubble.
//! - three vanishing invariants (`Y = m e^T + e m^T`, `det Y = 0`): any `b`
//!   with `e^T b = 0` and `m^T b = 1` solves `Y b = e`;
//!   `b = (e_a - e_c)/(m_a - m_c)` for the pair of masses that differ most. Two
//!   massless lines (an infrared divergence) or three equal masses are left to
//!   the sector code.
//! - `det Y = 0` at nonzero invariants: `Y` has a null vector `v` with
//!   `Σ_i v_i = 0` (rank two: `adj(Y) = c v v^T` and `e^T adj(Y) e = c (e·v)²`
//!   vanishes; rank one: its null plane meets `e^⊥`), proportional to `e × Y_r`
//!   for every row `Y_r` not proportional to `e`. `x^T Y x` is then constant
//!   along `v` on the simplex; integrating along `v` leaves the triangle's
//!   edges with linear weights, and for every `k` with `v_k ≠ 0`
//!   `C0 = Σ_{i≠k} w_i ∂B0^(i)/∂m_k`, `w_i = -v_i/v_k`, `Σ_i w_i = 1`,
//!   with `B0^(i)` the bubble `B0(p_jk; m_j, m_k)` and, from integration by
//!   parts of its Feynman-parameter integral (finite parts; the pole does not
//!   depend on the masses),
//!   `∂B0/∂m_k = [β (B0 + ln(m_k/μ²) - 2) + 2 m_j (ln(m_k/μ²) - ln(m_j/μ²))] / λ(p_jk, m_j, m_k)`,
//!   `β = m_k - m_j - p_jk`. Two identical propagators (`m_i = m_j`,
//!   `p_ij = 0`, `v = e_i - e_j`) give `C0 = ∂B0(p_ik; m_k', m_i)/∂m_i`. `k` is
//!   taken where `|v_k|` is largest (the weights then lie in `[0, 1]`), else
//!   at the first index whose terms are regular: `m_k ≠ 0` and
//!   `λ(p_jk, m_j, m_k) ≠ 0` for every `w_i ≠ 0`. Where none is regular (an
//!   infrared divergence, or a bubble at its threshold or pseudo-threshold)
//!   the sector code runs.
//!
//! The direct forms are not built for these points: the roots of their
//! quadratics coincide, and where a sector forms the discriminant `b² - 4ac`
//! from rounded coefficients (the one-mass sector, whose `b² - 4ac` is
//! `λ/(|m|² m)`) the roots split by rounding noise, every divided difference
//! is formed from that noise, and the arbitrary-precision result was wrong by
//! up to the whole value, at times with every bit tracked (numerica keeps the
//! tracked precision of a factor of exactly `±1` and of `log` next to one).
//! Two identical propagators returned 0 or values with a few tracked bits.
//!
//! The cases are decided by exact predicates on the represented input values
//! ([NativeFloat::exact_rational]), so that the decision does not depend on
//! the working precision (at 128 bits a rounded `det Y` of a zero-invariant
//! triangle is a rounding residue); the weights `w_i` are exact rationals.
use super::primitives::{C, NativeFloat, add, inverse, log, mul, truth};
use crate::cayley::{self, Ring};
use symbolica::domains::float::{Complex, SingleFloat};
use symbolica::domains::rational::Rational;

/// Exact complex rationals for the predicates.
type Exact = Complex<Rational>;

impl Ring for Exact {
    fn plus(&self, other: &Self) -> Self {
        self.clone() + other
    }

    fn minus(&self, other: &Self) -> Self {
        self.clone() - other
    }

    fn times(&self, other: &Self) -> Self {
        self.clone() * other
    }
}

impl<T: NativeFloat> Ring for C<T> {
    fn plus(&self, other: &Self) -> Self {
        add(self, other)
    }

    fn minus(&self, other: &Self) -> Self {
        sub(self, other)
    }

    fn times(&self, other: &Self) -> Self {
        mul(self, other)
    }
}

fn exact<T: NativeFloat>(z: &C<T>) -> Option<Exact> {
    Some(Complex::new(z.re.exact_rational()?, z.im.exact_rational()?))
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

/// The Källén function `λ(p, a, b) = (a - b - p)² - 4 b p` of a bubble.
fn kallen<Z: Ring>(p: &Z, a: &Z, b: &Z) -> Z {
    let beta = a.minus(b).minus(p);
    let bp = b.times(p);
    let four_bp = bp.plus(&bp).plus(&bp.plus(&bp));
    beta.times(&beta).minus(&four_bp)
}

/// The cases of `∂B0(p; b, a)/∂a` in dimensional regularisation, from
/// `-∫_0^1 σ D(σ)^(-1-ε) dσ`, `D(σ) = (1-σ) b + σ a - σ(1-σ) p` (module
/// documentation).
#[derive(Clone, Copy)]
enum Derivative {
    /// `a ≠ 0`, `λ(p, a, b) ≠ 0`: by parts, from `B0`.
    Regular,
    /// `a ≠ 0`, `λ = 0`, `b = 0` (`p = a`): `D = σ² a`, a soft pole.
    OnShell,
    /// `a = b`, `p = 0`: `D = a`.
    Constant,
    /// `a ≠ 0`, `λ = 0`, `b p ≠ 0`: `D = p (σ - r)²` with `r` off `[0, 1]`.
    PseudoThreshold,
    /// `a = 0`, `p = 0`, `b ≠ 0`: `D = (1-σ) b`.
    MasslessAtZero,
    /// `a = 0`, `b = p ≠ 0`: `D = p (1-σ)²`.
    MasslessOnShell,
    /// `a = 0`, `p ≠ 0`, `b ≠ p`: `D = (1-σ)(b - σ p)`.
    Massless,
    /// `a = b = p = 0`.
    Scaleless,
}

/// The case of `∂B0(p; b, a)/∂a`, or `None` at a threshold (`D` with a double
/// root inside `[0, 1]`: the derivative diverges).
fn derivative_case(p: &Exact, a: &Exact, b: &Exact) -> Option<Derivative> {
    if !a.is_zero() {
        if !kallen(p, a, b).is_zero() {
            return Some(Derivative::Regular);
        }
        if b.is_zero() {
            return Some(Derivative::OnShell);
        }
        if p.is_zero() {
            return Some(Derivative::Constant);
        }
        // q = (1 - r)/(-r) for the double root r of D; real and negative iff
        // r lies in (0, 1).
        let q = a.plus(p).minus(b) / &a.minus(b).minus(p);
        if q.im.is_zero() && q.re < Rational::zero() {
            return None;
        }
        return Some(Derivative::PseudoThreshold);
    }
    Some(match (b.is_zero(), p.is_zero()) {
        (true, true) => Derivative::Scaleless,
        (false, true) => Derivative::MasslessAtZero,
        _ if b.minus(p).is_zero() => Derivative::MasslessOnShell,
        _ => Derivative::Massless,
    })
}

/// The simple-pole coefficient of `∂B0(p; b, a)/∂a` in its case, exactly.
fn residue(case: Derivative, p: &Exact, a: &Exact, b: &Exact) -> Exact {
    let one = Complex::new(Rational::one(), Rational::zero());
    match case {
        Derivative::OnShell => one / &a.plus(a),
        Derivative::MasslessAtZero => one / b,
        Derivative::MasslessOnShell => -(one / &p.plus(p)),
        Derivative::Massless => one / &b.minus(p),
        _ => one.zero(),
    }
}

/// The exact plan of the mass-derivative reduction (`det Y = 0` at nonzero
/// invariants; module documentation): the index `k`, the exact weights `w_i`
/// and derivative cases of the nonzero terms, and the exact pole
/// `Σ_i w_i Res ∂B0^(i)/∂m_k` (zero for an infrared-finite triangle).
#[allow(clippy::type_complexity)]
fn derivative_plan(
    p: [&Exact; 3],
    m: [&Exact; 3],
) -> Option<(usize, [Option<(Exact, Derivative)>; 3], Exact)> {
    let v = (0..3)
        .map(|r| cayley::cross_row(p, m, r))
        .find(|v| v.iter().any(|z| !z.is_zero()))?;
    if !cayley::apply(p, m, &v).iter().all(SingleFloat::is_zero) {
        return None;
    }
    let largest = (0..3)
        .rev()
        .max_by(|&i, &j| v[i].norm_squared().cmp(&v[j].norm_squared()))?;
    let plan = |k: usize| {
        if v[k].is_zero() {
            return None;
        }
        let mut terms = [None, None, None];
        let mut pole = v[k].zero();
        for i in (0..3).filter(|&i| i != k && !v[i].is_zero()) {
            // B0^(i) joins propagators j and k.
            let j = 3 - i - k;
            let q = p[cayley::invariant_index(j, k)];
            let case = derivative_case(q, m[k], m[j])?;
            let weight = -(v[i].clone() / &v[k]);
            pole = pole.plus(&weight.times(&residue(case, q, m[k], m[j])));
            terms[i] = Some((weight, case));
        }
        Some((k, terms, pole))
    };
    [largest, 0, 1, 2].into_iter().find_map(plan)
}

/// `ln(z/mu²)` on the lower lip of the negative real axis (`z - i0`).
fn physical_log<T: NativeFloat>(z: &C<T>, mu: &C<T>) -> C<T> {
    let mut ratio = mul(z, &inverse(mu));
    if ratio.im.is_fully_zero() {
        ratio.im = -ratio.im.zero();
    }
    log(&ratio)
}

/// The finite part of `∂B0(p; b, a)/∂a` in its exactly decided case
/// ([Derivative]), the regular one from the finite `B0(p; b, a)`.
fn mass_derivative<T: NativeFloat>(
    case: Derivative,
    [p, a, b, mu]: [&C<T>; 4],
    b0: impl FnOnce(&[C<T>]) -> [C<T>; 3],
) -> C<T> {
    let real = |n: i64| C::new(a.re.one().from_i64(n), a.re.zero());
    let over = |x: &C<T>, y: &C<T>| mul(x, &inverse(y));
    match case {
        Derivative::Regular => {
            let bubble = b0(&[p.clone(), b.clone(), a.clone(), mu.clone()])[0].clone();
            let log_a = physical_log(a, mu);
            let beta = sub(&sub(a, b), p);
            let mut numerator = mul(&beta, &sub(&add(&bubble, &log_a), &real(2)));
            if truth(b) {
                let difference = mul(b, &sub(&log_a, &physical_log(b, mu)));
                numerator = add(&numerator, &add(&difference, &difference));
            }
            over(&numerator, &kallen(p, a, b))
        }
        Derivative::OnShell => neg(&over(&physical_log(a, mu), &add(a, a))),
        Derivative::Constant => neg(&inverse(&add(a, a))),
        Derivative::PseudoThreshold => {
            let upper = sub(&add(p, a), b);
            let q = over(&upper, &sub(&sub(a, b), p));
            sub(&over(&real(2), &upper), &over(&log(&q), p))
        }
        Derivative::MasslessAtZero => over(&sub(&real(1), &physical_log(b, mu)), b),
        Derivative::MasslessOnShell => add(&over(&physical_log(p, mu), &add(p, p)), &inverse(p)),
        Derivative::Massless => {
            let shifted = sub(b, p);
            let log_shifted = physical_log(&shifted, mu);
            let mut finite = neg(&over(&log_shifted, &shifted));
            if truth(b) {
                let ratio = sub(&physical_log(b, mu), &log_shifted);
                finite = add(&finite, &over(&mul(b, &ratio), &mul(&shifted, p)));
            }
            finite
        }
        Derivative::Scaleless => real(0),
    }
}

/// The C0 `[finite, pole, double pole]` at a vanishing Källén function of
/// the invariants, from the bubbles `b0([p², m_i², m_j², mu²])`, or `None`
/// where no reduction applies (module documentation).
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
    let values = input[..6].iter().map(exact).collect::<Option<Vec<_>>>()?;
    let (exact_p, exact_m) = (
        [&values[0], &values[1], &values[2]],
        [&values[3], &values[4], &values[5]],
    );
    let all_zero = exact_p.iter().all(|z| z.is_zero());
    let p = [&input[0], &input[1], &input[2]];
    let m = [&input[3], &input[4], &input[5]];
    let invariant = |i: usize, j: usize| p[cayley::invariant_index(i, j)];
    let mu = &input[6];
    let zero = C::new(input[0].re.zero(), input[0].re.zero());
    let (_, determinant) = cayley::cofactor_sums(exact_p, exact_m);
    if !all_zero && determinant.is_zero() {
        let (k, terms, pole) = derivative_plan(exact_p, exact_m)?;
        let from_exact = |z: &Exact| {
            C::new(
                input[0].re.from_rational(&z.re),
                input[0].re.from_rational(&z.im),
            )
        };
        let mut finite = zero.clone();
        for (i, term) in terms.iter().enumerate() {
            let Some((weight, case)) = term else {
                continue;
            };
            let j = 3 - i - k;
            let value = mass_derivative(*case, [invariant(j, k), m[k], m[j], mu], &mut b0);
            finite = add(&finite, &mul(&from_exact(weight), &value));
        }
        return Some([finite, from_exact(&pole), zero]);
    }
    let coefficients: [Option<C<T>>; 3] = if !all_zero {
        // b = adj(Y) e / det Y.
        let (sums, determinant) = cayley::cofactor_sums(p, m);
        let inverse_determinant = inverse(&determinant);
        sums.map(|sum| Some(mul(&sum, &inverse_determinant)))
    } else {
        let masses = &values[3..6];
        if masses.iter().filter(|z| z.is_zero()).count() >= 2 {
            return None;
        }
        // The pair of masses that differ most (b = (e_a - e_c)/(m_a - m_c)).
        let (first, second) = [(0, 1), (1, 2), (0, 2)]
            .into_iter()
            .filter(|&(i, j)| !(masses[i].clone() - &masses[j]).is_zero())
            .map(|(i, j)| (i, j, modulus_squared(&sub(m[i], m[j]))))
            .reduce(|best, next| if next.2 > best.2 { next } else { best })
            .map(|(i, j, _)| (i, j))?;
        let weight = inverse(&sub(m[first], m[second]));
        let mut coefficients = [None, None, None];
        coefficients[second] = Some(neg(&weight));
        coefficients[first] = Some(weight);
        coefficients
    };
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
        let (j, k) = cayley::others(i);
        let arguments = [
            invariant(j, k).clone(),
            m[j].clone(),
            m[k].clone(),
            mu.clone(),
        ];
        let bubble = b0(&arguments);
        finite = sub(&finite, &mul(coefficient, &bubble[0]));
    }
    // Σ b_i = 0 with Y b = e needs two nonzero coefficients.
    (terms >= 2).then(|| [finite, zero.clone(), zero])
}
