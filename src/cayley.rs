//! The modified Cayley matrix of a triangle, `Y_ij = m_i + m_j - p_ij`
//! (`p_ii = 0`; `p = [p12, p23, p13]`, the input order of C0), over any ring of
//! values: exact rationals (the predicates of `native::gram`), floating values
//! and expressions (`triangle::with_gram_reduction`).
//!
//! On the simplex `Σ_i x_i = 1` the Feynman-parameter denominator of C0 is
//! `x^T Y x / 2`, and `e^T adj(Y) e` (`e = (1, 1, 1)`) is proportional to the
//! Källén function `λ(p12, p23, p13)`.

/// Ring operations on the matrix entries.
pub(crate) trait Ring: Clone {
    fn plus(&self, other: &Self) -> Self;
    fn minus(&self, other: &Self) -> Self;
    fn times(&self, other: &Self) -> Self;
}

impl Ring for symbolica::atom::Atom {
    fn plus(&self, other: &Self) -> Self {
        self + other
    }

    fn minus(&self, other: &Self) -> Self {
        self - other
    }

    fn times(&self, other: &Self) -> Self {
        self * other
    }
}

/// The index in `[p12, p23, p13]` of the invariant between the distinct
/// propagators `i` and `j`.
pub(crate) fn invariant_index(i: usize, j: usize) -> usize {
    match (i.min(j), i.max(j)) {
        (0, 1) => 0,
        (1, 2) => 1,
        _ => 2,
    }
}

/// The two propagators other than `i` (the bubble without propagator `i`).
pub(crate) fn others(i: usize) -> (usize, usize) {
    match i {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    }
}

/// `Y_ij`.
pub(crate) fn entry<Z: Ring>(p: [&Z; 3], m: [&Z; 3], i: usize, j: usize) -> Z {
    if i == j {
        m[i].plus(m[i])
    } else {
        m[i].plus(m[j]).minus(p[invariant_index(i, j)])
    }
}

/// The cofactor sums `(adj(Y) e)_i` and `det Y`.
pub(crate) fn cofactor_sums<Z: Ring>(p: [&Z; 3], m: [&Z; 3]) -> ([Z; 3], Z) {
    let y = |i: usize, j: usize| entry(p, m, i, j);
    // Y = [[a, d, f], [d, b, e], [f, e, c]].
    let [a, b, c] = [y(0, 0), y(1, 1), y(2, 2)];
    let [d, e, f] = [y(0, 1), y(1, 2), y(0, 2)];
    let a11 = b.times(&c).minus(&e.times(&e));
    let a22 = a.times(&c).minus(&f.times(&f));
    let a33 = a.times(&b).minus(&d.times(&d));
    let a12 = f.times(&e).minus(&d.times(&c));
    let a13 = d.times(&e).minus(&b.times(&f));
    let a23 = d.times(&f).minus(&a.times(&e));
    let determinant = a.times(&a11).plus(&d.times(&a12)).plus(&f.times(&a13));
    (
        [
            a11.plus(&a12).plus(&a13),
            a12.plus(&a22).plus(&a23),
            a13.plus(&a23).plus(&a33),
        ],
        determinant,
    )
}

/// `Y_ra - Y_rb = m_a - m_b - p_ra + p_rb`, formed from the inputs rather than
/// from the rounded entries.
fn row_difference<Z: Ring>(p: [&Z; 3], m: [&Z; 3], r: usize, a: usize, b: usize) -> Z {
    let mut difference = m[a].minus(m[b]);
    if a != r {
        difference = difference.minus(p[invariant_index(r, a)]);
    }
    if b != r {
        difference = difference.plus(p[invariant_index(r, b)]);
    }
    difference
}

/// `e × Y_r`: orthogonal to `e` and to row `r`. Where `Y` has a null vector `v`
/// orthogonal to `e`, every row is orthogonal to `v`, so `e × Y_r` is
/// proportional to `v` for every row not proportional to `e`.
pub(crate) fn cross_row<Z: Ring>(p: [&Z; 3], m: [&Z; 3], r: usize) -> [Z; 3] {
    [
        row_difference(p, m, r, 2, 1),
        row_difference(p, m, r, 0, 2),
        row_difference(p, m, r, 1, 0),
    ]
}

/// `Y v`.
pub(crate) fn apply<Z: Ring>(p: [&Z; 3], m: [&Z; 3], v: &[Z; 3]) -> [Z; 3] {
    core::array::from_fn(|i| {
        entry(p, m, i, 0)
            .times(&v[0])
            .plus(&entry(p, m, i, 1).times(&v[1]))
            .plus(&entry(p, m, i, 2).times(&v[2]))
    })
}
