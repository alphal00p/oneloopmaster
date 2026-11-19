//! Three-point function at small external momenta without cancellation.
//!
//! When every external invariant is small against the squared masses, the
//! general C0 expression cancels like `max(m²)/|p²|`: binary64 keeps about
//! `16 + log10(|p²|/m²)` digits (5.7e-5 relative at `|p²| = 1e-12 m²` for
//! masses `(0, m, m)`, 1.5e-3 for `(m, m, m)`), nearly equal masses lose
//! `m²/|m1² − m2²|` even at zero momenta, and with masses `(0, m, m)` at
//! `p² = 0` both binary64 and double-double return zero. With Feynman
//! parameters on the triangle (`x0 + x1 + x2 = 1`, area element `dA`),
//! `D(x) = Σ x_i m_i²` and `Q(x) = Σ_{i<j} x_i x_j s_ij`,
//!
//! `C0 = −∫ dA / (D − Q) = −Σ_{n≥0} T_n`,   `T_n = ∫ dA Qⁿ / Dⁿ⁺¹`.
//!
//! With at most one massless propagator and `Re m_i² > 0` otherwise,
//! `|x_i x_j s_ij| ≤ |s_ij| Re D / (Re m_i² + Re m_j²)`, so `|Q| ≤ ρ |D|` with
//! `ρ = Σ_{i<j} |s_ij| / (Re m_i² + Re m_j²)` and the series converges like
//! `ρⁿ` (there is no infrared or collinear singularity in this domain). By the
//! Hermite–Genocchi formula each monomial of `Qⁿ` integrates to a confluent
//! divided difference of `g_n(z) = z^{n+1} ln z` over the squared masses,
//!
//! `∫ dA x^α / Dⁿ⁺¹ = c_n α! g_n[m0²^(α0+1), m1²^(α1+1), m2²^(α2+1)]`,
//! `c_n = (−1)ⁿ / ((n+1)! n!)`,   `|α| = 2n`,
//!
//! since `g_n^{(2n+2)}(z) = z^{−n−1} / c_n`. Equal masses are one node of
//! summed multiplicity. A node `c` of multiplicity `N` contributes its residue
//! `[h^{N−1}] g_n(c + h) Π_far (c − z_f + h)^{−M_f}`. A massless propagator's
//! node contributes nothing: there `α_i ≤ n`, and `g_n` and its first `n`
//! derivatives vanish at zero.
//!
//! Masses closer than a third of the larger one form a cluster expanded about
//! its largest member `c`, the other members entering through
//! `Π (1 − δ t)^{−M}` (offsets `δ`): `Σ_k H_k G_{N−1+k}`, truncated after
//! `K` shells (two offsets by the recurrence of the product). Without a
//! nonzero far node the shells are bounded by `C(n + k, n) rᵏ` (`r` the
//! largest offset over `|c|`) relative to the divided difference, so order
//! `n` needs only `ρⁿ` of the tolerance. A nonzero far node `z_f` adds the
//! growth of `(c − z_f + h)^{−M_f}`: the shells grow like
//! `C(k + P, P) (spread/radius)ᵏ`, with `P` the multiplicities of the largest
//! offset and the far node less two (the bound of the review of #21), and the
//! cluster's part can cancel against the far residue, so every order takes
//! the full tolerance unless the far node is [`FAR_SEPARATION`] times heavier
//! or lighter than the centre. Each cluster is evaluated in units of its
//! centre (`h = c u`), so that clusters of small masses do not overflow.
//! Inputs are scaled by a power of two so that the largest squared mass is of
//! order one.
//!
//! Residues of masses at moderate distances, and a cluster with a nearby far
//! node, lose up to about 100 units of the precision to cancellation. All
//! nonzero masses expanded as one cluster about the largest do not, so that
//! grouping is planned first and evaluated if its work fits
//! [`WORK_BUDGET`]; otherwise the grouping above, if it fits.
//!
//! The finite part is returned (the integral has no poles in this domain), or
//! `None` where the caller keeps the general expression: two or more massless
//! propagators, an invariant with an imaginary part, `ρ > 1/16`, more than
//! [`MAX_ORDER`] orders for the precision (binary64 `ρ ≲ 1/24`, double-double
//! `ρ ≲ 1/415`), or no grouping within [`MAX_CLUSTER_TERMS`] shells and
//! [`WORK_BUDGET`]. The work grows with the number of orders and with the
//! spread of the masses, so the budget declines at larger `ρ`, where the
//! general expression loses fewer digits: in the tests, declined points keep
//! 3e-11 relative in binary64 (ρ above 4e-4) and 2e-23 in double-double (ρ
//! above 3e-8).

use super::primitives::{C, NativeFloat, add, inverse, log, mul};

/// Largest order `n` of the momentum series: binary64 takes it up to
/// `ρ ≈ 1/24`, double-double up to `ρ ≈ 1/415`. Just beyond, the general
/// expression loses about `log10(1/ρ)` digits for comparable masses, and more
/// where two light masses are small against the third (about 1e-9 relative
/// in binary64 for squared masses 1e-8 of the third).
const MAX_ORDER: usize = 12;
/// Largest number of shells of a cluster's expansion about its largest member.
const MAX_CLUSTER_TERMS: usize = 256;
/// Largest number of terms of an evaluation: expansion shells, residue
/// products and table entries, summed over the monomials of every order (see
/// `Cluster::cost` and `Cluster::table_cost`). About 3 to 5 ns each in
/// binary64 and 50 ns in double-double, so the series costs at most about
/// 10 to 20 times the general expression.
const WORK_BUDGET: usize = 8192;
/// A nonzero far node this many times heavier or lighter than a cluster's
/// centre does not cancel against the cluster's part.
const FAR_SEPARATION: f64 = 4.0;

fn modulus<T: NativeFloat>(z: &C<T>) -> T {
    NativeFloat::hypot(&z.re, &z.im)
}

fn sub<T: NativeFloat>(a: &C<T>, b: &C<T>) -> C<T> {
    C::new(a.re.clone() - &b.re, a.im.clone() - &b.im)
}

fn scale<T: NativeFloat>(a: &C<T>, factor: &T) -> C<T> {
    C::new(a.re.clone() * factor, a.im.clone() * factor)
}

fn constant<T: NativeFloat>(prototype: &T, numerator: i64, denominator: i64) -> T {
    prototype.from_i64(numerator) / prototype.from_i64(denominator)
}

fn complex_zero<T: NativeFloat>(prototype: &T) -> C<T> {
    C::new(prototype.zero(), prototype.zero())
}

/// `2^-e` with `2^e` within a factor of two of `x > 0` (exact scaling).
fn power_of_two_inverse<T: NativeFloat>(x: &T) -> Option<T> {
    let value = x.to_f64();
    if !(value.is_finite() && value > 0.0) {
        return None;
    }
    let exponent = value.log2().floor() as i64;
    let mut factor = x.one();
    let step = x.from_i64(1 << 30);
    let two = x.from_i64(2);
    let (mut remaining, up) = if exponent < 0 {
        (-exponent, true)
    } else {
        (exponent, false)
    };
    while remaining > 0 {
        let (by, count) = if remaining >= 30 {
            (&step, 30)
        } else {
            (&two, 1)
        };
        factor = if up { factor * by } else { factor / by };
        remaining -= count;
    }
    Some(factor)
}

/// Terms `K` of a geometric tail of ratio `r` with `r^(K+1) / (1 − r)`
/// below `target / 8`, or `None` beyond `limit`.
fn geometric_terms<T: NativeFloat>(r: &T, target: &T, limit: usize) -> Option<usize> {
    if r.is_zero() {
        return Some(0);
    }
    let one = r.one();
    if !(r.clone() < one) {
        return None;
    }
    let target = target.clone() * (one - r) / r.from_i64(8);
    let mut bound = r.clone();
    let mut terms = 0;
    while bound > target {
        terms += 1;
        if terms > limit {
            return None;
        }
        bound = bound * r;
    }
    Some(terms)
}

/// Shells `K` with the tail `Σ_{k>K} C(k + P, P) rᵏ`, bounded by
/// `C(K + 1 + P, P) r^{K+1} / (1 − r)^{P+1}`, below `exp(log_target)`, or
/// `None` beyond `limit` (binary64 bookkeeping suffices for a length; the
/// bound is the one of the review of #21).
fn binomial_terms(r: f64, polynomial: usize, log_target: f64, limit: usize) -> Option<usize> {
    if r == 0.0 {
        return Some(0);
    }
    if !(r < 1.0) {
        return None;
    }
    let p = polynomial as f64;
    let goal = log_target.exp();
    // The bound at K = 0, then by its ratios (P ≤ 4 MAX_ORDER + 4 keeps it
    // far inside the binary64 range).
    let mut bound = (p + 1.0) * r / (1.0 - r).powi(polynomial as i32 + 1);
    for k in 0..=limit {
        if bound < goal {
            return Some(k);
        }
        bound *= r * (k as f64 + 2.0 + p) / (k as f64 + 2.0);
    }
    None
}

/// `x` before `y` among nonzero squared masses: larger modulus, then larger
/// real part, then larger imaginary part, so that a cluster's centre does not
/// depend on the propagator order.
fn larger<T: NativeFloat>(x: &C<T>, y: &C<T>) -> bool {
    let (mx, my) = (modulus(x), modulus(y));
    if mx != my {
        return mx > my;
    }
    if x.re != y.re {
        return x.re > y.re;
    }
    x.im > y.im
}

/// A distinct squared mass and its propagators.
struct Node<T: NativeFloat> {
    value: C<T>,
    members: Vec<usize>,
    zero: bool,
}

/// Nodes expanded together about `centre` (their largest member), in units of
/// the centre `c`: `offsets` are the other members `(node, (z − c)/c)`, `far`
/// the remaining nodes `(node, (c − z)/c)`, `terms[n]` the number of shells
/// at order `n`.
struct Cluster<T: NativeFloat> {
    centre: usize,
    inverse: C<T>,
    log: C<T>,
    offsets: Vec<(usize, C<T>)>,
    far: Vec<(usize, C<T>)>,
    terms: Vec<usize>,
}

/// Largest multiplicity of `node` at order `n`: `α_i ≤ n` for one
/// propagator, `α_i + α_j ≤ 2n` for two or three.
fn largest_multiplicity<T: NativeFloat>(node: &Node<T>, n: usize) -> usize {
    node.members.len() + if node.members.len() == 1 { n } else { 2 * n }
}

/// One cluster's tables at order `n`: `divided[M]` is the series of
/// `g_n(c (1 + u)) / (c^{n+1} (d + u)^M)` for the first far node (`divided[0]`
/// without one), `second[M][q]` the coefficients of `(d + u)^{−M}` for a
/// second far node, `offsets[M][k]` those of `(1 − δ t)^{−M}` for a single
/// offset, `terms` the number of shells.
struct Tables<T: NativeFloat> {
    divided: Vec<Vec<C<T>>>,
    second: Vec<Vec<C<T>>>,
    offsets: Vec<Vec<C<T>>>,
    terms: usize,
}

/// Taylor coefficients `[u^k] (1 + u)^p (ln c + ln(1 + u))`, those of
/// `g_n(c (1 + u)) / c^p` with `p = n + 1`, for `k = 0..=order`
/// (`log_c = ln c`, `reciprocal[k] = 1/k`).
fn g_series<T: NativeFloat>(log_c: &C<T>, p: usize, order: usize, reciprocal: &[T]) -> Vec<C<T>> {
    let prototype = &log_c.re;
    let mut harmonic = vec![prototype.zero()];
    for j in 1..=p {
        harmonic.push(harmonic[j - 1].clone() + &reciprocal[j]);
    }
    let mut series = vec![complex_zero(prototype); order + 1];
    let mut binomial = prototype.one();
    for k in 0..=order.min(p) {
        // binom(p, k) (ln c + H_p − H_{p−k}).
        let offset = harmonic[p].clone() - &harmonic[p - k];
        series[k] = C::new(
            (log_c.re.clone() + &offset) * &binomial,
            log_c.im.clone() * &binomial,
        );
        binomial = binomial * prototype.from_i64((p - k) as i64) * &reciprocal[k + 1];
    }
    if order > p {
        // p! (−1)^{k−p−1} (k−p−1)! / k!.
        series[p + 1] = C::new(reciprocal[p + 1].clone(), prototype.zero());
        for k in p + 1..order {
            let ratio = prototype.from_i64(-((k - p) as i64)) * &reciprocal[k + 1];
            series[k + 1] = scale(&series[k], &ratio);
        }
    }
    series
}

/// `series / (d + u)` in place.
fn divide_by_linear<T: NativeFloat>(series: &mut [C<T>], inverse_d: &C<T>) {
    series[0] = mul(&series[0], inverse_d);
    for j in 1..series.len() {
        series[j] = mul(&sub(&series[j], &series[j - 1]), inverse_d);
    }
}

/// `table[M][q] = [u^q] (d + u)^{−M}` for `M ≤ powers`, `q ≤ order`.
fn inverse_powers<T: NativeFloat>(
    d: &C<T>,
    powers: usize,
    order: usize,
    reciprocal: &[T],
) -> Vec<Vec<C<T>>> {
    let prototype = &d.re;
    let inverse_d = inverse(d);
    let mut leading = C::new(prototype.one(), prototype.zero());
    let mut table = Vec::with_capacity(powers + 1);
    for m in 0..=powers {
        let mut row = vec![complex_zero(prototype); order + 1];
        row[0] = leading.clone();
        for q in 0..order {
            let ratio = prototype.from_i64(-((m + q) as i64)) * &reciprocal[q + 1];
            row[q + 1] = scale(&mul(&row[q], &inverse_d), &ratio);
        }
        table.push(row);
        leading = mul(&leading, &inverse_d);
    }
    table
}

/// `table[M][k] = [t^k] (1 − δ t)^{−M}` for `M ≤ powers`, `k ≤ terms`.
fn offset_powers<T: NativeFloat>(
    delta: &C<T>,
    powers: usize,
    terms: usize,
    reciprocal: &[T],
) -> Vec<Vec<C<T>>> {
    let prototype = &delta.re;
    (0..=powers)
        .map(|m| {
            let mut row = vec![complex_zero(prototype); terms + 1];
            row[0] = C::new(prototype.one(), prototype.zero());
            for k in 0..terms {
                let ratio = prototype.from_i64((m + k) as i64) * &reciprocal[k + 1];
                row[k + 1] = scale(&mul(&row[k], delta), &ratio);
            }
            row
        })
        .collect()
}

impl<T: NativeFloat> Cluster<T> {
    /// Terms of one monomial's contribution at order `n` (shells and residue
    /// products), for the work estimate.
    fn cost(&self, nodes: &[Node<T>], n: usize) -> usize {
        match (self.offsets.len(), self.far.len()) {
            (0, 2) => largest_multiplicity(&nodes[self.centre], n),
            (0, _) => 1,
            (1, _) => self.terms[n] + 1,
            _ => 3 * (self.terms[n] + 1),
        }
    }

    /// Terms of the tables at order `n`, for the work estimate.
    fn table_cost(&self, nodes: &[Node<T>], n: usize) -> usize {
        let largest = |node: usize| largest_multiplicity(&nodes[node], n);
        let terms = self.terms[n];
        let length = 2 * n + 3 + terms;
        let mut cost = length;
        if let Some((k, _)) = self.far.first() {
            cost += largest(*k) * length;
        }
        if let Some((k, _)) = self.far.get(1) {
            cost += largest(*k) * largest(self.centre);
        }
        if let [(k, _)] = self.offsets.as_slice() {
            cost += (largest(*k) + 1) * (terms + 1);
        }
        cost
    }

    /// Tables at order `n`, each only as large as the multiplicities of its
    /// nodes require.
    fn tables(&self, nodes: &[Node<T>], n: usize, reciprocal: &[T]) -> Tables<T> {
        let terms = self.terms[n];
        let largest = |node: usize| largest_multiplicity(&nodes[node], n);
        // The multiplicities of all nodes add up to 2n + 3.
        let size = (largest(self.centre)
            + self.offsets.iter().map(|(k, _)| largest(*k)).sum::<usize>())
        .min(2 * n + 3);
        let gamma = g_series(&self.log, n + 1, size - 1 + terms, reciprocal);
        let mut divided = vec![gamma];
        if let Some((k, d)) = self.far.first() {
            let inverse_d = inverse(d);
            for m in 1..=largest(*k) {
                let mut next = divided[m - 1].clone();
                divide_by_linear(&mut next, &inverse_d);
                divided.push(next);
            }
        }
        let second = match self.far.get(1) {
            Some((k, d)) => inverse_powers(d, largest(*k), largest(self.centre), reciprocal),
            None => Vec::new(),
        };
        let offsets = match self.offsets.as_slice() {
            [(k, delta)] => offset_powers(delta, largest(*k), terms, reciprocal),
            _ => Vec::new(),
        };
        Tables {
            divided,
            second,
            offsets,
            terms,
        }
    }

    /// This cluster's part of the divided difference for node multiplicities
    /// `m`, in units of `c^{−n−1}`.
    fn contribution(
        &self,
        tables: &Tables<T>,
        m: &[usize],
        reciprocal: &[T],
        prototype: &T,
    ) -> C<T> {
        let centre = m[self.centre];
        let g = match self.far.first() {
            Some((k, _)) => &tables.divided[m[*k]],
            None => &tables.divided[0],
        };
        let mut value = complex_zero(prototype);
        match (self.offsets.len(), self.far.len()) {
            (0, 2) => {
                let e = &tables.second[m[self.far[1].0]];
                for q in 0..centre {
                    value = add(&value, &mul(&g[centre - 1 - q], &e[q]));
                }
            }
            (0, _) => value = g[centre - 1].clone(),
            (1, _) => {
                let (node, _) = self.offsets[0];
                let size = centre + m[node];
                let h = &tables.offsets[m[node]];
                for k in 0..=tables.terms {
                    value = add(&value, &mul(&h[k], &g[size - 1 + k]));
                }
            }
            _ => {
                // F(t) = (1 − a t)^{−A} (1 − b t)^{−B} satisfies
                // (k + 1) F_{k+1} = (k (a + b) + A a + B b) F_k
                //                   − (k − 1 + A + B) a b F_{k−1}.
                let ((na, a), (nb, b)) = (&self.offsets[0], &self.offsets[1]);
                let (ma, mb) = (m[*na], m[*nb]);
                let size = centre + ma + mb;
                let sum = add(a, b);
                let product = mul(a, b);
                let mut alpha = add(
                    &scale(a, &prototype.from_i64(ma as i64)),
                    &scale(b, &prototype.from_i64(mb as i64)),
                );
                let mut beta = scale(&product, &prototype.from_i64((ma + mb) as i64 - 1));
                let mut previous = complex_zero(prototype);
                let mut current = C::new(prototype.one(), prototype.zero());
                value = g[size - 1].clone();
                for k in 0..tables.terms {
                    let next = scale(
                        &sub(&mul(&alpha, &current), &mul(&beta, &previous)),
                        &reciprocal[k + 1],
                    );
                    previous = current;
                    current = next;
                    value = add(&value, &mul(&current, &g[size + k]));
                    alpha = add(&alpha, &sum);
                    beta = add(&beta, &product);
                }
            }
        }
        value
    }
}

/// How a cluster's expansion length follows from the order.
struct Lengths {
    /// Largest offset over the radius of convergence.
    ratio: f64,
    /// Offset nodes and nonzero far nodes.
    offsets: Vec<usize>,
    nonzero_far: Vec<usize>,
    /// Every nonzero far node [`FAR_SEPARATION`] times heavier or lighter.
    separated: bool,
}

impl Lengths {
    /// Shells at order `n`, or `None` beyond [`MAX_CLUSTER_TERMS`]: `full`
    /// is the log of the tolerance, `relaxed` that of the share of order `n`
    /// where nothing cancels against the cluster's part.
    fn at<T: NativeFloat>(
        &self,
        nodes: &[Node<T>],
        n: usize,
        full: f64,
        relaxed: f64,
    ) -> Option<usize> {
        if self.nonzero_far.is_empty() {
            // Shells bounded by C(n + k, n) rᵏ relative to the divided
            // difference.
            return binomial_terms(self.ratio, n, relaxed, MAX_CLUSTER_TERMS);
        }
        let top = |list: &[usize]| {
            list.iter()
                .map(|&k| largest_multiplicity(&nodes[k], n))
                .max()
                .unwrap_or(1)
                - 1
        };
        let polynomial = top(&self.offsets) + top(&self.nonzero_far);
        let target = if self.separated { relaxed } else { full };
        binomial_terms(self.ratio, polynomial, target, MAX_CLUSTER_TERMS)
    }
}

/// Clusters for `groups` (node indices; massless nodes are never expanded),
/// with their shells per order, or `None` where an expansion would be longer
/// than [`MAX_CLUSTER_TERMS`] or the evaluation would exceed
/// [`WORK_BUDGET`]. `counts[n]` is the number of monomials at order `n`,
/// `log_rho = ln ρ`.
fn plan<T: NativeFloat>(
    nodes: &[Node<T>],
    groups: &[Vec<usize>],
    order: usize,
    log_rho: f64,
    log_tolerance: f64,
    counts: &[usize],
) -> Option<Vec<Cluster<T>>> {
    let prototype = &nodes[0].value.re;
    let one = C::new(prototype.one(), prototype.zero());
    let full = log_tolerance - 8f64.ln();
    let mut clusters = Vec::new();
    let mut lengths: Vec<Option<Lengths>> = Vec::new();
    for members in groups {
        if nodes[members[0]].zero {
            continue;
        }
        let mut centre = members[0];
        for &k in &members[1..] {
            if larger(&nodes[k].value, &nodes[centre].value) {
                centre = k;
            }
        }
        let c = &nodes[centre].value;
        let unit = inverse(c);
        let offsets: Vec<(usize, C<T>)> = members
            .iter()
            .filter(|&&k| k != centre)
            .map(|&k| (k, mul(&sub(&nodes[k].value, c), &unit)))
            .collect();
        let far: Vec<(usize, C<T>)> = (0..nodes.len())
            .filter(|k| !members.contains(k))
            .map(|k| (k, mul(&sub(c, &nodes[k].value), &unit)))
            .collect();
        lengths.push((!offsets.is_empty()).then(|| {
            let mut radius = prototype.one();
            for (_, d) in &far {
                let distance = modulus(d);
                if distance < radius {
                    radius = distance;
                }
            }
            let mut spread = prototype.zero();
            for (_, delta) in &offsets {
                let size = modulus(delta);
                if size > spread {
                    spread = size;
                }
            }
            let nonzero: Vec<&(usize, C<T>)> =
                far.iter().filter(|(k, _)| !nodes[*k].zero).collect();
            Lengths {
                ratio: (spread / &radius).to_f64(),
                offsets: offsets.iter().map(|(k, _)| *k).collect(),
                nonzero_far: nonzero.iter().map(|(k, _)| *k).collect(),
                // A nonzero far node `z_f` at least `FAR_SEPARATION` times
                // heavier or lighter than the centre has a residue negligible
                // against the cluster's part (`|z_f / c| = |1 − d|`).
                separated: nonzero.iter().all(|(_, d)| {
                    let ratio = modulus(&sub(&one, d)).to_f64();
                    ratio >= FAR_SEPARATION || ratio <= 1.0 / FAR_SEPARATION
                }),
            }
        }));
        clusters.push(Cluster {
            centre,
            inverse: unit,
            log: complex_zero(prototype),
            offsets,
            far,
            terms: vec![0; order + 1],
        });
    }
    let mut work = 0;
    for (n, count) in counts.iter().enumerate() {
        // Order n needs only ρⁿ of the tolerance where nothing cancels against
        // a cluster's part.
        let share = if n == 0 { 0.0 } else { n as f64 * log_rho };
        let relaxed = full - share - ((order + 1) as f64).ln();
        for (cluster, rule) in clusters.iter_mut().zip(&lengths) {
            if let Some(rule) = rule {
                cluster.terms[n] = rule.at(nodes, n, full, relaxed)?;
            }
        }
        for cluster in &clusters {
            work += count * cluster.cost(nodes, n) + cluster.table_cost(nodes, n);
        }
        if work > WORK_BUDGET {
            return None;
        }
    }
    Some(clusters)
}

/// Finite part of `C0(p1², p2², p3²; m1², m2², m3²)` for small external
/// invariants, accurate to about `tolerance` relative (pass about the working
/// precision's unit roundoff), or `None` where the caller should keep the
/// general expression (see the module documentation). `input` is laid out as
/// for [`NativeEvaluator`](super::NativeEvaluator) with
/// [`ScalarIntegral::C0`](crate::ScalarIntegral::C0):
/// `[p1², p2², p3², m1², m2², m3², μ²]` (`μ²` is not used), with `p1²` the
/// invariant between propagators 1 and 2, `p2²` between 2 and 3 and `p3²`
/// between 1 and 3, and complex squared masses `m² − i m Γ`. Only binary64
/// and double-double are supported (`None` otherwise).
pub fn c0_small_momentum<T: NativeFloat>(input: &[C<T>], tolerance: &T) -> Option<C<T>> {
    if input.len() < 6 || T::FIXED_BITS.is_none() {
        return None;
    }
    let prototype = &input[0].re;
    let zero_real = prototype.zero();
    let is_zero = |z: &C<T>| z.re.is_zero() && z.im.is_zero();
    // The domain first, from the unscaled inputs (ρ is scale-invariant).
    let masses = [&input[3], &input[4], &input[5]];
    let mut massless = 0;
    for z in masses {
        if is_zero(z) {
            massless += 1;
        } else if !(z.re > zero_real) {
            return None;
        }
    }
    if massless > 1 || input[..3].iter().any(|s| !s.im.is_zero()) {
        return None;
    }
    // s01 = p1², s12 = p2², s02 = p3² (the evaluator's leg order).
    let pairs = [(0usize, 1usize), (1, 2), (0, 2)];
    let mut rho = zero_real.clone();
    let mut nonzero = 0;
    for (k, &(i, j)) in pairs.iter().enumerate() {
        let s = &input[k].re;
        if !s.is_zero() {
            nonzero += 1;
            let size = if *s < zero_real {
                -s.clone()
            } else {
                s.clone()
            };
            rho = rho + size / (masses[i].re.clone() + &masses[j].re);
        }
    }
    if !(rho <= constant(prototype, 1, 16)) {
        return None;
    }
    let order = geometric_terms(&rho, tolerance, MAX_ORDER)?;
    // Monomials per order: C(n + z − 1, z − 1) for z nonzero invariants.
    let counts: Vec<usize> = (0..=order)
        .map(|n| match nonzero {
            0 => usize::from(n == 0),
            1 => 1,
            2 => n + 1,
            _ => (n + 1) * (n + 2) / 2,
        })
        .collect();

    let largest = masses
        .iter()
        .map(|z| modulus(z))
        .fold(zero_real.clone(), |a, b| if b > a { b } else { a });
    let factor = power_of_two_inverse(&largest)?;
    let a: Vec<C<T>> = masses.iter().map(|z| scale(z, &factor)).collect();
    let s = [
        input[0].re.clone() * &factor,
        input[1].re.clone() * &factor,
        input[2].re.clone() * &factor,
    ];

    // Distinct squared masses.
    let mut nodes: Vec<Node<T>> = Vec::new();
    for (i, z) in a.iter().enumerate() {
        match nodes
            .iter_mut()
            .find(|node| node.value.re == z.re && node.value.im == z.im)
        {
            Some(node) => node.members.push(i),
            None => nodes.push(Node {
                value: z.clone(),
                members: vec![i],
                zero: is_zero(z),
            }),
        }
    }
    let third = constant(prototype, 1, 3);
    let close = |x: &Node<T>, y: &Node<T>| {
        if x.zero || y.zero {
            return false;
        }
        let (mx, my) = (modulus(&x.value), modulus(&y.value));
        let bigger = if mx > my { mx } else { my };
        modulus(&sub(&x.value, &y.value)) < bigger * &third
    };
    // Groups: the closest close pair, joined by the third node when it is
    // close to either; every other node alone.
    let count = nodes.len();
    let mut best: Option<(usize, usize, T)> = None;
    for i in 0..count {
        for j in i + 1..count {
            if close(&nodes[i], &nodes[j]) {
                let distance = modulus(&sub(&nodes[i].value, &nodes[j].value));
                if best.as_ref().is_none_or(|(_, _, d)| distance < *d) {
                    best = Some((i, j, distance));
                }
            }
        }
    }
    let mut default: Vec<Vec<usize>> = Vec::new();
    match best {
        None => default.extend((0..count).map(|k| vec![k])),
        Some((i, j, _)) => match (0..count).find(|&k| k != i && k != j) {
            None => default.push(vec![i, j]),
            Some(k) if close(&nodes[k], &nodes[i]) || close(&nodes[k], &nodes[j]) => {
                default.push(vec![i, j, k]);
            }
            Some(k) => {
                default.push(vec![i, j]);
                default.push(vec![k]);
            }
        },
    }
    // All nonzero masses as one cluster first: it does not cancel against
    // residues of masses at moderate distances or against a nearby far node,
    // but needs more shells as the masses spread. The grouping above second.
    let nonzero_nodes: Vec<usize> = (0..count).filter(|&k| !nodes[k].zero).collect();
    let mut candidates: Vec<Vec<Vec<usize>>> = Vec::new();
    if nonzero_nodes.len() >= 2
        && !default
            .iter()
            .any(|group| group.len() == nonzero_nodes.len())
    {
        candidates.push(vec![nonzero_nodes]);
    }
    candidates.push(default);
    let log_rho = if rho.is_zero() {
        0.0
    } else {
        rho.to_f64().ln()
    };
    let log_tolerance = tolerance.to_f64().ln();
    let mut clusters = candidates
        .iter()
        .find_map(|groups| plan(&nodes, groups, order, log_rho, log_tolerance, &counts))?;
    for cluster in &mut clusters {
        cluster.log = log(&nodes[cluster.centre].value);
    }
    let longest = clusters
        .iter()
        .flat_map(|cluster| cluster.terms.iter().copied())
        .max()
        .unwrap_or(0);
    let mut reciprocal = vec![prototype.zero()];
    for k in 1..=2 * order + longest + 8 {
        reciprocal.push(prototype.one() / prototype.from_i64(k as i64));
    }

    // Factorials, their inverses and powers of the invariants.
    let mut factorial = vec![prototype.one()];
    let mut inverse_factorial = vec![prototype.one()];
    for k in 1..=2 * order + 2 {
        factorial.push(factorial[k - 1].clone() * prototype.from_i64(k as i64));
        inverse_factorial.push(inverse_factorial[k - 1].clone() * &reciprocal[k]);
    }
    let powers: Vec<Vec<T>> = s
        .iter()
        .map(|x| {
            let mut column = vec![prototype.one()];
            for k in 1..=order {
                column.push(column[k - 1].clone() * x);
            }
            column
        })
        .collect();

    let mut total = complex_zero(prototype);
    let mut multiplicity = vec![0usize; count];
    let mut sums = vec![complex_zero(prototype); clusters.len()];
    for n in 0..=order {
        let p = n + 1;
        let tables: Vec<Tables<T>> = clusters
            .iter()
            .map(|cluster| cluster.tables(&nodes, n, &reciprocal))
            .collect();
        for sum in sums.iter_mut() {
            *sum = complex_zero(prototype);
        }
        for k01 in 0..=n {
            for k12 in 0..=n - k01 {
                let k02 = n - k01 - k12;
                let monomial = powers[0][k01].clone() * &powers[1][k12] * &powers[2][k02];
                if monomial.is_zero() {
                    continue;
                }
                let alpha = [k01 + k02, k01 + k12, k12 + k02];
                // Multiplicity of each node: Σ (α_i + 1) over its propagators.
                for (slot, node) in multiplicity.iter_mut().zip(&nodes) {
                    *slot = node.members.iter().map(|&i| alpha[i] + 1).sum();
                }
                // (−1)ⁿ α0! α1! α2! / ((n+1)! k01! k12! k02!).
                let mut weight = factorial[alpha[0]].clone()
                    * &factorial[alpha[1]]
                    * &factorial[alpha[2]]
                    * &inverse_factorial[p]
                    * &inverse_factorial[k01]
                    * &inverse_factorial[k12]
                    * &inverse_factorial[k02];
                if n % 2 == 1 {
                    weight = -weight;
                }
                let coefficient = weight * monomial;
                for ((sum, cluster), table) in sums.iter_mut().zip(&clusters).zip(&tables) {
                    let value = cluster.contribution(table, &multiplicity, &reciprocal, prototype);
                    *sum = add(sum, &scale(&value, &coefficient));
                }
            }
        }
        // Back from units of the centre: times c^{−n−1}, one factor at a time
        // (the intermediates lie between the two magnitudes).
        for (sum, cluster) in sums.iter().zip(&clusters) {
            let mut value = sum.clone();
            for _ in 0..p {
                value = mul(&value, &cluster.inverse);
            }
            total = add(&total, &value);
        }
    }
    let result = scale(&total, &(-factor));
    (result.re.is_finite() && result.im.is_finite()).then_some(result)
}
