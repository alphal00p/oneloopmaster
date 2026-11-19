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
//! `[h^{N−1}] g_n(c + h) Π_far (c − z_f + h)^{−M_f}`. Masses closer than a
//! third of the larger one form a cluster expanded about its largest member
//! `c`, the other members entering through `Π (1 − δ t)^{−M}` (offsets `δ`):
//! `Σ_k H_k G_{N−1+k}`. A massless propagator's node contributes nothing:
//! there `α_i ≤ n`, and `g_n` and its first `n` derivatives vanish at zero.
//! Inputs are scaled by a power of two so that the largest squared mass is of
//! order one.
//!
//! The finite part is returned (the integral has no poles in this domain), or
//! `None` where the caller keeps the general expression: two or more massless
//! propagators, `ρ > 1/16`, more than [`MAX_ORDER`] orders for the precision
//! (binary64 `ρ ≲ 1/24`, double-double `ρ ≲ 1/415`), or a cluster whose
//! expansion converges too slowly.

use super::primitives::{C, NativeFloat, add, inverse, log, mul};

/// Largest order `n` of the momentum series: binary64 takes it up to
/// `ρ ≈ 1/24` (beyond, the general expression loses at most about 5 bits),
/// double-double up to `ρ ≈ 1/415`.
const MAX_ORDER: usize = 12;
/// Largest number of terms of a cluster's expansion about its largest member.
const MAX_CLUSTER_TERMS: usize = 256;

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

/// Terms `K` of a tail `Σ_{j>K} C(j + P, P) r^j` bounded by
/// `C(K + 1 + P, P) r^{K+1} / (1 − r)^{P+1}` below `target / 8` (past the
/// binomial's peak), or `None` beyond `limit`. Binary64 bookkeeping suffices
/// for a length.
fn binomial_terms(r: f64, polynomial: usize, target: f64, limit: usize) -> Option<usize> {
    if r == 0.0 {
        return Some(0);
    }
    if !(r < 1.0) {
        return None;
    }
    let p = polynomial as f64;
    let goal = (target / 8.0).ln() + (p + 1.0) * (1.0 - r).ln();
    // ln C(k + 1 + P, P) at k = 0.
    let mut log_binomial = (p + 1.0).ln();
    let mut k = 0usize;
    loop {
        let past_peak = (k as f64 + 2.0 + p) * r < k as f64 + 2.0;
        if past_peak && log_binomial + (k as f64 + 1.0) * r.ln() < goal {
            return Some(k);
        }
        k += 1;
        if k > limit {
            return None;
        }
        log_binomial += ((k as f64 + 1.0 + p) / (k as f64 + 1.0)).ln();
    }
}

/// A distinct squared mass and its propagators.
struct Node<T: NativeFloat> {
    value: C<T>,
    members: Vec<usize>,
    zero: bool,
}

/// Nodes expanded together about `centre` (their largest member); `offsets`
/// are the other members `(node, z − centre)`, `far` the remaining nodes
/// `(node, centre − z)`, `terms` the length of the expansion.
struct Cluster<T: NativeFloat> {
    centre: usize,
    log: C<T>,
    offsets: Vec<(usize, C<T>)>,
    far: Vec<(usize, C<T>)>,
    terms: usize,
}

/// Largest multiplicity of `node` at order `n`: `α_i ≤ n` for one
/// propagator, `α_i + α_j ≤ 2n` for two or three.
fn largest_multiplicity<T: NativeFloat>(node: &Node<T>, n: usize) -> usize {
    node.members.len() + if node.members.len() == 1 { n } else { 2 * n }
}

/// One cluster's tables at order `n`: `divided[M]` is the series of
/// `g_n(c + h) / (d + h)^M` for the first far node (`divided[0]` without
/// one), `second[M][q]` the coefficients of `(d + h)^{−M}` for a second far
/// node, `offsets[o][M][k]` those of `(1 − δ t)^{−M}`, `terms` the length of
/// the cluster expansion.
struct Tables<T: NativeFloat> {
    divided: Vec<Vec<C<T>>>,
    second: Vec<Vec<C<T>>>,
    offsets: Vec<Vec<Vec<C<T>>>>,
    terms: usize,
}

/// Taylor coefficients `[h^k] g_n(c + h)`, `g_n(z) = z^p ln z`, `p = n + 1`,
/// for `k = 0..=order` (`log_c = ln c`, `reciprocal[k] = 1/k`).
fn g_series<T: NativeFloat>(
    c: &C<T>,
    log_c: &C<T>,
    p: usize,
    order: usize,
    reciprocal: &[T],
) -> Vec<C<T>> {
    let prototype = &c.re;
    let mut powers = vec![C::new(prototype.one(), prototype.zero())];
    for j in 1..=p {
        powers.push(mul(&powers[j - 1], c));
    }
    let mut harmonic = vec![prototype.zero()];
    for j in 1..=p {
        harmonic.push(harmonic[j - 1].clone() + &reciprocal[j]);
    }
    let mut series = vec![complex_zero(prototype); order + 1];
    let mut binomial = prototype.one();
    for k in 0..=order.min(p) {
        // binom(p, k) c^{p−k} (ln c + H_p − H_{p−k}).
        let offset = harmonic[p].clone() - &harmonic[p - k];
        let bracket = C::new(log_c.re.clone() + &offset, log_c.im.clone());
        series[k] = scale(&mul(&powers[p - k], &bracket), &binomial);
        binomial = binomial * prototype.from_i64((p - k) as i64) * &reciprocal[k + 1];
    }
    if order > p {
        // p! (−1)^{k−p−1} (k−p−1)! / k! c^{p−k}.
        let inverse_c = inverse(c);
        series[p + 1] = scale(&inverse_c, &reciprocal[p + 1]);
        for k in p + 1..order {
            let ratio = prototype.from_i64(-((k - p) as i64)) * &reciprocal[k + 1];
            series[k + 1] = scale(&mul(&series[k], &inverse_c), &ratio);
        }
    }
    series
}

/// `series / (d + h)` in place.
fn divide_by_linear<T: NativeFloat>(series: &mut [C<T>], inverse_d: &C<T>) {
    series[0] = mul(&series[0], inverse_d);
    for j in 1..series.len() {
        series[j] = mul(&sub(&series[j], &series[j - 1]), inverse_d);
    }
}

/// `table[M][q] = [h^q] (d + h)^{−M}` for `M ≤ powers`, `q ≤ order`.
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
    /// Tables at order `n`, each only as large as the multiplicities of its
    /// nodes require.
    fn tables(&self, nodes: &[Node<T>], n: usize, reciprocal: &[T]) -> Tables<T> {
        let largest = |node: usize| largest_multiplicity(&nodes[node], n);
        // The multiplicities of all nodes add up to 2n + 3.
        let size = (largest(self.centre)
            + self.offsets.iter().map(|(k, _)| largest(*k)).sum::<usize>())
        .min(2 * n + 3);
        let gamma = g_series(
            &nodes[self.centre].value,
            &self.log,
            n + 1,
            size - 1 + self.terms,
            reciprocal,
        );
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
        let offsets = self
            .offsets
            .iter()
            .map(|(k, delta)| offset_powers(delta, largest(*k), self.terms, reciprocal))
            .collect();
        Tables {
            divided,
            second,
            offsets,
            terms: self.terms,
        }
    }

    /// This cluster's part of the divided difference for node multiplicities `m`.
    fn contribution(&self, tables: &Tables<T>, m: &[usize], prototype: &T) -> C<T> {
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
                let h = &tables.offsets[0][m[node]];
                for k in 0..=tables.terms {
                    value = add(&value, &mul(&h[k], &g[size - 1 + k]));
                }
            }
            _ => {
                let ((a, _), (b, _)) = (&self.offsets[0], &self.offsets[1]);
                let size = centre + m[*a] + m[*b];
                let (ha, hb) = (&tables.offsets[0][m[*a]], &tables.offsets[1][m[*b]]);
                for i in 0..=tables.terms {
                    let mut inner = complex_zero(prototype);
                    for j in 0..=tables.terms - i {
                        inner = add(&inner, &mul(&hb[j], &g[size - 1 + i + j]));
                    }
                    value = add(&value, &mul(&ha[i], &inner));
                }
            }
        }
        value
    }
}

/// Finite part of `C0(p1², p2², p3²; m1², m2², m3²)` (`input` as for the
/// evaluator, `mu²` last and unused) for small external invariants, accurate
/// to about `tolerance` relative, or `None` where the caller should keep the
/// general expression (see the module documentation).
pub(crate) fn c0_small_momentum<T: NativeFloat>(input: &[C<T>], tolerance: &T) -> Option<C<T>> {
    let prototype = &input[0].re;
    let zero_real = prototype.zero();
    let is_zero = |z: &C<T>| z.re.is_zero() && z.im.is_zero();
    let masses = [&input[3], &input[4], &input[5]];
    if masses.iter().filter(|z| is_zero(z)).count() > 1 {
        return None;
    }
    if masses.iter().any(|z| !is_zero(z) && !(z.re > zero_real)) {
        return None;
    }
    let largest = masses
        .iter()
        .map(|z| modulus(z))
        .fold(zero_real.clone(), |a, b| if b > a { b } else { a });
    let factor = power_of_two_inverse(&largest)?;
    let a: Vec<C<T>> = masses.iter().map(|z| scale(z, &factor)).collect();
    // s01 = p1², s12 = p2², s02 = p3² (the evaluator's leg order).
    let s = [
        input[0].re.clone() * &factor,
        input[1].re.clone() * &factor,
        input[2].re.clone() * &factor,
    ];
    let pairs = [(0usize, 1usize), (1, 2), (0, 2)];
    let mut rho = zero_real.clone();
    for (k, &(i, j)) in pairs.iter().enumerate() {
        if !s[k].is_zero() {
            let size = if s[k] < zero_real {
                -s[k].clone()
            } else {
                s[k].clone()
            };
            rho = rho + size / (a[i].re.clone() + &a[j].re);
        }
    }
    if !(rho <= constant(prototype, 1, 16)) {
        return None;
    }
    let order = geometric_terms(&rho, tolerance, MAX_ORDER)?;

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
        let larger = if mx > my { mx } else { my };
        modulus(&sub(&x.value, &y.value)) < larger * &third
    };
    // Groups: the closest close pair, joined by the third node when it is
    // close to either; every other node alone.
    let count = nodes.len();
    let mut group: Vec<usize> = (0..count).collect();
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
    if let Some((i, j, _)) = best {
        group[j] = i;
        if let Some(k) = (0..count).find(|&k| k != i && k != j) {
            if close(&nodes[k], &nodes[i]) || close(&nodes[k], &nodes[j]) {
                group[k] = i;
            }
        }
    }
    let mut clusters: Vec<Cluster<T>> = Vec::new();
    for leader in 0..count {
        let members: Vec<usize> = (0..count).filter(|&k| group[k] == leader).collect();
        if members.is_empty() || nodes[members[0]].zero {
            continue;
        }
        let centre = *members
            .iter()
            .max_by(|&&x, &&y| {
                modulus(&nodes[x].value)
                    .partial_cmp(&modulus(&nodes[y].value))
                    .unwrap()
            })
            .unwrap();
        let c = &nodes[centre].value;
        let offsets: Vec<(usize, C<T>)> = members
            .iter()
            .filter(|&&k| k != centre)
            .map(|&k| (k, sub(&nodes[k].value, c)))
            .collect();
        let far: Vec<(usize, C<T>)> = (0..count)
            .filter(|k| !members.contains(k))
            .map(|k| (k, sub(c, &nodes[k].value)))
            .collect();
        let mut radius = modulus(c);
        for (_, d) in &far {
            let distance = modulus(d);
            if distance < radius {
                radius = distance;
            }
        }
        let mut spread = zero_real.clone();
        for (_, delta) in &offsets {
            let size = modulus(delta);
            if size > spread {
                spread = size;
            }
        }
        // Every order needs the full length: a cluster's part of the divided
        // difference can cancel against the far nodes' residues. The terms do
        // not decay like (spread/radius)^k: the coefficients of (1 - δt)^{-M}
        // and those of a nonzero far node's (d + h)^{-M} grow like
        // C(k + M - 1, k), with multiplicities M up to 2n + 2 at the top order.
        let terms = if offsets.is_empty() {
            0
        } else {
            let top = |members: &mut dyn Iterator<Item = usize>| {
                members.map(|k| largest_multiplicity(&nodes[k], order)).max().unwrap_or(1) - 1
            };
            let polynomial = top(&mut offsets.iter().map(|(k, _)| *k))
                + top(&mut far.iter().map(|(k, _)| *k).filter(|&k| !nodes[k].zero));
            binomial_terms(
                (spread / &radius).to_f64(),
                polynomial,
                tolerance.to_f64(),
                MAX_CLUSTER_TERMS,
            )?
        };
        clusters.push(Cluster {
            centre,
            log: log(c),
            offsets,
            far,
            terms,
        });
    }
    let longest = clusters
        .iter()
        .map(|cluster| cluster.terms)
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
    for n in 0..=order {
        let p = n + 1;
        let tables: Vec<Tables<T>> = clusters
            .iter()
            .map(|cluster| cluster.tables(&nodes, n, &reciprocal))
            .collect();
        let mut term = complex_zero(prototype);
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
                let mut value = complex_zero(prototype);
                for (cluster, table) in clusters.iter().zip(&tables) {
                    value = add(
                        &value,
                        &cluster.contribution(table, &multiplicity, prototype),
                    );
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
                term = add(&term, &scale(&value, &(weight * monomial)));
            }
        }
        total = add(&total, &term);
    }
    let result = scale(&total, &(-factor));
    (result.re.is_finite() && result.im.is_finite()).then_some(result)
}
