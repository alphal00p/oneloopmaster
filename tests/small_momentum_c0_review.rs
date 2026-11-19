#![cfg(feature = "generated-evaluators")]
//! Review of PR #21 (small-momentum C0): reproduction of the defects, through the public API only.
//!
//! Run (prints the tables; the `regression_*` tests FAIL on 53d2507 by design):
//!
//! ```text
//! cargo test --release --features generated-evaluators --test small_momentum_c0_review -- --nocapture --test-threads=1
//! ```
//!
//! Optional: `REVIEW_POINTS=<n>` random points per class (default 4000), `REVIEW_THREADS=<k>` worker threads for
//! the random scan (default 1; more than one thread needs `SYMBOLICA_LICENSE` in the environment, as does
//! `--test-threads` above 1).
//!
//! Conventions (those of `NativeEvaluator` for `ScalarIntegral::C0`): input
//! `[p1², p2², p3², m1², m2², m3², μ²]`, all `Complex`; `p1²` is the invariant between propagators 1 and 2
//! (`s01` in `small_momentum_c0.rs`), `p2²` between 2 and 3 (`s12`), `p3²` between 1 and 3 (`s02`); squared
//! masses, complex as `m² − i m Γ` (imaginary part ≤ 0); invariants real; output `[finite, 1/ε, 1/ε²]`.
//!
//! Errors: `max_slot |value − reference| / |reference finite|`. The reference is the general expression in
//! `Float` at 1024 bits (`NativeEvaluator::<Float>`, which never takes the series). It loses at most a few
//! times `log10(m²/|p²|)` digits (about 2·log10 for (0, m, m)) out of 308 at these points; `reference_check`
//! confirms it agrees with the same evaluation at 768 bits to better than 1e-150 relative at every point of
//! the table (measured: 5.8e-190).
//!
//! "series" = the PR's default evaluator; "general" = the same evaluator after `set_small_momentum(false)`.
use oneloop::{NativeEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, DoubleFloat, Float, RealLike};

type Point = [(f64, f64); 7];

struct Evaluators {
    f64_series: NativeEvaluator<f64>,
    f64_general: NativeEvaluator<f64>,
    dd_series: NativeEvaluator<DoubleFloat>,
    dd_general: NativeEvaluator<DoubleFloat>,
    reference: NativeEvaluator<Float>,
}

impl Evaluators {
    fn new() -> Self {
        let mut f64_general = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
        f64_general.set_small_momentum(false);
        let mut dd_general = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::C0).unwrap();
        dd_general.set_small_momentum(false);
        Self {
            f64_series: NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap(),
            f64_general,
            dd_series: NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::C0).unwrap(),
            dd_general,
            reference: reference_evaluator(1024),
        }
    }

    /// [f64 series, f64 general, dd series, dd general] relative errors against 1024 bits.
    fn errors(&mut self, p: &Point) -> [f64; 4] {
        let r = evaluate_reference(&mut self.reference, p, 1024);
        [
            rel_f64(&eval_f64(&mut self.f64_series, p), &r),
            rel_f64(&eval_f64(&mut self.f64_general, p), &r),
            rel_dd(&eval_dd(&mut self.dd_series, p), &r),
            rel_dd(&eval_dd(&mut self.dd_general, p), &r),
        ]
    }
}

fn reference_evaluator(bits: u32) -> NativeEvaluator<Float> {
    let mut e = NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::C0, bits).unwrap();
    e.set_small_momentum(false);
    e
}

fn reference(p: &Point, bits: u32) -> Vec<Complex<Float>> {
    evaluate_reference(&mut reference_evaluator(bits), p, bits)
}

fn evaluate_reference(e: &mut NativeEvaluator<Float>, p: &Point, bits: u32) -> Vec<Complex<Float>> {
    let x: Vec<_> = p
        .iter()
        .map(|&(re, im)| Complex::new(Float::with_val(bits, re), Float::with_val(bits, im)))
        .collect();
    let mut out = vec![Complex::new(Float::new(bits), Float::new(bits)); 3];
    e.evaluate(&x, &mut out).unwrap();
    out
}

fn eval_f64(e: &mut NativeEvaluator<f64>, p: &Point) -> Vec<Complex<f64>> {
    let x: Vec<_> = p.iter().map(|&(re, im)| Complex::new(re, im)).collect();
    let mut out = vec![Complex::new(0.0, 0.0); 3];
    e.evaluate(&x, &mut out).unwrap();
    out
}

fn eval_dd(e: &mut NativeEvaluator<DoubleFloat>, p: &Point) -> Vec<Complex<DoubleFloat>> {
    let x: Vec<_> = p
        .iter()
        .map(|&(re, im)| Complex::new(DoubleFloat::from(re), DoubleFloat::from(im)))
        .collect();
    let mut out = vec![Complex::new(DoubleFloat::from(0.0), DoubleFloat::from(0.0)); 3];
    e.evaluate(&x, &mut out).unwrap();
    out
}

fn wide(x: &DoubleFloat) -> Float {
    let d = x.clone().into_inner();
    Float::with_val(1024, d.hi()) + Float::with_val(1024, d.lo())
}

fn rel_f64(v: &[Complex<f64>], r: &[Complex<Float>]) -> f64 {
    let scale = r[0].re.to_f64().hypot(r[0].im.to_f64());
    let mut w = 0.0f64;
    for s in 0..3 {
        let dr = Float::with_val(1024, v[s].re) - r[s].re.clone();
        let di = Float::with_val(1024, v[s].im) - r[s].im.clone();
        let e = dr.to_f64().hypot(di.to_f64());
        w = if e.is_nan() { f64::INFINITY } else { w.max(e) };
    }
    w / scale
}

fn rel_dd(v: &[Complex<DoubleFloat>], r: &[Complex<Float>]) -> f64 {
    let scale = r[0].re.to_f64().hypot(r[0].im.to_f64());
    let mut w = 0.0f64;
    for s in 0..3 {
        let dr = wide(&v[s].re) - r[s].re.clone();
        let di = wide(&v[s].im) - r[s].im.clone();
        let e = dr.to_f64().hypot(di.to_f64());
        w = if e.is_nan() { f64::INFINITY } else { w.max(e) };
    }
    w / scale
}

/// The PR's series parameter ρ = Σ |s_ij| / (Re m_i² + Re m_j²).
fn rho(p: &Point) -> f64 {
    [(0, 3, 4), (1, 4, 5), (2, 3, 5)]
        .iter()
        .map(|&(s, i, j)| if p[s].0 == 0.0 { 0.0 } else { p[s].0.abs() / (p[i].0 + p[j].0) })
        .sum()
}

/// Squared masses `m` and invariants `(0.4, −0.35, 0.25)·(Re m_i² + Re m_j²)·ρ` for (p1², p2², p3²); μ² = 7.3.
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

fn r(x: f64) -> (f64, f64) {
    (x, 0.0)
}

/// The cases of the review's table. The first five are the blocking defect (a cluster of two nearby masses
/// with a third, nonzero mass as far node); the last two are silent fallbacks to the general expression.
fn table_cases() -> Vec<(&'static str, Point)> {
    vec![
        ("A (1, 0.66655, 0.45), rho 1/30", scaled([r(1.0), r(0.66655), r(0.45)], 1.0 / 30.0)),
        ("B (1, 0.66655, 0.45), rho 1/500", scaled([r(1.0), r(0.66655), r(0.45)], 1.0 / 500.0)),
        ("C (1, 0.7, 0.45), rho 1/30", scaled([r(1.0), r(0.7), r(0.45)], 1.0 / 30.0)),
        (
            "D (675.04, 445.21, 322.54)",
            [
                r(0.0),
                r(-11.506190046397727),
                r(-14.539770758302355),
                r(675.0403603627077),
                r(445.21451896598046),
                r(322.5358104582457),
                r(7.3),
            ],
        ),
        ("E SM (mZ2, mW2, mH2), rho 1/30", scaled([r(8315.0), r(6467.2), r(15625.0)], 1.0 / 30.0)),
        (
            "F two tiny masses (overflow)",
            [
                r(0.0),
                r(-5.107859805229796e-24),
                r(0.0),
                r(0.037034422757829706),
                r(6.334868353316832e-9),
                r(5.15208355328473e-9),
                r(7.3),
            ],
        ),
        (
            "G wide widths (MAX_CLUSTER_TERMS)",
            [
                r(2.4365910182777774e-7),
                r(-2.232458914626797e-7),
                r(-3.042667720433323e-7),
                (4084365.994061994, -234743.79550801977),
                (2835925.136725336, -1027301.8922336465),
                r(2774783.2466520304),
                r(7.3),
            ],
        ),
    ]
}

fn print_errors(name: &str, p: &Point, e: [f64; 4]) {
    println!(
        "{name:34} rho {:9.3e} | {:9.2e} {:9.2e} | {:9.2e} {:9.2e}",
        rho(p),
        e[0],
        e[1],
        e[2],
        e[3]
    );
}

const HEADER: &str = "case                               rho           | f64 series f64 general | dd series  dd general";

#[test]
fn review_failing_cases() {
    let mut ev = Evaluators::new();
    println!("\n{HEADER}");
    for (name, p) in table_cases() {
        print_errors(name, &p, ev.errors(&p));
    }
    println!("\nfull inputs [p1², p2², p3², m1², m2², m3², μ²] as (re, im):");
    for (name, p) in table_cases() {
        println!("{name}: {p:?}");
    }
}

#[test]
fn reference_check() {
    let mut worst = 0.0f64;
    for (_, p) in table_cases() {
        let (a, b) = (reference(&p, 1024), reference(&p, 768));
        let scale = a[0].re.to_f64().hypot(a[0].im.to_f64());
        for s in 0..3 {
            let d = (a[s].re.clone() - b[s].re.clone()).to_f64().hypot((a[s].im.clone() - b[s].im.clone()).to_f64());
            worst = worst.max(d / scale);
        }
    }
    println!("\nreference 1024 vs 768 bits, worst relative difference over the table: {worst:.2e}");
    assert!(worst < 1e-150);
}

/// Masses (1, 1−d, 0.45) at ρ = 1/30 across d = 1/3: below, the cluster is {1, 1−d} with far node 0.45;
/// above, {1−d, 0.45} with far node 1.
#[test]
fn review_regime_boundary_scan() {
    let mut ev = Evaluators::new();
    println!("\nmasses (1, 1-d, 0.45), rho 1/30\n{:>8} | f64 series f64 general | dd series  dd general", "d");
    let (mut below, mut above) = ([0.0f64; 4], [0.0f64; 4]);
    let n = 400;
    for k in 0..=n {
        let d = 0.30 + 0.06 * k as f64 / n as f64;
        let p = scaled([r(1.0), r(1.0 - d), r(0.45)], 1.0 / 30.0);
        let e = ev.errors(&p);
        let side = if d < 1.0 / 3.0 { &mut below } else { &mut above };
        for j in 0..4 {
            side[j] = side[j].max(e[j]);
        }
        if k % 25 == 0 || (d - 1.0 / 3.0).abs() < 2e-4 {
            println!("{d:8.5} | {:9.2e} {:9.2e} | {:9.2e} {:9.2e}", e[0], e[1], e[2], e[3]);
        }
    }
    println!("worst d < 1/3 | {:9.2e} {:9.2e} | {:9.2e} {:9.2e}", below[0], below[1], below[2], below[3]);
    println!("worst d > 1/3 | {:9.2e} {:9.2e} | {:9.2e} {:9.2e}", above[0], above[1], above[2], above[3]);
}

/// xorshift64, written out so that the scan is reproducible without an RNG crate.
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

const CLASSES: [&str; 12] = [
    "(0,m,m)",
    "(m,m,m)",
    "separated",
    "near-equal",
    "pair ~1/3 apart",
    "chain",
    "hierarchical",
    "two light",
    "(0,m,M)",
    "widths <=0.05m",
    "widths <=0.6m",
    "uniform random",
];

/// Squared masses of a class, in a random order, times a random overall scale 10^[-6, 8).
fn masses(rng: &mut Rng, class: usize) -> [(f64, f64); 3] {
    let m = match class {
        0 => [r(0.0), r(1.0), r(1.0)],
        1 => [r(1.0), r(1.0), r(1.0)],
        2 => loop {
            let a = rng.log_range(1.5, 30.0);
            let b = rng.log_range(1.5, 30.0);
            if (a / b).max(b / a) > 1.5 {
                break [r(1.0), r(a), r(b)];
            }
        },
        3 => [r(1.0), r(1.0 - rng.log_range(1e-14, 0.5)), r(1.0 - rng.log_range(1e-14, 0.5))],
        4 => {
            let d = rng.range(0.30, 0.36);
            let c = match rng.pick(3) {
                0 => 0.0,
                1 => rng.log_range(0.05, 20.0),
                _ => 1.0 - d * rng.range(0.4, 1.6),
            };
            [r(1.0), r(1.0 - d), r(c.max(0.0))]
        }
        5 => {
            let d1 = rng.range(0.2, 0.333);
            let d2 = rng.range(0.2, 0.333);
            [r(1.0), r(1.0 - d1), r((1.0 - d1) * (1.0 - d2))]
        }
        6 => [r(rng.log_range(1e-10, 0.1)), r(1.0), r(rng.log_range(1.0, 10.0))],
        7 => [r(rng.log_range(1e-8, 0.1)), r(rng.log_range(1e-8, 0.1)), r(1.0)],
        8 => [r(0.0), r(1.0), r(rng.log_range(0.2, 5.0))],
        9 | 10 => {
            let width = if class == 9 { 0.05 } else { 0.6 };
            let sub = [1, 2, 3, 4, 8][rng.pick(5)];
            let mut base = masses(rng, sub);
            for z in base.iter_mut() {
                if z.0 > 0.0 && rng.uniform() < 0.8 {
                    // m² − i m Γ with Γ/m up to `width`
                    z.1 = -z.0 * rng.range(0.0, width);
                }
            }
            base
        }
        _ => {
            let mut v = [r(rng.uniform()), r(rng.uniform()), r(rng.uniform())];
            if rng.uniform() < 0.2 {
                v[rng.pick(3)] = r(0.0);
            }
            v
        }
    };
    let order = [[0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1], [2, 1, 0], [1, 0, 2]][rng.pick(6)];
    let scale = 10f64.powf(rng.range(-6.0, 8.0));
    order.map(|k| (m[k].0 * scale, m[k].1 * scale))
}

/// Invariants with Σ |s_ij|/(Re m_i² + Re m_j²) = target: one, two or three nonzero, random signs.
fn invariants(rng: &mut Rng, ms: &[(f64, f64); 3], target: f64) -> [f64; 3] {
    let pattern = rng.pick(7);
    let mut w = [rng.uniform(), rng.uniform(), rng.uniform()];
    match pattern {
        0 => w = [w[0], 0.0, 0.0],
        1 => w = [0.0, w[1], 0.0],
        2 => w = [0.0, 0.0, w[2]],
        3 => w[rng.pick(3)] = 0.0,
        _ => {}
    }
    for x in w.iter_mut() {
        if rng.uniform() < 0.5 {
            *x = -*x;
        }
    }
    let den = [ms[0].0 + ms[1].0, ms[1].0 + ms[2].0, ms[0].0 + ms[2].0];
    let raw: f64 = w.iter().map(|x| x.abs()).sum();
    if raw == 0.0 {
        return [0.0; 3];
    }
    let t = target / raw;
    [w[0] * den[0] * t, w[1] * den[1] * t, w[2] * den[2] * t]
}

fn env(name: &str, default: usize) -> usize {
    std::env::var(name).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

/// Random points inside the series domain (ρ ≤ 1/16), fixed seed. Counts of relative errors above 1e-10 and
/// 1e-6 (binary64) and above 1e-26 and 1e-20 (double-double), series vs general expression, and the worst
/// binary64 series point of each class.
#[test]
fn review_random_in_domain_scan() {
    let per_class = env("REVIEW_POINTS", 4000);
    let threads = env("REVIEW_THREADS", 1).max(1);
    let mut rng = Rng(0x5EED_C0C0_2026_1119);
    let mut points: Vec<(usize, Point)> = Vec::new();
    for class in 0..CLASSES.len() {
        let mut kept = 0;
        while kept < per_class {
            let ms = masses(&mut rng, class);
            let target = if rng.pick(10) == 0 { 0.0 } else { rng.log_range(1e-16, 1.0 / 16.0) };
            let s = invariants(&mut rng, &ms, target);
            let p: Point = [r(s[0]), r(s[1]), r(s[2]), ms[0], ms[1], ms[2], r(7.3)];
            if rho(&p) <= 1.0 / 16.0 {
                points.push((class, p));
                kept += 1;
            }
        }
    }
    let chunk = points.len().div_ceil(threads);
    let errors: Vec<[f64; 4]> = std::thread::scope(|scope| {
        let handles: Vec<_> = points
            .chunks(chunk)
            .map(|part| {
                scope.spawn(move || {
                    let mut ev = Evaluators::new();
                    part.iter().map(|(_, p)| ev.errors(p)).collect::<Vec<_>>()
                })
            })
            .collect();
        handles.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    println!(
        "\n{per_class} points per class, rho <= 1/16, seed 0x5EEDC0C020261119\n{:16} | f64 series >1e-10 >1e-6 | f64 general >1e-10 >1e-6 | dd series >1e-26 >1e-20 | dd general >1e-26 >1e-20",
        "class"
    );
    let mut worst: Vec<(f64, Point)> = vec![(0.0, [r(0.0); 7]); CLASSES.len()];
    let mut counts = vec![[0usize; 8]; CLASSES.len()];
    for ((class, p), e) in points.iter().zip(&errors) {
        let c = &mut counts[*class];
        let th = [1e-10, 1e-6, 1e-10, 1e-6, 1e-26, 1e-20, 1e-26, 1e-20];
        let idx = [0, 0, 1, 1, 2, 2, 3, 3];
        for k in 0..8 {
            if !(e[idx[k]] <= th[k]) {
                c[k] += 1;
            }
        }
        if !(e[0] <= worst[*class].0) {
            worst[*class] = (e[0], *p);
        }
    }
    for (k, name) in CLASSES.iter().enumerate() {
        let c = counts[k];
        println!(
            "{name:16} |              {:5} {:5} |               {:5} {:5} |            {:5} {:5} |             {:5} {:5}",
            c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]
        );
    }
    println!("\nworst binary64 series point per class (relative error, rho, full input):");
    for (k, name) in CLASSES.iter().enumerate() {
        println!("{name}: {:.2e} rho {:.3e} {:?}", worst[k].0, rho(&worst[k].1), worst[k].1);
    }
}

// ---------------------------------------------------------------------------------------------------------
// REGRESSION TESTS. They FAIL on 53d2507 by design and should pass once the cluster expansion is fixed.
// Threshold: never clearly worse than the general expression, which is about 1e-14 (binary64) and 1e-28
// (double-double) at these points.
// ---------------------------------------------------------------------------------------------------------

/// Fails on 53d2507 (cases A, C, D, E). Passes with the binomial-aware cluster length.
#[test]
fn regression_cluster_with_far_node_binary64() {
    let mut ev = Evaluators::new();
    let cases = table_cases();
    for (name, p) in &cases[..5] {
        if name.starts_with('B') {
            continue;
        }
        let e = ev.errors(p);
        assert!(e[0] <= 1e-12, "{name}: binary64 series {:e} (general {:e})", e[0], e[1]);
    }
}

/// Fails on 53d2507 (case B). Passes with the binomial-aware cluster length.
#[test]
fn regression_cluster_with_far_node_double_double() {
    let mut ev = Evaluators::new();
    let (name, p) = &table_cases()[1];
    let e = ev.errors(p);
    assert!(e[2] <= 1e-27, "{name}: double-double series {:e} (general {:e})", e[2], e[3]);
}

/// Fails on 53d2507 AND with the length fix alone: these need per-cluster scaling (F) and a better cluster
/// choice or a longer expansion (G).
#[test]
fn regression_in_domain_fallbacks() {
    let mut ev = Evaluators::new();
    let cases = table_cases();
    let (f, pf) = &cases[5];
    let (g, pg) = &cases[6];
    let ef = ev.errors(pf);
    let eg = ev.errors(pg);
    assert!(
        ef[2] <= 1e-27 && eg[0] <= 1e-12,
        "{f}: double-double {:e}; {g}: binary64 {:e}",
        ef[2],
        eg[0]
    );
}
