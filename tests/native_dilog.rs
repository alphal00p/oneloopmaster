//! Native arbitrary-precision complex dilogarithm against Symbolica's generic
//! `polylog(2, z)` at the same working precision (and, for double-double,
//! against the former 160-bit MPFR route).
use oneloop::{NativeFloat, double_dilog_via_float, symbolica_dilog_float};
use symbolica::domains::float::{Complex, DoubleFloat, Float, RealLike};

/// Arguments covering every branch of the argument reduction: inside and
/// outside the unit disk, `Re z` on both sides of 1/2, the unit circle, the
/// neighbourhood of 1, both cuts of the real axis and exact zero imaginary
/// parts (including `z = -1` and real `z > 1`).
fn arguments() -> Vec<(f64, f64)> {
    let mut points = vec![
        (0.0, 0.0),
        (0.5, 0.0),
        (1.0, 0.0),
        (-1.0, 0.0),
        (2.0, 0.0),
        (-3.5, 0.0),
        (1e-12, 0.0),
        (1.0 + 1e-9, 0.0),
        (1.0 - 1e-9, 1e-12),
        (0.5, 0.8660254037844386),
        (1.0, 1e-300),
        (0.0, 1.0),
        (-0.3, -0.7),
        (1e6, -3.0),
        (-2e5, 1e5),
    ];
    // Deterministic pseudo-random points in |z| < 4.
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..200 {
        points.push((8.0 * next() - 4.0, 8.0 * next() - 4.0));
    }
    points
}

fn relative(a: &Complex<Float>, b: &Complex<Float>) -> f64 {
    let difference = (a.re.clone() - &b.re).to_f64().hypot((a.im.clone() - &b.im).to_f64());
    let scale = b.re.to_f64().hypot(b.im.to_f64());
    if scale == 0.0 { difference } else { difference / scale }
}

#[test]
fn native_dilog_matches_symbolica_at_working_precision() {
    for bits in [106u32, 128, 256, 512] {
        let mut worst = 0.0f64;
        for (re, im) in arguments() {
            let z = Complex::new(Float::with_val(bits, re), Float::with_val(bits, im));
            let native = <Float as NativeFloat>::dilog(&z);
            let reference = symbolica_dilog_float(&z);
            let error = relative(&native, &reference);
            assert!(
                error <= 2f64.powi(-(bits as i32) + 10),
                "{bits} bits, z = ({re}, {im}): relative error {error:e}"
            );
            worst = worst.max(error);
        }
        println!("{bits} bits: worst relative difference {worst:e}");
    }
}

#[test]
fn native_double_double_dilog_matches_the_mpfr_route() {
    let mut worst = 0.0f64;
    for (re, im) in arguments() {
        let z = Complex::new(DoubleFloat::from(re), DoubleFloat::from(im));
        let native = <DoubleFloat as NativeFloat>::dilog(&z);
        let reference = double_dilog_via_float(&z);
        let lift = |w: &Complex<DoubleFloat>| {
            let pair = |x: DoubleFloat| {
                let inner = x.into_inner();
                Float::with_val(160, inner.hi()) + Float::with_val(160, inner.lo())
            };
            Complex::new(pair(w.re), pair(w.im))
        };
        let error = relative(&lift(&native), &lift(&reference));
        assert!(error <= 2f64.powi(-96), "z = ({re}, {im}): relative error {error:e}");
        worst = worst.max(error);
    }
    println!("double-double: worst relative difference {worst:e}");
}
