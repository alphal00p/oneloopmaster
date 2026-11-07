//! C0 at a vanishing Källén function of its invariants (a vanishing Gram
//! determinant), in all six orders of its arguments
//! (`tests/data/c0_gram_singular.txt`).
//!
//! There the roots of the sectors' quadratics coincide. The one-mass sector
//! forms its discriminant `b² - 4ac = λ/(|m|² m)` from rounded coefficients,
//! so the roots split by rounding noise and its divided differences are formed
//! from that noise: native `Float` gave -5.871 with 512 tracked bits for
//! `C0(0, p, p; m, 0, 0) = -4.928` (the first row), and wrong values at every
//! precision, mostly with few tracked bits. Zero-invariant triangles were
//! wrong by up to 1e154 in the orders with the massless line between the two
//! masses, and the massless triangle was NaN in two of its six orders. Native
//! `Float` now reduces these triangles exactly to bubbles (`native::gram`).
//!
//! References: the closed form `-[(1+a) ln(1+a) - a ln a]/(-p)`, `a = -p/m`,
//! of the one-mass rows and the mpmath Feynman-parameter integral (50 digits).
//! Native `Float` at 128 to 1024 bits must match them within the precision it
//! tracks (8 bits of slack; the references bound the check at 1e-48), track
//! at least `bits - 32` bits, agree across the orders and return vanishing
//! poles. Binary64 and double-double keep their kernels: on the one-mass rows
//! binary64 is right and double-double keeps only binary64 accuracy (its
//! generated code splits the roots by its own rounding); with
//! `set_gram_reduction` they reach their precision on every row.
//!
//! The interpreted expressions (`EvaluationBackend::Expression`, the
//! arbitrary-precision masters of NinjaRS's `SymbolicaMasters<Float>`) reduce
//! the same triangles to bubbles (`triangle::with_gram_reduction`): checked at
//! 128 and 256 bits.
#![cfg(feature = "generated-evaluators")]
use oneloop::{EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, DoubleFloat, Float, FloatLike, RealLike};

type C = Complex<f64>;

/// Working precisions of the native `Float` route.
const BITS: [u32; 5] = [128, 192, 256, 512, 1024];
/// Relative accuracy of the 50-digit references.
const REFERENCE: f64 = 1e-48;
/// Bits a `Float` value may lose to the reduction's own cancellations.
const SLACK: u32 = 32;

struct Row {
    args: Vec<C>,
    reference: (String, String),
    class: String,
}

fn rows() -> Vec<Row> {
    include_str!("data/c0_gram_singular.txt")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let tokens: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(tokens.len(), 13, "malformed row {line}");
            let v: Vec<f64> = tokens[..10].iter().map(|x| x.parse().unwrap()).collect();
            let mut args: Vec<C> = v[1..4].iter().map(|&x| C::new(x, 0.0)).collect();
            args.extend(v[4..10].chunks_exact(2).map(|z| C::new(z[0], z[1])));
            args.push(C::new(v[0], 0.0));
            Row {
                args,
                reference: (tokens[10].to_string(), tokens[11].to_string()),
                class: tokens[12].to_string(),
            }
        })
        .collect()
}

/// The triangle in its six propagator relabellings.
fn images(args: &[C]) -> Vec<Vec<C>> {
    let invariant = |i: usize, j: usize| match (i.min(j), i.max(j)) {
        (0, 1) => args[0],
        (1, 2) => args[1],
        _ => args[2],
    };
    [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
        [0, 2, 1],
        [2, 1, 0],
        [1, 0, 2],
    ]
    .map(|order: [usize; 3]| {
        let mut image: Vec<C> = (0..3)
            .map(|k| invariant(order[k], order[(k + 1) % 3]))
            .collect();
        image.extend(order.iter().map(|&k| args[3 + k]));
        image.push(args[6]);
        image
    })
    .to_vec()
}

fn on_stack(action: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(256 * 1024 * 1024)
        .spawn(action)
        .unwrap()
        .join()
        .unwrap();
}

fn report(label: &str, failures: Vec<String>) {
    assert!(
        failures.is_empty(),
        "{label}: {} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn reference(row: &Row, bits: u32) -> Complex<Float> {
    Complex::new(
        Float::parse(&row.reference.0, Some(bits)).unwrap(),
        Float::parse(&row.reference.1, Some(bits)).unwrap(),
    )
}

fn distance(value: &Complex<Float>, other: &Complex<Float>, scale: f64) -> f64 {
    let error = (value.re.clone() - &other.re)
        .to_f64()
        .hypot((value.im.clone() - &other.im).to_f64())
        / scale;
    if error.is_nan() { f64::INFINITY } else { error }
}

/// The precision of the dominant component.
fn tracked(z: &Complex<Float>) -> u32 {
    if z.re.to_f64().abs() >= z.im.to_f64().abs() {
        z.re.get_precision()
    } else {
        z.im.get_precision()
    }
}

/// The relative uncertainty that `bits` tracked bits claim, with 8 bits of
/// slack for numerica's estimate.
fn claimed(bits: u32) -> f64 {
    2f64.powi(8 - bits as i32)
}

#[test]
fn float_matches_the_references_in_every_order() {
    on_stack(|| {
        let mut failures = Vec::new();
        for bits in BITS {
            let mut native =
                NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::C0, bits).unwrap();
            for (line, row) in rows().iter().enumerate() {
                let reference = reference(row, bits + 64);
                let scale = reference.re.to_f64().hypot(reference.im.to_f64());
                let mut first: Option<Complex<Float>> = None;
                for (number, image) in images(&row.args).into_iter().enumerate() {
                    let input: Vec<_> = image
                        .iter()
                        .map(|z| {
                            Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im))
                        })
                        .collect();
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(Float::new(bits), Float::new(bits))
                    });
                    native.evaluate(&input, &mut output).unwrap();
                    let value = output[0].clone();
                    let precision = tracked(&value);
                    let error = distance(&value, &reference, scale);
                    let spread = first
                        .as_ref()
                        .map_or(0.0, |first| distance(&value, first, scale));
                    let poles = output[1..]
                        .iter()
                        .all(|pole| pole.re.is_fully_zero() && pole.im.is_fully_zero());
                    if !(error <= REFERENCE.max(claimed(precision)))
                        || !(spread <= claimed(precision))
                        || precision + SLACK < bits
                        || !poles
                    {
                        failures.push(format!(
                            "row {line} ({}) image {number} at {bits} bits: {} {} (error {error:e}, \
                             spread {spread:e}, {precision} tracked bits, poles vanish: {poles})",
                            row.class, value.re, value.im
                        ));
                    }
                    first.get_or_insert(value);
                }
            }
        }
        report("native Float", failures);
    });
}

/// Binary64 and double-double: their default kernels on the one-mass rows,
/// and with the reduction on every row.
#[test]
fn binary64_and_double_float() {
    on_stack(|| {
        let mut failures = Vec::new();
        for reduced in [false, true] {
            let mut binary64 = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
            let mut double = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::C0).unwrap();
            binary64.set_gram_reduction(reduced);
            double.set_gram_reduction(reduced);
            // Double-double keeps binary64 accuracy without the reduction.
            let bounds = if reduced {
                [1e-13, 1e-29]
            } else {
                [1e-13, 1e-14]
            };
            for (line, row) in rows().iter().enumerate() {
                if !reduced && !row.class.starts_with("one-mass") {
                    continue;
                }
                let reference = reference(row, 256);
                let scale = reference.re.to_f64().hypot(reference.im.to_f64());
                for (number, image) in images(&row.args).into_iter().enumerate() {
                    let mut output = [C::new(f64::NAN, f64::NAN); 3];
                    binary64.evaluate(&image, &mut output).unwrap();
                    let value = Complex::new(
                        Float::with_val(256, output[0].re),
                        Float::with_val(256, output[0].im),
                    );
                    let error = distance(&value, &reference, scale);
                    let input: Vec<_> = image
                        .iter()
                        .map(|z| Complex::new(DoubleFloat::from(z.re), DoubleFloat::from(z.im)))
                        .collect();
                    let zero = Complex::new(DoubleFloat::from(0.0), DoubleFloat::from(0.0));
                    let mut output = [zero; 3];
                    double.evaluate(&input, &mut output).unwrap();
                    let wide = |x: DoubleFloat| {
                        let pair = x.into_inner();
                        Float::with_val(256, pair.hi()) + Float::with_val(256, pair.lo())
                    };
                    let value = Complex::new(wide(output[0].re), wide(output[0].im));
                    let double_error = distance(&value, &reference, scale);
                    if !(error <= bounds[0]) || !(double_error <= bounds[1]) {
                        failures.push(format!(
                            "row {line} ({}) image {number}, reduction {reduced}: binary64 {error:e}, \
                             double-double {double_error:e}",
                            row.class
                        ));
                    }
                }
            }
        }
        report("binary64 and double-double", failures);
    });
}

/// Off the vanishing Källén function the `Float` route is the sector code's,
/// bit for bit: `p3 = p2 (1 + 2^-40)` next to the first row.
#[test]
fn float_is_unchanged_next_to_the_singular_point() {
    on_stack(|| {
        let bits = 256;
        let row = &rows()[0];
        let mut input: Vec<_> = row
            .args
            .iter()
            .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
            .collect();
        let factor = Float::with_val(bits, 1) + Float::with_val(bits, 2f64.powi(-40));
        input[2] = Complex::new(input[1].re.clone() * &factor, Float::with_val(bits, 0));
        let mut values = Vec::new();
        for reduced in [true, false] {
            let mut native =
                NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::C0, bits).unwrap();
            native.set_gram_reduction(reduced);
            let mut output = core::array::from_fn::<_, 3, _>(|_| {
                Complex::new(Float::new(bits), Float::new(bits))
            });
            native.evaluate(&input, &mut output).unwrap();
            values.push(output);
        }
        for (with, without) in values[0].iter().zip(&values[1]) {
            assert_eq!(with.re.to_rational(), without.re.to_rational());
            assert_eq!(with.im.to_rational(), without.im.to_rational());
            assert_eq!(with.re.get_precision(), without.re.get_precision());
        }
    });
}

/// The interpreted expressions at 128 and 256 bits (module documentation).
#[test]
fn expressions_match_the_references_in_every_order() {
    on_stack(|| {
        let mut failures = Vec::new();
        for bits in [128, 256] {
            let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
                ScalarIntegral::C0,
                bits,
                EvaluationBackend::Expression,
            )
            .unwrap();
            for (line, row) in rows().iter().enumerate() {
                let reference = reference(row, bits + 64);
                let scale = reference.re.to_f64().hypot(reference.im.to_f64());
                for (number, image) in images(&row.args).into_iter().enumerate() {
                    let input: Vec<_> = image
                        .iter()
                        .map(|z| {
                            Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im))
                        })
                        .collect();
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(Float::new(bits), Float::new(bits))
                    });
                    expression.evaluate(&input, &mut output).unwrap();
                    let value = output[0].clone();
                    let precision = tracked(&value);
                    let error = distance(&value, &reference, scale);
                    if !(error <= REFERENCE.max(claimed(precision))) || precision + SLACK < bits {
                        failures.push(format!(
                            "row {line} ({}) image {number} at {bits} bits: {} {} (error {error:e}, \
                             {precision} tracked bits)",
                            row.class, value.re, value.im
                        ));
                    }
                }
            }
        }
        report("interpreted expressions", failures);
    });
}
