//! Ahead-of-time generic Rust implementations of the native master expressions.
//!
//! The generated functions contain ordinary arithmetic and structured Rust
//! branches. They neither build nor execute a Symbolica/SymJIT evaluator.
#[cfg(feature = "rust-codegen")]
pub mod batch;
#[cfg_attr(not(feature = "generated-evaluators"), allow(dead_code))]
pub mod primitives;
mod dilog;
pub mod double_double;
mod small_momentum;
mod small_momentum_c0;
pub use small_momentum::b0_momentum_shift;
pub use small_momentum_c0::c0_small_momentum;
mod tensor;
pub use tensor::{NativeTensorEvaluator, TensorIntegral};
#[cfg(feature = "generated-evaluators")]
mod simple;
#[cfg(feature = "generated-evaluators")]
mod handwritten;
#[cfg(feature = "generated-evaluators")]
mod gram;
#[cfg(feature = "generated-evaluators")]
mod pseudo_triangle;
// The emitter has its own stable formatting and mechanically named SSA locals.
#[rustfmt::skip]
#[allow(clippy::all)]
#[cfg(feature = "generated-evaluators")]
mod generated_a0;
#[rustfmt::skip]
#[allow(clippy::all)]
#[cfg(feature = "generated-evaluators")]
mod generated_b0;
#[rustfmt::skip]
#[allow(clippy::all)]
#[cfg(feature = "generated-evaluators")]
mod generated_db0;
#[rustfmt::skip]
#[allow(clippy::all)]
#[cfg(feature = "generated-evaluators")]
mod generated_c0;
#[rustfmt::skip]
#[allow(clippy::all)]
#[cfg(feature = "generated-evaluators")]
mod generated_d0;
#[rustfmt::skip]
#[allow(clippy::all)]
#[cfg(feature = "generated-evaluators")]
mod generated_constants;

use crate::ScalarIntegral;
pub use primitives::NativeFloat;
use primitives::*;
#[cfg(any(feature = "generated-evaluators", feature = "rust-codegen"))]
mod constants;
#[cfg(any(feature = "generated-evaluators", feature = "rust-codegen"))]
pub use constants::{ConstantContext, ConstantSpec};

#[cfg(not(feature = "generated-evaluators"))]
mod disabled;
#[cfg(not(feature = "generated-evaluators"))]
pub use disabled::NativeEvaluator;

#[derive(Clone)]
#[cfg(feature = "generated-evaluators")]
struct Context<T: NativeFloat> {
    constants: ConstantContext<T>,
    simple: simple::Constants<T>,
    hand: handwritten::Constants<T>,
}

#[cfg(feature = "generated-evaluators")]
impl<T: NativeFloat> Context<T> {
    fn new(bits: u32) -> Result<Self, String> {
        let prototype = T::zero_at(bits);
        let constants = ConstantContext::new(generated_constants::CONSTANTS, bits)?;
        Ok(Self {
            constants,
            simple: simple::Constants::new(&prototype),
            hand: handwritten::Constants::new(&prototype),
        })
    }

    #[inline(always)]
    fn constant(&self, index: usize) -> &C<T> {
        self.constants.constant(index)
    }
}

/// A reusable direct Rust scalar master implementation.
///
/// `T` may be `f64`, Symbolica's `DoubleFloat`, or arbitrary-precision `Float`.
/// Inputs use the usual squared quantities and have `mu_squared` last; outputs
/// are `[finite, simple pole, double pole]`. This backend uses no runtime
/// evaluator or JIT. Constants are prepared once at construction.
#[derive(Clone)]
#[cfg(feature = "generated-evaluators")]
pub struct NativeEvaluator<T: NativeFloat> {
    family: ScalarIntegral,
    poles: crate::PoleSelection,
    bits: u32,
    context: Context<T>,
    arguments: Vec<C<T>>,
    /// Whether the hand-written C0 and D0 sector kernels run.
    handwritten: [bool; 2],
    /// Whether C0's massless sector (OneLOop trif0) takes its hand-written
    /// kernel where `handwritten[0]` is off (the `Float` default).
    massless_c0: bool,
    /// Whether a C0 at a vanishing Källén function is reduced to bubbles
    /// (`gram.rs`; the `Float` default).
    gram_c0: bool,
    /// Whether the hand-written C0 also takes three massive propagators with a
    /// negative Källén function (opt-in, see `set_handwritten_negative_kallen_c0`).
    negative_kallen_c0: bool,
    /// Relative accuracy target of the small-momentum B0 series (2^-(bits+4)).
    tolerance: T,
    small_momentum: bool,
    /// Binary64 D0: the double-double evaluator of the points next to a
    /// two-particle threshold ([`d0_next_to_threshold`]), built on first use.
    threshold_d0: Option<Box<NativeEvaluator<symbolica::domains::float::DoubleFloat>>>,
}

#[cfg(feature = "generated-evaluators")]
impl<T: NativeFloat> NativeEvaluator<T> {
    /// Construct at the numeric type's default working precision.
    pub fn new(family: ScalarIntegral) -> Result<Self, String> {
        Self::with_binary_precision(family, T::DEFAULT_BITS)
    }

    /// Construct with an explicit bit count. Fixed-precision types require
    /// their native bit count; use `Float` for arbitrary working precision.
    pub fn with_binary_precision(family: ScalarIntegral, bits: u32) -> Result<Self, String> {
        if bits == 0 {
            return Err("binary precision must be positive".into());
        }
        if let Some(fixed) = T::FIXED_BITS
            && bits != fixed
        {
            return Err(format!(
                "this numeric type has fixed precision {fixed} bits"
            ));
        }
        Ok(Self {
            family,
            poles: crate::PoleSelection::BOTH,
            bits,
            context: Context::new(bits)?,
            arguments: (0..family.arity())
                .map(|_| C::new(T::zero_at(bits), T::zero_at(bits)))
                .collect(),
            // Binary64 and double-double use the hand-written sector kernels.
            // `Float` keeps the generated expression code: its exact-zero tests
            // and precision bookkeeping follow the expression's own grouping.
            // Two exceptions for `Float`, whose generated code is wrong there
            // with every bit tracked: C0 at a vanishing Källén function is
            // reduced to bubbles (`gram.rs`), and the massless triangle (trif0)
            // takes the hand-written transcription of the current expression
            // (the generated direct forms depend on the argument order next to
            // light-like legs). Binary64 and double-double keep their kernels.
            handwritten: [T::FIXED_BITS.is_some(); 2],
            massless_c0: T::FIXED_BITS.is_none(),
            gram_c0: T::FIXED_BITS.is_none(),
            negative_kallen_c0: false,
            tolerance: {
                let half = T::zero_at(bits).one() / T::zero_at(bits).from_i64(2);
                (0..bits + 4).fold(T::zero_at(bits).one(), |acc, _| acc * &half)
            },
            small_momentum: true,
            threshold_d0: None,
        })
    }

    /// Select the hand-written sector kernels (the binary64 and double-double
    /// default) or the generated expression code for every sector (also for
    /// `Float`'s massless triangle). For validation only.
    #[doc(hidden)]
    pub fn set_handwritten(&mut self, enabled: bool) {
        self.handwritten = [enabled; 2];
        self.massless_c0 = enabled;
    }

    /// Reduce C0 at a vanishing Källén function of its invariants to bubbles
    /// (the `Float` default) or evaluate it with the sector code (the binary64
    /// and double-double default). For validation only.
    #[doc(hidden)]
    pub fn set_gram_reduction(&mut self, enabled: bool) {
        self.gram_c0 = enabled;
    }

    /// Opt into the hand-written 't Hooft--Veltman C0 for three massive
    /// propagators with a negative Källén function (off by default; it applies
    /// wherever the hand-written C0 kernel runs). By default that sector uses
    /// the generated expression code. In binary64 the two differ by up to
    /// 2.0e-12 relative (80,000 points; against 128 bits the hand-written path
    /// is 2.8e-12 in the worst class, the generated code 1.2e-12), and the
    /// hand-written path is about 1.9x faster.
    pub fn set_handwritten_negative_kallen_c0(&mut self, enabled: bool) {
        self.negative_kallen_c0 = enabled;
    }

    /// Integral family represented by this implementation.
    pub fn family(&self) -> ScalarIntegral {
        self.family
    }

    /// Working precision in bits (not decimal digits).
    pub fn binary_precision(&self) -> u32 {
        self.bits
    }

    /// Use the small-momentum B0 and C0 representations (the default) or the
    /// general expressions everywhere. For validation only.
    #[doc(hidden)]
    pub fn set_small_momentum(&mut self, enabled: bool) {
        self.small_momentum = enabled;
    }

    /// Retain the chosen UV/IR poles, leaving the finite result unchanged.
    pub fn with_pole_selection(mut self, selection: crate::PoleSelection) -> Self {
        self.set_pole_selection(selection);
        self
    }

    /// Change pole selection for subsequent calls on this evaluator.
    pub fn set_pole_selection(&mut self, selection: crate::PoleSelection) {
        self.poles = selection;
    }

    pub fn pole_selection(&self) -> crate::PoleSelection {
        self.poles
    }

    /// Evaluate one point. Shape/domain errors leave output unchanged.
    pub fn evaluate(&mut self, input: &[C<T>], output: &mut [C<T>]) -> Result<(), String> {
        self.validate_shape(input.len(), output.len(), 1)?;
        self.validate_point(input)?;
        self.evaluate_row(input, output);
        Ok(())
    }

    /// Evaluate flat row-major points. Validate every row before writing any
    /// output; empty batches are valid. No hidden worker threads are spawned.
    pub fn evaluate_batch(
        &mut self,
        input: &[C<T>],
        output: &mut [C<T>],
        rows: usize,
    ) -> Result<(), String> {
        self.validate_shape(input.len(), output.len(), rows)?;
        for point in input.chunks_exact(self.family.arity()) {
            self.validate_point(point)?;
        }
        for (point, result) in input
            .chunks_exact(self.family.arity())
            .zip(output.chunks_exact_mut(3))
        {
            self.evaluate_row(point, result);
        }
        Ok(())
    }

    fn evaluate_row(&mut self, input: &[C<T>], output: &mut [C<T>]) {
        let threshold = self.threshold_values(input);
        let input = if T::FIXED_BITS.is_some() {
            input
        } else {
            for (working, value) in self.arguments.iter_mut().zip(input) {
                working.re = value.re.at_precision(self.bits);
                working.im = value.im.at_precision(self.bits);
            }
            &self.arguments
        };
        let mut values = match threshold {
            Some(values) => values,
            None => self.family_values(input),
        };
        // A finite C0/D0 point whose kernels return NaN or infinity: invariants
        // within 2^-30 of the largest argument are taken as exactly zero and
        // the point is evaluated again (e.g. binary64 D0 with two opposite
        // legs at p^2 ~ 1e-10 against m^2 = 3e4 and s ~ 1e6: both binary64
        // kernels give NaN; the zero-invariant value agrees with 256 bits to
        // 1e-15). The change is below 2^-30 of the scale.
        let invariants = match self.family {
            ScalarIntegral::C0 => 3,
            ScalarIntegral::D0 => 6,
            _ => 0,
        };
        if invariants > 0 && values.iter().any(|z| !(z.re.is_finite() && z.im.is_finite())) {
            let size = |z: &C<T>| z.re.clone() * &z.re + z.im.clone() * &z.im;
            let mut scale = size(&input[0]);
            for z in input.iter().skip(1) {
                let value = size(z);
                if value > scale {
                    scale = value;
                }
            }
            let limit = scale * &input[0].re.from_rational(&(1, 1u64 << 60).into());
            let mut snapped = input.to_vec();
            let mut changed = false;
            for z in snapped.iter_mut().take(invariants) {
                if truth(z) && size(z) <= limit {
                    *z = C::new(z.re.zero(), z.re.zero());
                    changed = true;
                }
            }
            if changed {
                let retry = self.family_values(&snapped);
                if retry.iter().all(|z| z.re.is_finite() && z.im.is_finite()) {
                    values = retry;
                }
            }
        }
        // All finite terms have been calculated before separating the poles.
        if self.poles != crate::PoleSelection::BOTH
            && values.iter().all(|z| z.re.is_finite() && z.im.is_finite())
        {
            self.poles
                .apply(crate::PoleCoefficient::Scalar(self.family), &mut values)
                .expect("scalar pole coefficients have no unsupported ranks");
        }
        for (out, value) in output.iter_mut().zip(values) {
            *out = value;
        }
    }

    /// A binary64 D0 point next to a two-particle threshold
    /// ([`d0_next_to_threshold`]): the double-double kernels' values, rounded
    /// once. `None` for every other point, and where double-double is not
    /// finite (the binary64 kernels and their repair then decide).
    fn threshold_values(&mut self, input: &[C<T>]) -> Option<[C<T>; 3]> {
        use core::any::Any;
        use symbolica::domains::float::DoubleFloat;
        if self.family != ScalarIntegral::D0 || T::FIXED_BITS != Some(53) {
            return None;
        }
        let point: &[C<T>; 11] = input.try_into().ok()?;
        let point = (point as &dyn Any).downcast_ref::<[C<f64>; 11]>()?;
        if !d0_next_to_threshold(point) {
            return None;
        }
        if self.threshold_d0.is_none() {
            let wide = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::D0).ok()?;
            self.threshold_d0 = Some(Box::new(wide));
        }
        let wide = self.threshold_d0.as_mut()?;
        // The same kernels as this evaluator's (see the setters).
        wide.handwritten = self.handwritten;
        wide.massless_c0 = self.massless_c0;
        wide.gram_c0 = self.gram_c0;
        wide.negative_kallen_c0 = self.negative_kallen_c0;
        wide.small_momentum = self.small_momentum;
        let rounded = wide_d0(wide, point)?;
        (&rounded as &dyn Any).downcast_ref::<[C<T>; 3]>().cloned()
    }

    /// The family's kernels for one point (no repair).
    fn family_values(&self, input: &[C<T>]) -> [C<T>; 3] {
        match self.family {
        ScalarIntegral::A0 if !truth(&input[0]) => {
            let mut values = simple::zeros(&input[0].re);
            values[1] = input[0].clone();
            values
        }
        ScalarIntegral::A0 => generated_a0::evaluate(input, &self.context),
        ScalarIntegral::B0 => self.b0_values(input),
        ScalarIntegral::DB0 => generated_db0::evaluate(input, &self.context),
        ScalarIntegral::C0 => self
            .small_momentum_c0(input)
            .or_else(|| {
                self.gram_c0
                    .then(|| gram::c0(input, |bubble| self.b0_values(bubble)))
                    .flatten()
            })
            .or_else(|| simple::massless_triangle(input, &self.context.simple))
            .or_else(|| {
                if self.handwritten[0] {
                    handwritten::c0(input, &self.context.hand, self.negative_kallen_c0)
                } else if self.massless_c0 {
                    handwritten::massless_c0(input, &self.context.hand)
                } else {
                    None
                }
            })
            .unwrap_or_else(|| generated_c0::evaluate(input, &self.context)),
        // The hand-written sectors leave a soft line's on-shell pair on
        // a channel to the generated code, which evaluates the labelling
        // with both on the legs (crate::labelling). A soft three-mass box
        // with a vanishing channel, which the generated code would exchange
        // into the finite formula, or with one of its massive pairs at its
        // pseudo-threshold, takes the infrared sector directly
        // (handwritten::soft_box).
        ScalarIntegral::D0 => handwritten::soft_box::d0_any_labelling(input, &self.context.hand)
            .or_else(|| pseudo_triangle::d0(input, self.bits))
            .or_else(|| simple::massless_onshell_box(input, &self.context.simple))
            .or_else(|| {
                self.handwritten[1]
                    .then(|| handwritten::d0(input, &self.context.hand))
                    .flatten()
            })
            .unwrap_or_else(|| {
                let relabelled: Option<[C<T>; 11]> = crate::labelling::soft_order(input)
                    .map(|order| core::array::from_fn(|k| input[order[k]].clone()));
                let input = relabelled.as_ref().map_or(input, |point| &point[..]);
                generated_d0::evaluate(input, &self.context)
            }),
    }
    }

    /// B0 `[finite, pole, double pole]` of `[p², m1², m2², mu²]`.
    fn b0_values(&self, input: &[C<T>]) -> [C<T>; 3] {
        if input[..3].iter().all(|z| !truth(z)) {
            return simple::zeros(&input[0].re);
        }
        self.small_momentum_b0(input)
            .unwrap_or_else(|| generated_b0::evaluate(input, &self.context))
    }

    fn validate_shape(&self, inputs: usize, outputs: usize, rows: usize) -> Result<(), String> {
        let expected_inputs = rows
            .checked_mul(self.family.arity())
            .ok_or("input size overflow")?;
        let expected_outputs = rows.checked_mul(3).ok_or("output size overflow")?;
        if inputs != expected_inputs || outputs != expected_outputs {
            return Err(format!(
                "expected {expected_inputs} inputs and {expected_outputs} outputs for {rows} rows; received {inputs} and {outputs}"
            ));
        }
        Ok(())
    }

    /// C0 for small external invariants (binary64 and double-double): the
    /// general expression cancels like `max(m²)/|p²|` there, so evaluate the
    /// cancellation-free momentum series of `small_momentum_c0` (finite part;
    /// the integral has no poles in its domain). `Float` keeps the general
    /// expression, whose working precision absorbs the loss.
    fn small_momentum_c0(&self, input: &[C<T>]) -> Option<[C<T>; 3]> {
        if !self.small_momentum || T::FIXED_BITS.is_none() {
            return None;
        }
        let finite = small_momentum_c0::c0_small_momentum(input, &self.tolerance)?;
        let mut values = simple::zeros(&input[0].re);
        values[0] = finite;
        Some(values)
    }

    /// B0 for `0 < |p²| <= max(m²)/16`: the general expression cancels like
    /// `max(m²)/p²` there, so evaluate it at `p² = 0` exactly and add the
    /// cancellation-free momentum shift of `small_momentum`.
    fn small_momentum_b0(&self, input: &[C<T>]) -> Option<[C<T>; 3]> {
        if !self.small_momentum || !truth(&input[0]) {
            return None;
        }
        let shift = small_momentum::b0_momentum_shift(&input[0], &input[1], &input[2], &self.tolerance)?;
        let zero = C::new(input[0].re.zero(), input[0].re.zero());
        let at_zero = [zero, input[1].clone(), input[2].clone(), input[3].clone()];
        let mut values = generated_b0::evaluate(&at_zero, &self.context);
        values[0] = add(&values[0], &shift);
        Some(values)
    }

    fn validate_point(&self, input: &[C<T>]) -> Result<(), String> {
        validate_point(self.family, input)
    }
}

/// A binary64 D0 point that is a soft three-mass box with a vanishing channel
/// or one of its massive pairs at its pseudo-threshold
/// (`handwritten::soft_box`), for the SymJIT entry points, whose prebuilt
/// asset predates both; `None` for every other point.
#[cfg(feature = "generated-evaluators")]
pub(crate) fn soft_box_d0(point: &[C<f64>]) -> Option<[C<f64>; 3]> {
    if !handwritten::soft_box::candidate(point) {
        return None;
    }
    handwritten::soft_box::d0_any_labelling(point, &handwritten::Constants::new(&0.0))
}

/// Whether a binary64 D0 point has an invariant next to the two-particle
/// threshold `p = (√m_i + √m_j)²` of the two massive propagators it joins: with
/// `a = m_i + m_j - p` (squared masses), `Re a < 0` and the Källén function
/// `λ = a² - 4 m_i m_j` below `2^-12 |a|²` (`p` within about 6e-5 of the
/// threshold, relatively). The pairs are the legs `p1..p4` and the channels
/// `s12, s23` with the propagators they join.
///
/// The sector formulas take the masses' square roots (`k = a / (√m_i √m_j)`
/// and `√(k² - 4)`, or the contour's edge discriminants), so next to a
/// threshold the binary64 kernels lose about `2^-53 |a| / √|λ|`, and the
/// three-mass sectors about `2^-53 a² / |λ|`; at a threshold of masses
/// without exact square roots `√(k² - 4)` is `√(2^-53)` instead of zero. Two
/// onia at threshold (two quark pairs of mass² 2.78515625 at P/2 each, exact
/// in binary64) had boxes 5.3e-9 and 4.0e-9 off at their exact thresholds,
/// 5e-10 one part in 1e15 off them, 2e-11 at 1e-12. Double-double is exact
/// there (to 3e-32). Pseudo-thresholds (`Re a > 0`, among them equal masses
/// at `p = 0`) are not singular on the physical sheet, and binary64 keeps them
/// (within 2.4e-13 in the same scan).
#[cfg(feature = "generated-evaluators")]
fn d0_next_to_threshold(point: &[C<f64>; 11]) -> bool {
    // (invariant, propagator i, propagator j)
    const PAIRS: [(usize, usize, usize); 6] =
        [(0, 0, 1), (1, 1, 2), (2, 2, 3), (3, 3, 0), (4, 0, 2), (5, 1, 3)];
    PAIRS.iter().any(|&(k, i, j)| {
        let (mi, mj, p) = (point[6 + i], point[6 + j], point[k]);
        let zero = |z: C<f64>| z.re == 0.0 && z.im == 0.0;
        if zero(mi) || zero(mj) {
            return false;
        }
        let (are, aim) = (mi.re + mj.re - p.re, mi.im + mj.im - p.im);
        if are >= 0.0 {
            return false;
        }
        let (mre, mim) = (mi.re * mj.re - mi.im * mj.im, mi.re * mj.im + mi.im * mj.re);
        let (lre, lim) = (are * are - aim * aim - 4.0 * mre, 2.0 * are * aim - 4.0 * mim);
        let a2 = are * are + aim * aim;
        // |λ|² < 2^-24 |a|⁴
        lre * lre + lim * lim < a2 * a2 * (1.0 / 16_777_216.0)
    })
}

/// The double-double D0 kernels' values at a binary64 point, rounded once;
/// `None` where they are not finite.
#[cfg(feature = "generated-evaluators")]
fn wide_d0(
    wide: &NativeEvaluator<symbolica::domains::float::DoubleFloat>,
    point: &[C<f64>; 11],
) -> Option<[C<f64>; 3]> {
    use symbolica::domains::float::{DoubleFloat, RealLike};
    let lifted: Vec<C<DoubleFloat>> =
        point.iter().map(|z| C::new(DoubleFloat::from(z.re), DoubleFloat::from(z.im))).collect();
    let rounded = wide.family_values(&lifted).map(|z| C::new(z.re.to_f64(), z.im.to_f64()));
    rounded.iter().all(|z| z.re.is_finite() && z.im.is_finite()).then_some(rounded)
}

/// A binary64 D0 point next to a two-particle threshold
/// ([`d0_next_to_threshold`]), for the SymJIT and expression entry points: the
/// native double-double kernels' values, rounded once. `None` for every other
/// point.
#[cfg(feature = "generated-evaluators")]
pub(crate) fn threshold_d0(point: &[C<f64>]) -> Option<[C<f64>; 3]> {
    use std::cell::RefCell;
    use symbolica::domains::float::DoubleFloat;
    thread_local! {
        static WIDE: RefCell<Option<NativeEvaluator<DoubleFloat>>> = const { RefCell::new(None) };
    }
    let point: &[C<f64>; 11] = point.try_into().ok()?;
    if !d0_next_to_threshold(point) {
        return None;
    }
    WIDE.with(|cell| {
        let mut wide = cell.borrow_mut();
        if wide.is_none() {
            *wide = Some(NativeEvaluator::new(ScalarIntegral::D0).ok()?);
        }
        wide_d0(wide.as_ref()?, point)
    })
}

/// A binary64 D0 point with three massive propagators pairwise at their
/// pseudo-thresholds (`pseudo_triangle`), for the SymJIT and expression entry
/// points: the native closed form's value. `None` for every other point.
#[cfg(feature = "generated-evaluators")]
pub(crate) fn pseudo_triangle_d0(point: &[C<f64>]) -> Option<[C<f64>; 3]> {
    pseudo_triangle::d0(point, 53)
}

/// [pseudo_triangle_d0] for the interpreted `Float` at `bits`.
#[cfg(feature = "generated-evaluators")]
pub(crate) fn pseudo_triangle_d0_float(
    point: &[C<symbolica::domains::float::Float>],
    bits: u32,
) -> Option<[C<symbolica::domains::float::Float>; 3]> {
    pseudo_triangle::d0(point, bits)
}

/// Whether a binary64 D0 point has four real masses and a massive pair at its
/// pseudo-threshold: an invariant `p` between lines `i` and `j` with
/// `(m_i - m_j)² - 2 p (m_i + m_j) + p² = 0` exactly and `m_i + m_j - p` not
/// negative (squared masses), the test of `box_integral::pair_thresholds`.
#[cfg(feature = "generated-evaluators")]
fn d0_four_masses_at_pseudo_threshold(point: &[C<f64>; 11]) -> bool {
    // (invariant, propagator i, propagator j), as in d0_next_to_threshold.
    const PAIRS: [(usize, usize, usize); 6] = [
        (0, 0, 1),
        (1, 1, 2),
        (2, 2, 3),
        (3, 3, 0),
        (4, 0, 2),
        (5, 1, 3),
    ];
    if point[6..10].iter().any(|m| m.im != 0.0 || m.re == 0.0) {
        return false;
    }
    PAIRS.iter().any(|&(k, i, j)| {
        let (mi, mj, p) = (point[6 + i].re, point[6 + j].re, point[k]);
        let (difference, sum) = (mi - mj, mi + mj);
        let kallen = difference * difference - 2.0 * p.re * sum + p.re * p.re;
        p.im == 0.0 && kallen == 0.0 && sum - p.re >= 0.0
    })
}

/// A binary64 D0 point with four real masses and a massive pair at its
/// pseudo-threshold ([`d0_four_masses_at_pseudo_threshold`]), for the SymJIT
/// entry points: the native kernels' values. The prebuilt asset predates
/// `box_integral::chart_pair`: the rounding of the pair's `k` decided whether
/// it takes the four-mass chart's `p23` or `p12` slot, and some labellings
/// came out wrong. `None` for every other point.
#[cfg(feature = "generated-evaluators")]
pub(crate) fn pseudo_threshold_d0(point: &[C<f64>]) -> Option<[C<f64>; 3]> {
    use std::cell::RefCell;
    thread_local! {
        static NATIVE: RefCell<Option<NativeEvaluator<f64>>> = const { RefCell::new(None) };
    }
    let point: &[C<f64>; 11] = point.try_into().ok()?;
    if !d0_four_masses_at_pseudo_threshold(point) {
        return None;
    }
    NATIVE.with(|cell| {
        let mut native = cell.borrow_mut();
        if native.is_none() {
            *native = Some(NativeEvaluator::new(ScalarIntegral::D0).ok()?);
        }
        let values = native.as_ref()?.family_values(point);
        values
            .iter()
            .all(|z| z.re.is_finite() && z.im.is_finite())
            .then_some(values)
    })
}

/// Shared domain contract for the explicit manual backends. Call only after
/// checking shape, so every point has the family's complete argument list.
pub(crate) fn validate_point<T: NativeFloat>(
    family: ScalarIntegral,
    input: &[C<T>],
) -> Result<(), String> {
    if input.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err("all input components must be finite".into());
    }
    let momenta = match family {
        ScalarIntegral::A0 => 0,
        ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
        ScalarIntegral::C0 => 3,
        ScalarIntegral::D0 => 6,
    };
    if input[..momenta].iter().any(|z| !z.im.is_zero()) {
        return Err("external squared invariants must be real".into());
    }
    if input[momenta..input.len() - 1]
        .iter()
        .any(|z| z.im > z.im.zero())
    {
        return Err("squared masses require nonpositive imaginary parts".into());
    }
    let mu = &input[input.len() - 1];
    if !mu.im.is_zero() || mu.re <= mu.re.zero() {
        return Err("mu_squared must be positive and real".into());
    }
    Ok(())
}

#[cfg(all(test, feature = "generated-evaluators"))]
mod complex_mass_tests {
    use super::double_double;
    use crate::{NativeEvaluator, ScalarIntegral};
    use symbolica::domains::float::{Complex, DoubleFloat};

    fn bits(z: &Complex<DoubleFloat>) -> [u64; 4] {
        let (re, im) = (z.re.into_inner(), z.im.into_inner());
        [re.hi().to_bits(), re.lo().to_bits(), im.hi().to_bits(), im.lo().to_bits()]
    }

    /// Masters with complex squared masses (m^2 - i m Gamma, as in the
    /// complex-mass scheme) and with real ones, through the generated
    /// double-double expressions: the fast paths must not change a bit.
    #[test]
    fn complex_mass_masters_are_bitwise_the_full_kernel() {
        let real = |x: f64| Complex::new(DoubleFloat::from(x), DoubleFloat::from(0.0));
        // (mass, width): top, W, Z, Higgs-like, massless and stable.
        let masses = [(1.73, 0.0135), (0.804, 0.0209), (0.912, 0.0250), (1.25, 0.004), (0.0, 0.0), (1.73, 0.0)];
        let squared = |(m, w): (f64, f64)| Complex::new(DoubleFloat::from(m * m), DoubleFloat::from(-m * w));
        let invariant = |k: usize| {
            // Deterministic invariants in [-12, 12], some on thresholds.
            let x = ((k as u64 * 2_654_435_761) % 10_007) as f64 / 10_007.0;
            if k % 7 == 0 { 0.0 } else { 24.0 * x - 12.0 }
        };
        let mut evaluated = 0;
        for (family, momenta, propagators) in [(ScalarIntegral::D0, 6, 4), (ScalarIntegral::C0, 3, 3), (ScalarIntegral::B0, 1, 2)] {
            let mut evaluator = NativeEvaluator::<DoubleFloat>::new(family).expect("double-double evaluator");
            for point in 0..24 {
                let mut input: Vec<_> = (0..momenta).map(|i| real(invariant(point * 11 + i + 1))).collect();
                input.extend((0..propagators).map(|i| squared(masses[(point + 3 * i) % masses.len()])));
                input.push(real(1.0));
                let mut output = vec![real(0.0); 3];
                let mut full = output.clone();
                if evaluator.evaluate(&input, &mut output).is_err() {
                    continue;
                }
                double_double::FULL_KERNEL.with(|f| f.set(true));
                let result = evaluator.evaluate(&input, &mut full);
                double_double::FULL_KERNEL.with(|f| f.set(false));
                result.expect("same domain");
                for (fast, full) in output.iter().zip(&full) {
                    assert_eq!(bits(fast), bits(full), "{family:?} point {point}");
                }
                evaluated += 1;
            }
        }
        assert!(evaluated >= 60, "{evaluated}");
    }
}

#[cfg(all(test, feature = "generated-evaluators"))]
mod threshold_tests {
    use super::d0_next_to_threshold;
    use symbolica::domains::float::Complex;

    fn point(values: [f64; 11]) -> [Complex<f64>; 11] {
        values.map(|x| Complex::new(x, 0.0))
    }

    /// Normal thresholds of two massive propagators only: not their
    /// pseudo-thresholds (equal masses at a vanishing leg), a massless member,
    /// a point 1e-4 away or a mass with a width.
    #[test]
    fn normal_thresholds_of_massive_pairs() {
        const M: f64 = 2.78515625;
        let rest = [-43.0, -32.0, 29.5, -57.25];
        let leg = |p1: f64, [m1, m2, m3, m4]: [f64; 4]| {
            point([p1, M, rest[0], rest[1], rest[2], rest[3], m1, m2, m3, m4, 4.0])
        };
        assert!(d0_next_to_threshold(&leg(4.0 * M, [M, M, 0.0, 0.0])));
        assert!(d0_next_to_threshold(&leg(4.0 * M * (1.0 - 1e-9), [M, M, 0.0, 0.0])));
        assert!(!d0_next_to_threshold(&leg(4.0 * M * (1.0 + 1e-4), [M, M, 0.0, 0.0])));
        assert!(!d0_next_to_threshold(&leg(0.0, [M, M, 0.0, 0.0])));
        assert!(!d0_next_to_threshold(&leg(4.0 * M, [M, 0.0, M, 0.0])));
        // The channel s12 joins propagators 1 and 3.
        let channel = point([1.0, 2.0, 3.0, 5.0, 4.0 * M, -7.0, M, 0.0, M, 0.0, 4.0]);
        assert!(d0_next_to_threshold(&channel));
        let mut wide = leg(4.0 * M, [M, M, 0.0, 0.0]);
        wide[7].im = -0.1;
        assert!(!d0_next_to_threshold(&wide));
    }
}
