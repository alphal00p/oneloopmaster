//! Ahead-of-time generic Rust implementations of the native master expressions.
//!
//! The generated functions contain ordinary arithmetic and structured Rust
//! branches. They neither build nor execute a Symbolica/SymJIT evaluator.
mod primitives;
mod simple;
// The emitter has its own stable formatting and mechanically named SSA locals.
#[rustfmt::skip]
#[allow(clippy::all)]
mod generated_a0;
#[rustfmt::skip]
#[allow(clippy::all)]
mod generated_b0;
#[rustfmt::skip]
#[allow(clippy::all)]
mod generated_db0;
#[rustfmt::skip]
#[allow(clippy::all)]
mod generated_c0;
#[rustfmt::skip]
#[allow(clippy::all)]
mod generated_d0;
#[rustfmt::skip]
#[allow(clippy::all)]
mod generated_constants;

use crate::ScalarIntegral;
pub use primitives::NativeFloat;
use primitives::*;
use symbolica::domains::{integer::Integer, rational::Rational};

/// Exact initialization data emitted alongside the generated arithmetic.
/// No strings or rational conversions are used in a warmed scalar call.
#[allow(dead_code)] // The generator also supports fixed Li2 constants in future formulas.
enum ConstantSpec {
    Rational { re: &'static str, im: &'static str },
    Pi,
    Polylog2 { re: &'static str, im: &'static str },
}

#[derive(Clone)]
struct Context<T: NativeFloat> {
    constants: Vec<C<T>>,
    simple: simple::Constants<T>,
}

impl<T: NativeFloat> Context<T> {
    fn new(bits: u32) -> Result<Self, String> {
        let prototype = T::zero_at(bits);
        let rational = |text: &str| -> Result<T, String> {
            let (num, den) = text.split_once('/').unwrap_or((text, "1"));
            let num = num.parse::<Integer>().map_err(|e| e.to_string())?;
            let den = den.parse::<Integer>().map_err(|e| e.to_string())?;
            if den == 0 {
                return Err("invalid generated rational denominator".into());
            }
            Ok(prototype.from_rational(&Rational::from((num, den))))
        };
        let constants = generated_constants::CONSTANTS
            .iter()
            .map(|spec| {
                Ok(match spec {
                    ConstantSpec::Rational { re, im } => C::new(rational(re)?, rational(im)?),
                    ConstantSpec::Pi => C::new(prototype.pi(), prototype.zero()),
                    ConstantSpec::Polylog2 { re, im } => {
                        polylog2(&C::new(rational(re)?, rational(im)?))
                    }
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(Self {
            constants,
            simple: simple::Constants::new(&prototype),
        })
    }

    #[inline(always)]
    fn constant(&self, index: usize) -> &C<T> {
        &self.constants[index]
    }
}

/// A reusable direct Rust scalar master implementation.
///
/// `T` may be `f64`, Symbolica's `DoubleFloat`, or arbitrary-precision `Float`.
/// Inputs use the usual squared quantities and have `mu_squared` last; outputs
/// are `[finite, simple pole, double pole]`. This backend uses no runtime
/// evaluator or JIT. Constants are prepared once at construction.
#[derive(Clone)]
pub struct NativeEvaluator<T: NativeFloat> {
    family: ScalarIntegral,
    bits: u32,
    context: Context<T>,
    arguments: Vec<C<T>>,
}

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
            bits,
            context: Context::new(bits)?,
            arguments: (0..family.arity())
                .map(|_| C::new(T::zero_at(bits), T::zero_at(bits)))
                .collect(),
        })
    }

    /// Integral family represented by this implementation.
    pub fn family(&self) -> ScalarIntegral {
        self.family
    }

    /// Working precision in bits (not decimal digits).
    pub fn binary_precision(&self) -> u32 {
        self.bits
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
        let input = if T::FIXED_BITS.is_some() {
            input
        } else {
            for (working, value) in self.arguments.iter_mut().zip(input) {
                working.re = value.re.at_precision(self.bits);
                working.im = value.im.at_precision(self.bits);
            }
            &self.arguments
        };
        let values = match self.family {
            ScalarIntegral::A0 if !truth(&input[0]) => {
                let mut values = simple::zeros(&input[0].re);
                values[1] = input[0].clone();
                values
            }
            ScalarIntegral::A0 => generated_a0::evaluate(input, &self.context),
            ScalarIntegral::B0 if input[..3].iter().all(|z| !truth(z)) => {
                simple::zeros(&input[0].re)
            }
            ScalarIntegral::B0 => generated_b0::evaluate(input, &self.context),
            ScalarIntegral::DB0 => generated_db0::evaluate(input, &self.context),
            ScalarIntegral::C0 => simple::massless_triangle(input, &self.context.simple)
                .unwrap_or_else(|| generated_c0::evaluate(input, &self.context)),
            ScalarIntegral::D0 => simple::massless_onshell_box(input, &self.context.simple)
                .unwrap_or_else(|| generated_d0::evaluate(input, &self.context)),
        };
        for (out, value) in output.iter_mut().zip(values) {
            *out = value;
        }
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

    fn validate_point(&self, input: &[C<T>]) -> Result<(), String> {
        validate_point(self.family, input)
    }
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
