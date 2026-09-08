//! Portable SymJIT caches and row-major batched numerical evaluation.
use crate::{ScalarIntegral, masters};
use std::sync::{Mutex, OnceLock};
use symbolica::{
    domains::{float::Complex, rational::Rational},
    evaluate::{ExpressionEvaluator, JITCompilationSettings, JITCompiledEvaluator},
};

type C = Complex<f64>;

pub(crate) fn portable_environment() -> Result<(), String> {
    // Config::default accepts architecture-pinned compiler types through this
    // override, and SymJIT stores that choice in its serialized IR. Require the
    // ordinary Native default for this explicitly portable wrapper.
    if std::env::var_os("SYMJIT_TOML").is_some() {
        return Err(
            "unset SYMJIT_TOML when building or loading portable OneLOop evaluators".into(),
        );
    }
    if std::fs::exists("symjit.toml")
        .map_err(|error| format!("cannot check the default symjit.toml override: {error}"))?
    {
        return Err(
            "remove the working-directory symjit.toml override when building or loading portable OneLOop evaluators".into(),
        );
    }
    Ok(())
}

// Deliberately version-bound: portable IR is not a stable cross-version ABI.
const FORMAT: &[u8] = b"oneloop-jit-v1:fb845d34bda8ccf1fedef6544d3aa46dc24944e3:native-fixes-v1:symjit-2.24.1-patched-v3:strict-o2-simd1\0";

fn cache_payload(data: &[u8]) -> Result<(&[u8], usize, usize), String> {
    let payload = data
        .strip_prefix(FORMAT)
        .ok_or("incompatible evaluator cache")?;
    let dimensions = payload.get(..16).ok_or("truncated evaluator cache")?;
    let dimension = |offset| {
        usize::try_from(u64::from_le_bytes(
            dimensions[offset..offset + 8].try_into().unwrap(),
        ))
        .map_err(|_| "evaluator dimension does not fit this platform".to_owned())
    };
    Ok((&payload[16..], dimension(0)?, dimension(8)?))
}

/// Branch-preserving SymJIT O2 settings with row-major SIMD batching.
/// The patched backend falls back to scalar lanes on divergent branches.
pub fn jit_settings() -> JITCompilationSettings {
    JITCompilationSettings::new()
        .optimization_level(2)
        .direct_translation(true)
        .with_option("fastmath", "false")
        .with_option("fast_complex", "false")
        .with_option("use_threads", "false")
        .with_option("use_simd", "true")
        .with_option("simd_branch", "false")
        .with_option("enable_simd512", "false")
}

/// A reusable evaluator for arbitrary combinations of native expressions.
/// Inputs and outputs of a batch are row-major, one complete point per row.
#[derive(Clone)]
pub struct JitEvaluator {
    inner: JITCompiledEvaluator<C>,
    inputs: usize,
    outputs: usize,
}

impl JitEvaluator {
    pub(crate) fn compile(
        exact: &ExpressionEvaluator<Complex<Rational>>,
        inputs: usize,
        outputs: usize,
    ) -> Result<Self, String> {
        portable_environment()?;
        let inner = exact.jit_compile::<C>(jit_settings())?;
        if inner.input_count() != inputs || inner.output_count() != outputs {
            return Err(
                "compiled evaluator dimensions disagree with requested expression shape".into(),
            );
        }
        Ok(Self {
            inner,
            inputs,
            outputs,
        })
    }

    /// Number of scalar complex arguments per point.
    pub fn input_count(&self) -> usize {
        self.inputs
    }

    /// Number of complex outputs per point.
    pub fn output_count(&self) -> usize {
        self.outputs
    }

    /// Evaluate one point, checking the exact input and output dimensions.
    pub fn evaluate(&mut self, args: &[C], output: &mut [C]) -> Result<(), String> {
        self.validate(args.len(), output.len(), 1)?;
        self.inner.evaluate(args, output);
        Ok(())
    }

    /// Evaluate `rows` points from a flat row-major input array. The output is
    /// also row-major. Empty batches are valid and do not enter the backend.
    pub fn evaluate_batch(
        &mut self,
        args: &[C],
        output: &mut [C],
        rows: usize,
    ) -> Result<(), String> {
        self.validate(args.len(), output.len(), rows)?;
        if rows != 0 {
            self.inner.batch_evaluate(args, output, rows);
        }
        Ok(())
    }

    fn validate(&self, args: usize, output: usize, rows: usize) -> Result<(), String> {
        let expected_in = rows.checked_mul(self.inputs).ok_or("input size overflow")?;
        let expected_out = rows
            .checked_mul(self.outputs)
            .ok_or("output size overflow")?;
        if args != expected_in || output != expected_out {
            return Err(format!(
                "expected {expected_in} inputs and {expected_out} outputs for {rows} rows; received {args} and {output}"
            ));
        }
        Ok(())
    }

    /// Serialize portable IR together with every nested function definition.
    /// Treat this as trusted executable content, not an untrusted data format.
    pub fn to_bytes(&self) -> Result<Vec<u8>, String> {
        let mut data = FORMAT.to_vec();
        data.extend_from_slice(&(self.inputs as u64).to_le_bytes());
        data.extend_from_slice(&(self.outputs as u64).to_le_bytes());
        data.extend(self.inner.export_portable()?);
        Ok(data)
    }

    /// Load a version-matched trusted cache. Native machine code is regenerated
    /// for the host; no symbolic expression construction is required.
    /// Register any custom numerical Symbol callbacks before loading. The
    /// `SYMJIT_TOML` must be unset and the working directory must not contain
    /// `symjit.toml`, since these can override the portable compiler settings.
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        crate::initialize()?;
        Self::from_bytes_raw(data)
    }

    fn from_bytes_raw(data: &[u8]) -> Result<Self, String> {
        portable_environment()?;
        let (payload, inputs, outputs) = cache_payload(data)?;
        let inner = JITCompiledEvaluator::import_portable(payload, jit_settings())?;
        // Dimensions in the outer header cannot authorize shorter slices than
        // the executable consumes. SymJIT enters raw generated machine code.
        if inner.input_count() != inputs || inner.output_count() != outputs {
            return Err("evaluator cache dimensions disagree with executable".into());
        }
        Ok(Self {
            inputs,
            outputs,
            inner,
        })
    }
}

/// Three-coefficient evaluator for one scalar family. The last argument of
/// every point is the positive real squared renormalization scale `mu_squared`.
#[derive(Clone)]
pub struct ScalarEvaluator {
    family: ScalarIntegral,
    evaluator: JitEvaluator,
}

impl ScalarEvaluator {
    /// Clone the eagerly loaded embedded backend, without decoding or compiling
    /// its intermediate code again. The returned instance has its own workspace.
    pub fn prebuilt(family: ScalarIntegral) -> Result<Self, String> {
        let _ = embedded(family)?;
        Self::cached(family)
    }

    /// Clone a shared backend prepared at startup. Without `prebuilt`, startup
    /// builds all five backends from the native expressions instead.
    pub fn cached(family: ScalarIntegral) -> Result<Self, String> {
        cache(family)?
            .lock()
            .map_err(|_| "OneLOop JIT cache poisoned".into())
            .map(|backend| backend.evaluator.clone())
    }

    /// Rebuild O2 code from the current transparent native expressions.
    /// As with manual expression construction, use an adequately sized stack.
    pub fn rebuild(family: ScalarIntegral) -> Result<Self, String> {
        crate::initialize()?;
        Self::rebuild_raw(family)
    }

    fn rebuild_raw(family: ScalarIntegral) -> Result<Self, String> {
        portable_environment()?;
        Ok(Self {
            family,
            evaluator: JitEvaluator::compile(masters::exact_evaluator(family), family.arity(), 3)?,
        })
    }

    /// Integral family represented by this evaluator.
    pub fn family(&self) -> ScalarIntegral {
        self.family
    }

    /// Evaluate a single point into `[finite, simple pole, double pole]`.
    pub fn evaluate(&mut self, args: &[C], output: &mut [C]) -> Result<(), String> {
        self.evaluator.evaluate(args, output)
    }

    /// Evaluate row-major points, writing three consecutive coefficients per row.
    pub fn evaluate_batch(
        &mut self,
        args: &[C],
        output: &mut [C],
        rows: usize,
    ) -> Result<(), String> {
        self.evaluator.evaluate_batch(args, output, rows)
    }

    /// Serialize this family cache, including its family and dimension checks.
    pub fn to_bytes(&self) -> Result<Vec<u8>, String> {
        let mut bytes = vec![self.family as u8];
        bytes.extend(self.evaluator.to_bytes()?);
        Ok(bytes)
    }

    /// Load trusted bytes for a specific family; rejects mismatched metadata.
    pub fn from_bytes(family: ScalarIntegral, bytes: &[u8]) -> Result<Self, String> {
        crate::initialize()?;
        Self::from_bytes_raw(family, bytes)
    }

    fn from_bytes_raw(family: ScalarIntegral, bytes: &[u8]) -> Result<Self, String> {
        if bytes.first().copied() != Some(family as u8) {
            return Err("wrong scalar family in evaluator cache".into());
        }
        let (_, inputs, outputs) = cache_payload(&bytes[1..])?;
        if inputs != family.arity() || outputs != 3 {
            return Err("wrong dimensions in scalar evaluator cache".into());
        }
        let evaluator = JitEvaluator::from_bytes_raw(&bytes[1..])?;
        if evaluator.input_count() != family.arity() || evaluator.output_count() != 3 {
            return Err("wrong dimensions in scalar evaluator cache".into());
        }
        Ok(Self { family, evaluator })
    }
}

#[cfg(feature = "prebuilt")]
fn embedded(family: ScalarIntegral) -> Result<&'static [u8], String> {
    Ok(match family {
        ScalarIntegral::A0 => include_bytes!("../assets/evaluators/a0.bin"),
        ScalarIntegral::B0 => include_bytes!("../assets/evaluators/b0.bin"),
        ScalarIntegral::DB0 => include_bytes!("../assets/evaluators/db0.bin"),
        ScalarIntegral::C0 => include_bytes!("../assets/evaluators/c0.bin"),
        ScalarIntegral::D0 => include_bytes!("../assets/evaluators/d0.bin"),
    })
}

#[cfg(not(feature = "prebuilt"))]
fn embedded(_family: ScalarIntegral) -> Result<&'static [u8], String> {
    Err("embedded evaluators are disabled; use ScalarEvaluator::rebuild".into())
}

// A master exposes scalar Laurent tags, while its backend computes all three
// coefficients together. Keep four FIFO groups for tag-major SIMD callbacks:
// all lanes request one tag before requesting the next. Each tag is consumed
// once per group, including when several lanes have identical arguments.
// This avoids three backend calls per point without memoizing repeated points
// in throughput surveys. Bitwise keys preserve signed branch-axis zeros.
struct CoefficientGroup {
    arguments: [C; 11],
    values: [C; 3],
    used_tags: u8,
}

impl CoefficientGroup {
    fn new() -> Self {
        Self {
            arguments: [C::new(0., 0.); 11],
            values: [C::new(0., 0.); 3],
            used_tags: 0,
        }
    }

    fn reusable(&self, arguments: &[C], tag: usize) -> bool {
        self.used_tags != 0
            && self.used_tags & (1 << tag) == 0
            && arguments
                .iter()
                .zip(&self.arguments)
                .all(|(a, b)| a.re.to_bits() == b.re.to_bits() && a.im.to_bits() == b.im.to_bits())
    }

    fn invalidate(&mut self) {
        self.used_tags = 0;
    }
}

struct CoefficientGroups {
    groups: [CoefficientGroup; 4],
    oldest: usize,
}

impl CoefficientGroups {
    fn new() -> Self {
        Self {
            groups: core::array::from_fn(|_| CoefficientGroup::new()),
            oldest: 0,
        }
    }

    fn invalidate(&mut self) {
        for group in &mut self.groups {
            group.invalidate();
        }
        self.oldest = 0;
    }

    fn evaluate(&mut self, arguments: &[C], tag: usize, evaluate: impl FnOnce(&mut [C; 3])) -> C {
        for offset in 0..self.groups.len() {
            let index = (self.oldest + offset) % self.groups.len();
            let group = &mut self.groups[index];
            if group.reusable(arguments, tag) {
                group.used_tags |= 1 << tag;
                return group.values[tag];
            }
        }

        let index = self.oldest;
        self.oldest = (self.oldest + 1) % self.groups.len();
        let group = &mut self.groups[index];
        // Invalidate the evicted entry before calling the backend. If it
        // panics, this entry stays invalid and the shared mutex is poisoned;
        // no partially updated output is published. Other entries are intact.
        group.invalidate();
        let mut output = [C::new(0., 0.); 3];
        evaluate(&mut output);
        group.arguments[..arguments.len()].copy_from_slice(arguments);
        group.values = output;
        group.used_tags = 1 << tag;
        group.values[tag]
    }
}

struct PreparedBackend {
    evaluator: ScalarEvaluator,
    group: CoefficientGroups,
}

type Backends = [Mutex<PreparedBackend>; 5];
static CACHED: OnceLock<Backends> = OnceLock::new();

pub(crate) fn initialize_all() -> Result<(), String> {
    let mut backends = Vec::with_capacity(5);
    for family in [
        ScalarIntegral::A0,
        ScalarIntegral::B0,
        ScalarIntegral::DB0,
        ScalarIntegral::C0,
        ScalarIntegral::D0,
    ] {
        #[cfg(feature = "prebuilt")]
        let evaluator = ScalarEvaluator::from_bytes_raw(family, embedded(family)?);
        #[cfg(not(feature = "prebuilt"))]
        let evaluator = ScalarEvaluator::rebuild_raw(family);
        backends.push(Mutex::new(PreparedBackend {
            evaluator: evaluator.map_err(|error| {
                format!("could not initialize {} evaluator: {error}", family.name())
            })?,
            group: CoefficientGroups::new(),
        }));
    }
    let backends = backends
        .try_into()
        .map_err(|_| "incorrect OneLOop backend inventory")?;
    CACHED
        .set(backends)
        .map_err(|_| "OneLOop backends were already initialized".into())
}

fn cache(family: ScalarIntegral) -> Result<&'static Mutex<PreparedBackend>, String> {
    crate::initialize()?;
    Ok(&CACHED.get().ok_or("OneLOop backends are not initialized")?[family as usize])
}

/// Evaluate one point using the already prepared shared family backend.
/// Inputs end in `mu_squared`; outputs are finite, simple pole and double pole.
pub fn evaluate(family: ScalarIntegral, args: &[C], output: &mut [C]) -> Result<(), String> {
    let mut backend = cache(family)?
        .lock()
        .map_err(|_| "OneLOop JIT cache poisoned")?;
    backend.group.invalidate();
    backend.evaluator.evaluate(args, output)
}

/// Evaluate flat row-major points using the already prepared shared backend.
/// No expression construction or evaluator compilation takes place here.
pub fn evaluate_batch(
    family: ScalarIntegral,
    args: &[C],
    output: &mut [C],
    rows: usize,
) -> Result<(), String> {
    let mut backend = cache(family)?
        .lock()
        .map_err(|_| "OneLOop JIT cache poisoned")?;
    backend.group.invalidate();
    backend.evaluator.evaluate_batch(args, output, rows)
}

pub(crate) fn evaluate_cached(family: ScalarIntegral, tag: usize, args: &[C]) -> C {
    let mut backend = cache(family)
        .expect("initialized OneLOop backends")
        .lock()
        .expect("OneLOop JIT cache poisoned");
    let PreparedBackend { evaluator, group } = &mut *backend;
    group.evaluate(args, tag, |output| {
        evaluator
            .evaluate(args, output)
            .expect("validated master arguments");
    })
}

/// Explicitly rebuild and replace the evaluator used by a family's f64 Symbol
/// hooks. Other existing manual evaluators are unaffected.
pub fn rebuild_cached_evaluator(family: ScalarIntegral) -> Result<(), String> {
    let evaluator = ScalarEvaluator::rebuild(family)?;
    let mut backend = cache(family)?
        .lock()
        .map_err(|_| "OneLOop JIT cache poisoned")?;
    backend.evaluator = evaluator;
    backend.group.invalidate();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counted_evaluate(
        group: &mut CoefficientGroups,
        arguments: &[C],
        tag: usize,
        calls: &mut usize,
    ) -> C {
        group.evaluate(arguments, tag, |output| {
            *calls += 1;
            for (index, value) in output.iter_mut().enumerate() {
                *value = C::new(*calls as f64, index as f64);
            }
        })
    }

    #[test]
    fn coefficient_groups_execute_once_per_point_in_every_tag_order() {
        let orders = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        let args = [C::new(2., 0.), C::new(1., 0.)];
        for first in orders {
            for second in orders {
                let mut group = CoefficientGroups::new();
                let mut calls = 0;
                for (point, order) in [first, second].into_iter().enumerate() {
                    for tag in order {
                        assert_eq!(
                            counted_evaluate(&mut group, &args, tag, &mut calls),
                            C::new((point + 1) as f64, tag as f64),
                            "orders={first:?}/{second:?}, point={point}, tag={tag}"
                        );
                        assert_eq!(calls, point + 1);
                    }
                }
            }
        }
    }

    #[test]
    fn coefficient_groups_repeated_tags_always_start_another_execution() {
        let args = [C::new(2., 0.), C::new(1., 0.)];
        for tag in 0..3 {
            let mut group = CoefficientGroups::new();
            let mut calls = 0;
            for expected_calls in 1..=4 {
                assert_eq!(
                    counted_evaluate(&mut group, &args, tag, &mut calls),
                    C::new(expected_calls as f64, tag as f64)
                );
                assert_eq!(calls, expected_calls);
            }
            for other_tag in (0..3).filter(|other| *other != tag) {
                assert_eq!(
                    counted_evaluate(&mut group, &args, other_tag, &mut calls),
                    C::new(1., other_tag as f64)
                );
                assert_eq!(calls, 4);
            }
        }
    }

    #[test]
    fn coefficient_groups_execute_again_for_signed_zeros_and_changed_points() {
        let args = [C::new(0., 0.); 11];
        // Exercise both component keys at every argument position, including
        // the last slot used by the largest scalar family.
        for index in 0..args.len() {
            for imaginary in [false, true] {
                let mut group = CoefficientGroups::new();
                let mut calls = 0;
                counted_evaluate(&mut group, &args, 0, &mut calls);
                let mut signed_axis = args;
                if imaginary {
                    signed_axis[index].im = -0.;
                } else {
                    signed_axis[index].re = -0.;
                }
                assert_eq!(
                    counted_evaluate(&mut group, &signed_axis, 1, &mut calls),
                    C::new(2., 1.)
                );
                assert_eq!(calls, 2);
                assert_eq!(
                    counted_evaluate(&mut group, &args, 2, &mut calls),
                    C::new(1., 2.)
                );
                assert_eq!(calls, 2);
                let mut other_point = args;
                other_point[index].re = 4.;
                assert_eq!(
                    counted_evaluate(&mut group, &other_point, 0, &mut calls),
                    C::new(3., 0.)
                );
                assert_eq!(calls, 3);
            }
        }
    }

    #[test]
    fn coefficient_groups_invalidate_explicitly_and_before_backend_failure() {
        let args = [C::new(2., 0.), C::new(1., 0.)];
        let mut group = CoefficientGroups::new();
        let mut calls = 0;
        counted_evaluate(&mut group, &args, 0, &mut calls);
        group.invalidate();
        assert_eq!(
            counted_evaluate(&mut group, &args, 1, &mut calls),
            C::new(2., 1.)
        );
        assert_eq!(calls, 2);

        // Fill the ring, leaving the previous point in the next victim slot.
        for index in 1..4 {
            let mut lane = args;
            lane[0].re += index as f64;
            counted_evaluate(&mut group, &lane, 0, &mut calls);
        }
        assert_eq!(calls, 5);

        let mut other_point = args;
        other_point[1].re = 4.;
        let failure = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            group.evaluate(&other_point, 2, |output| {
                calls += 1;
                output[0] = C::new(123., 456.);
                panic!("simulated backend failure after partial output");
            });
        }));
        assert!(failure.is_err());
        assert_eq!(calls, 6);
        assert_eq!(group.groups[0].used_tags, 0);
        assert!(group.groups[1..].iter().all(|entry| entry.used_tags == 1));
        // Tag 2 of the old point was unused: reusing it here would conceal a
        // missing invalidation on the failed transition.
        assert_eq!(
            counted_evaluate(&mut group, &args, 2, &mut calls),
            C::new(7., 2.)
        );
        assert_eq!(calls, 7);
    }

    #[test]
    fn coefficient_groups_tag_major_lanes_reuse_oldest_matching_point() {
        let orders = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        for lanes in [1, 2, 4] {
            for identical in [false, true] {
                for first in orders {
                    for second in orders {
                        let mut group = CoefficientGroups::new();
                        let mut calls = 0;
                        for (batch, order) in [first, second].into_iter().enumerate() {
                            for tag in order {
                                for lane in 0..lanes {
                                    let args = [
                                        C::new(if identical { 2. } else { 2. + lane as f64 }, 0.),
                                        C::new(1., 0.),
                                    ];
                                    assert_eq!(
                                        counted_evaluate(&mut group, &args, tag, &mut calls),
                                        C::new((batch * lanes + lane + 1) as f64, tag as f64),
                                        "lanes={lanes}, identical={identical}, orders={first:?}/{second:?}, batch={batch}, lane={lane}, tag={tag}"
                                    );
                                }
                                assert_eq!(calls, (batch + 1) * lanes);
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn coefficient_groups_evict_fifo_without_reusing_consumed_tags() {
        let mut group = CoefficientGroups::new();
        let mut calls = 0;
        for point in 0..5 {
            counted_evaluate(&mut group, &[C::new(point as f64, 0.)], 0, &mut calls);
        }
        assert_eq!(calls, 5);
        // Point zero was evicted; its missing tag must execute again.
        assert_eq!(
            counted_evaluate(&mut group, &[C::new(0., 0.)], 1, &mut calls),
            C::new(6., 1.)
        );
        assert_eq!(calls, 6);
        // Point two is still buffered and its tag remains unused.
        assert_eq!(
            counted_evaluate(&mut group, &[C::new(2., 0.)], 1, &mut calls),
            C::new(3., 1.)
        );
        assert_eq!(calls, 6);
        assert_eq!(
            counted_evaluate(&mut group, &[C::new(2., 0.)], 1, &mut calls),
            C::new(7., 1.)
        );
        assert_eq!(calls, 7);
    }

    #[test]
    fn coefficient_groups_manual_invalidation_clears_all_lanes() {
        for identical in [false, true] {
            let mut group = CoefficientGroups::new();
            let mut calls = 0;
            let args = core::array::from_fn::<_, 4, _>(|lane| {
                [C::new(if identical { 2. } else { lane as f64 }, 0.)]
            });
            for lane in &args {
                counted_evaluate(&mut group, lane, 0, &mut calls);
            }
            group.invalidate();
            assert!(group.groups.iter().all(|entry| entry.used_tags == 0));
            for (lane, arguments) in args.iter().enumerate() {
                assert_eq!(
                    counted_evaluate(&mut group, arguments, 1, &mut calls),
                    C::new((5 + lane) as f64, 1.)
                );
            }
            assert_eq!(calls, 8);
        }
    }
}
