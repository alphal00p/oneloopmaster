//! Public scalar master Symbols and their numeric evaluation hooks.
//!
//! A call has the form `oneloopmaster::M(tag, arguments...)`, with Laurent
//! tags `0`, `-1`, and `-2`. Numeric arguments follow the corresponding
//! lowercase constructor, with the squared renormalization scale last.
//! Native definitions supplied by [`OneLoopExpressions`] take precedence
//! over the hooks and remain visible to Symbolica's evaluator compiler.

use super::*;
use std::{
    cell::Cell,
    collections::BTreeMap,
    sync::{Mutex, OnceLock},
};
use symbolica::{
    atom::{AtomView, EvaluationInfo},
    domains::{
        float::{Complex, DoubleFloat, Float},
        rational::Rational,
    },
    evaluate::ExpressionEvaluator,
};

/// Scalar integral family. Each evaluator returns finite, simple-pole and
/// double-pole coefficients, with `mu_squared` as its final numeric argument.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(usize)]
pub enum ScalarIntegral {
    A0,
    B0,
    DB0,
    C0,
    D0,
}

type Family = ScalarIntegral;

impl ScalarIntegral {
    /// Conventional master name.
    pub fn name(self) -> &'static str {
        ["A0", "B0", "dB0", "C0", "D0"][self as usize]
    }

    /// Number of numeric arguments, including the squared scale but not the
    /// Laurent tag used by the Symbol interface.
    pub fn arity(self) -> usize {
        [2, 4, 4, 7, 11][self as usize]
    }

    fn series(self, a: &[Atom]) -> LaurentSeries {
        match self {
            Self::A0 => a0(&a[0], &a[1]),
            Self::B0 => b0(&a[0], &a[1], &a[2], &a[3]),
            Self::DB0 => db0(&a[0], &a[1], &a[2], &a[3]),
            Self::C0 => {
                expressions::triangle_series([&a[0], &a[1], &a[2]], [&a[3], &a[4], &a[5]], &a[6])
            }
            Self::D0 => expressions::box_series(
                [&a[0], &a[1], &a[2], &a[3], &a[4], &a[5]],
                [&a[6], &a[7], &a[8], &a[9]],
                &a[10],
            ),
        }
    }

    fn parameters(self) -> Vec<Symbol> {
        (0..self.arity())
            .map(|index| symbolica::symbol!(format!("__olo_master_{}_{}", self.name(), index)))
            .collect()
    }

    fn check_arity(self, count: usize) {
        assert_eq!(
            count,
            self.arity(),
            "oneloopmaster::{} expects {} numeric arguments after its Laurent tag; received {}",
            self.name(),
            self.arity(),
            count
        );
    }
}

fn coefficient_index(family: Family, tags: &[AtomView<'_>]) -> usize {
    assert_eq!(
        tags.len(),
        1,
        "oneloopmaster::{} requires one Laurent tag",
        family.name()
    );
    if tags[0] == 0 {
        0
    } else if tags[0] == -1 {
        1
    } else if tags[0] == -2 {
        2
    } else {
        panic!(
            "invalid oneloopmaster::{} Laurent tag {}; expected 0, -1, or -2",
            family.name(),
            tags[0]
        );
    }
}

// A hook only evaluates native definitions, which never call a public master.
// Detect accidental recursion before acquiring any cache or OnceLock: this
// turns a future dependency cycle into a useful error rather than a deadlock.
thread_local! { static IN_MASTER: Cell<bool> = const { Cell::new(false) }; }
struct EvaluationGuard;
impl EvaluationGuard {
    fn enter() -> Self {
        IN_MASTER.with(|active| {
            assert!(!active.get(), "reentrant OneLOop master evaluation");
            active.set(true);
        });
        Self
    }
}
impl Drop for EvaluationGuard {
    fn drop(&mut self) {
        IN_MASTER.with(|active| active.set(false));
    }
}

type ExactEvaluator = ExpressionEvaluator<Complex<Rational>>;
type PrecisionEvaluator = crate::NativeEvaluator<Float>;

pub(crate) fn exact_evaluator(family: ScalarIntegral) -> &'static ExactEvaluator {
    initialization::ensure_symbolica_state();
    static CACHE: [OnceLock<ExactEvaluator>; 5] = [const { OnceLock::new() }; 5];
    CACHE[family as usize].get_or_init(|| {
        let args = family
            .parameters()
            .into_iter()
            .map(Symbol::to_atom)
            .collect::<Vec<_>>();
        OneLoopExpressions::new()
            .evaluator(family.series(&args).coefficients(), &args)
            .expect("compile native OneLOop master definitions")
    })
}

pub(crate) fn jit_evaluator(family: ScalarIntegral) -> Result<JitEvaluator, String> {
    let args = family
        .parameters()
        .into_iter()
        .map(Symbol::to_atom)
        .collect::<Vec<_>>();
    OneLoopExpressions::new().jit_evaluator(family.series(&args).coefficients(), &args)
}

fn scalar_f64(family: Family, tag: usize, args: &[Complex<f64>]) -> Complex<f64> {
    family.check_arity(args.len());
    let _guard = EvaluationGuard::enter();
    crate::backend::native_coefficient(family, tag, args)
}

fn scalar_double_float(
    family: Family,
    tag: usize,
    args: &[Complex<DoubleFloat>],
) -> Complex<DoubleFloat> {
    family.check_arity(args.len());
    let _guard = EvaluationGuard::enter();
    static CACHE: [OnceLock<Mutex<crate::NativeEvaluator<DoubleFloat>>>; 5] =
        [const { OnceLock::new() }; 5];
    let mut evaluator = CACHE[family as usize]
        .get_or_init(|| {
            Mutex::new(crate::NativeEvaluator::new(family).expect("native DoubleFloat constants"))
        })
        .lock()
        .expect("OneLOop DoubleFloat cache poisoned");
    let mut output = core::array::from_fn::<_, 3, _>(|_| Complex::new(0.0.into(), 0.0.into()));
    let result = evaluator.evaluate(args, &mut output);
    drop(evaluator);
    // Domain rejection must not poison a reusable shared evaluator when the
    // caller catches the numeric callback's deliberate panic.
    result.expect("valid master arguments");
    output[tag]
}

fn scalar_float(family: Family, tag: usize, args: &[Complex<Float>]) -> Complex<Float> {
    family.check_arity(args.len());
    let _guard = EvaluationGuard::enter();
    // The requested precision reaches an external function in its arguments;
    // include every real and imaginary component, including the scale.
    let precision = args
        .iter()
        .flat_map(|value| [value.re.prec(), value.im.prec()])
        .max()
        .expect("every master has numeric arguments");
    static CACHE: [OnceLock<Mutex<BTreeMap<u32, PrecisionEvaluator>>>; 5] =
        [const { OnceLock::new() }; 5];
    let mut cache = CACHE[family as usize]
        .get_or_init(|| Mutex::new(BTreeMap::new()))
        .lock()
        .expect("OneLOop precision evaluator cache poisoned");
    let evaluator = cache.entry(precision).or_insert_with(|| {
        crate::NativeEvaluator::with_binary_precision(family, precision)
            .expect("native arbitrary-precision constants")
    });
    let mut output = core::array::from_fn::<_, 3, _>(|_| {
        Complex::new(Float::new(precision), Float::new(precision))
    });
    let result = evaluator.evaluate(args, &mut output);
    drop(cache);
    result.expect("valid master arguments");
    output[tag].clone()
}

fn info(family: Family) -> EvaluationInfo {
    // Symbolica routes a call containing only its tag through the constant
    // interface. Such a call is malformed for every master; report the same
    // arity contract there as in the numeric callbacks.
    EvaluationInfo::constant(move |tags, _precision| {
        coefficient_index(family, tags);
        Err(format!(
            "oneloopmaster::{} expects {} numeric arguments after its Laurent tag; received 0",
            family.name(),
            family.arity()
        ))
    })
    .with_tags(1)
    .register_tagged::<Complex<f64>>(move |tags| {
        let tag = coefficient_index(family, tags);
        Box::new(move |args: &[Complex<f64>]| scalar_f64(family, tag, args))
    })
    .register_tagged::<Complex<Float>>(move |tags| {
        let tag = coefficient_index(family, tags);
        Box::new(move |args: &[Complex<Float>]| scalar_float(family, tag, args))
    })
    .register_tagged::<Complex<DoubleFloat>>(move |tags| {
        let tag = coefficient_index(family, tags);
        Box::new(move |args: &[Complex<DoubleFloat>]| scalar_double_float(family, tag, args))
    })
}

/// Initializes and returns `oneloopmaster::A0(tag, m², mu²)`.
///
/// All master accessors attach direct Rust callbacks for `Complex<f64>`,
/// `Complex<DoubleFloat>` and `Complex<Float>`. Invalid tags, domains or numeric arities panic with a descriptive
/// message because Symbolica's numeric callback interface returns a number,
/// not a `Result`. The native map performs its own function validation.
#[allow(non_snake_case)]
pub fn A0() -> Symbol {
    initialization::ensure_symbolica_state();
    static SYMBOL: OnceLock<Symbol> = OnceLock::new();
    *SYMBOL.get_or_init(|| symbolica::symbol!("oneloopmaster::A0", eval = info(Family::A0)))
}

/// Initializes and returns `oneloopmaster::B0(tag, p², m0², m1², mu²)`.
#[allow(non_snake_case)]
pub fn B0() -> Symbol {
    initialization::ensure_symbolica_state();
    static SYMBOL: OnceLock<Symbol> = OnceLock::new();
    *SYMBOL.get_or_init(|| symbolica::symbol!("oneloopmaster::B0", eval = info(Family::B0)))
}

/// Initializes and returns `oneloopmaster::dB0(tag, p², m0², m1², mu²)`.
#[allow(non_snake_case)]
pub fn dB0() -> Symbol {
    initialization::ensure_symbolica_state();
    static SYMBOL: OnceLock<Symbol> = OnceLock::new();
    *SYMBOL.get_or_init(|| symbolica::symbol!("oneloopmaster::dB0", eval = info(Family::DB0)))
}

/// Initializes and returns `oneloopmaster::C0(tag, p1², p2², p3², m1², m2², m3², mu²)`.
#[allow(non_snake_case)]
pub fn C0() -> Symbol {
    initialization::ensure_symbolica_state();
    static SYMBOL: OnceLock<Symbol> = OnceLock::new();
    *SYMBOL.get_or_init(|| symbolica::symbol!("oneloopmaster::C0", eval = info(Family::C0)))
}

/// Initializes and returns the box master, with the arguments of [`crate::d0`]
/// flattened as six momenta, four squared masses, and the squared scale.
#[allow(non_snake_case)]
pub fn D0() -> Symbol {
    initialization::ensure_symbolica_state();
    static SYMBOL: OnceLock<Symbol> = OnceLock::new();
    *SYMBOL.get_or_init(|| symbolica::symbol!("oneloopmaster::D0", eval = info(Family::D0)))
}

/// Adds the fifteen tagged native master definitions to a complete native map.
/// Pure sector-call constructors keep registration independent of map creation.
pub(crate) fn register_masters(map: &mut FunctionMap) {
    for (family, master) in [
        (Family::A0, A0()),
        (Family::B0, B0()),
        (Family::DB0, dB0()),
        (Family::C0, C0()),
        (Family::D0, D0()),
    ] {
        let vars = family.parameters();
        let args = vars
            .iter()
            .copied()
            .map(Symbol::to_atom)
            .collect::<Vec<_>>();
        for (tag, body) in [0, -1, -2]
            .into_iter()
            .zip(family.series(&args).into_coefficients())
        {
            map.add_tagged_function_with_options(
                master,
                vec![Atom::num(tag)],
                vars.clone(),
                body,
                native_function_options(),
            )
            .expect("unique native OneLOop master definition");
        }
    }
}
