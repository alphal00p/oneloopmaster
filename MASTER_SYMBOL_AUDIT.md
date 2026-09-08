# Public master Symbols: API audit

## Current implementation and release checks — 2026-09-08

`src/masters.rs` now defines the five requested `oneloopmaster::` Symbols, each
with a leading Laurent-power tag and `EvaluationInfo` hooks for `Complex<f64>`
and `Complex<Float>`. The shared function map also registers all fifteen native
coefficient bodies under these same Symbols. Hooks reuse cached evaluators of
those bodies; there is no separate numerical OneLOop implementation. Every
master has the squared renormalization scale as its last numeric argument,
after the leading Laurent tag and the family's kinematic arguments.

The current engine is Symbolica dev `fb845d34bda8ccf1fedef6544d3aa46dc24944e3`
with the dependency patches recorded in [patches/README.md](patches/README.md),
including its native binary64 Li2 specialization. A Symbolica `initialize!`
hook, ordered after `symbolica::special_functions`, registers all master/helper
symbols and prepares all five binary64 SymJIT backends together. Default manual
calls and binary64 master hooks use that same prepared cache. The default
embedded portable IR is restored at startup; explicit rebuilding remains
available. Rust `initialize()` forces startup and returns configuration errors;
Python module import invokes it eagerly.

The generic exact evaluator is cached once per family when needed; numeric
`Complex<Float>` evaluators are cached by requested precision, not eagerly for
every possible precision. Exact rational constants are converted directly at
that precision. Initialization enters Symbolica state before family-local cache
locks, and numeric reentry is rejected before cache locking. Invalid tags and
arities are explicit.

The current release suite passes **57 tests**, with zero failures and four ignored
test entries (three opt-in diagnostics and the subprocess initialization helper),
including all three tests in `tests/masters.rs`. Those three cover all 879 acceptance
coefficients through both hook/map routes and all Laurent tags, constant/variable
calls, independent 128/256-bit elementary controls, mixed precision/tiny-width
inputs, and malformed calls. The previously recorded 395.58-second serial debug
run belongs to the prior engine, not the timing of this release rerun.
An earlier run failed at the three-mass box on
`tests/data/parity.txt:217`, returning NaN instead of
`0.03220049403325557 + 0.07401441373129154i` after enabling the contour route.
An exact zero-endpoint dilogarithm limit fixes that failure in both the focused
31-point three-mass checks and the complete master acceptance rerun. The fixture
test collects all disagreements rather than stopping at the first one.

The added `tests/master_callbacks_jit.rs` regression also passes. It compiles bare
master calls into an outer SymJIT evaluator without a native function map and
compares all three outputs against both the Fortran fixtures and prepared family
backends. Original and portable-restored evaluators cover scalar calls and batch
sizes 1/3/4/5/31, repeated inputs, mixed rows and partial tails. Canonical tag order
covers all **293 acceptance rows**; each of the other five tag orders covers eight
finite/pole controls per family, not the whole table. Family inventory and
executable dimensions are asserted. This is execution-path coverage, not a new
claim about unsampled analytic regions.

No global parity claim follows from the interface being implemented. The expanded
fixed-f64 audit still fails on 14 small-momentum B0/dB0 coefficients, although all
342 expanded points pass at 256 bits. See [SCALAR_PARITY_AUDIT.md](SCALAR_PARITY_AUDIT.md).
Portable evaluator export includes nested native definitions; custom numerical
callbacks must already be registered when importing a manually composed cache.
Expression export without the required definitions or hooks and additional
numeric domains remain separate integration concerns.

Mixed-precision inputs retain their original precision; choosing a cache at the
maximum input-component precision does not restore missing input digits. The
precision cache currently retains an evaluator for each distinct requested
precision without eviction. Accuracy beyond the underlying transcendental
implementation's roughly 900-bit stopping threshold is not certified.
The binary64 backend computes all three coefficients together. Distinct Laurent
tags at bitwise-identical arguments reuse that result; signed zeros remain part
of the key. A repeated tag starts a fresh group even at identical inputs, so
repeated points are not memoized away in throughput measurements. A bounded
four-entry FIFO supports tag-major SIMD callback lanes; different inputs use
different entries. Manual shared-cache calls and explicit replacement invalidate
all entries. Seven counting-backend tests check this behavior, including all tag
permutations, successive 1/2/4-lane groups, FIFO eviction and selected-entry
invalidation before a backend failure. Float hooks still evaluate all three
coefficients on each call. Same-family hooks serialize on their evaluator cache;
interleaved points or threads can reduce binary64 grouping, and an unexpected
numerical panic can poison that cache.
The transparent mapped route remains preferable for whole amplitudes. These
are operational limitations, not additional scalar algorithms or proof of a
formula mismatch.

Prepared scalar JIT, public master calls compiled with the native function map,
and bare master callbacks are distinct benchmark routes. The bare-callback survey
has both an interpreted scalar outer evaluator (`hooks`) and an outer JIT with
batching (`hooks-jit`), without a native map. The latter still invokes each scalar
family callback per lane, rather than merging callbacks into a family SIMD call.
Measured performance does not meet the requested
global 1.5-times-Fortran target; see the dated [performance report](performance/2026-09-08/README.md)
for completed measurements and their scope. Shared-host community runtime
validation remains pending.

The following is the preserved **2026-09-07 baseline audit**, before these
implementation changes; its missing-interface statements describe that baseline.

## Implementation at the 2026-09-07 baseline

- `a0`, `b0` and `db0` in `src/two_point.rs` return expanded native `Atom` bodies.
- `OneLoopExpressions::c0` and `d0` in `src/expressions.rs` return native
  conditional expressions calling internal per-sector Symbols.
- Those sector Symbols have transparent `FunctionMap` definitions. None of the
  public masters has a Symbol-attached `EvaluationInfo` hook.
- Keeping the function map makes the current C0/D0 expressions evaluable;
  exporting only their printed function calls loses the definitions.

## Symbolica source and documentation checked

The inspected dev revision is `0b57776bf911faeea7e28ea133706fb03740ffeb`, with
the separately recorded local patches. The official
[numerical evaluation documentation](https://symbolica.io/docs/numerical_evaluation.html)
distinguishes nested symbolic definitions from external evaluation hooks.

| API | Source in the inspected Symbolica checkout | Meaning |
| --- | --- | --- |
| `FunctionMap::add_tagged_function_with_options` | `src/evaluate/function_map.rs:188` | A Symbol and symbolic tags map to a native `Atom` body. |
| `EvaluationInfo::with_tags` / `register_tagged` | `src/atom.rs:471`, `src/atom.rs:527` | Numeric-domain callbacks specialized by leading symbolic tags. |
| `SymbolBuilder::with_evaluation_info` or `symbol!(..., eval = ...)` | `src/atom.rs:1051` | Attaches the numerical evaluator to the Symbol itself. |
| `SymbolBuilder::with_normalization_function` | `src/atom.rs:934` | Rewrites during normalization; not a substitute for a numeric evaluator. |
| Function-map lookup before evaluation-hook lookup | `src/evaluate/tree.rs:1121` | A native definition takes precedence over the external callback when provided. |

## Recommendation at the 2026-09-07 baseline, since implemented

1. Define once-initialized public A0/B0/dB0/C0/D0 Symbols with a leading
   Laurent-power tag (`0`, `-1`, `-2`). Keep existing formula constructors for
   inspection and compatibility.
2. Register each coefficient's existing native expression under its master
   Symbol in the shared map, using tagged definitions and non-inlined calls.
   This retains compact notation and transparent definitions.
3. If bare-Symbol evaluation is required, attach a dedicated `EvaluationInfo`
   bridge that evaluates those same native definitions. Do not introduce a
   separate Fortran or Rust numerical OneLOop algorithm.
4. Prefer the mapped route for amplitudes: it exposes the native bodies to
   Symbolica. A numeric callback alone is opaque to the outer optimizer and
   does not automatically acquire every numeric/error-propagating domain.

The bridge must support `Complex<f64>` and `Complex<Float>` deliberately. Real
external invariants do not imply a real integral value. Arbitrary-precision
constants must be mapped at the input precision, without an intermediate f64
rounding step. The callback signature has no `Result` channel, so invalid tags,
arity and unsupported numeric domains require explicit API/error decisions.
Avoid reentrant initialization of the definition cache and do not eagerly
expand all master calls in a normalization hook.

This design can reconcile compact master Symbols with the existing native-body
requirement, but a fallback numerical hook is still opaque when used without
the map. It must not be advertised as equivalent to transparent native
evaluation or automatic error-controlled precision.

## Required acceptance checks

- All five families and all three Laurent tags, comparing direct Symbol-hook
  evaluation, mapped evaluation, the existing expression bodies and the
  independent Fortran fixtures.
- Constant-only calls, repeated initialization, repeated/mixed master calls,
  and function-map precedence over the hook.
- Complex masses, timelike inputs with imaginary output, independent scale
  changes, and 128-/256-bit inputs without loss through f64 conversion.
- Invalid tag/arity handling, supported numeric domains, and serialization /
  export behavior when native definitions or external hooks are unavailable.

Adding this interface will not repair the 44 coefficient mismatches documented
in [SCALAR_PARITY_AUDIT.md](SCALAR_PARITY_AUDIT.md). Both audit criteria remain open.
