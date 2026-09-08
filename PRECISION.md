# Arbitrary-precision evaluation

The same exact Symbolica expression graph is available as a native expression
interpreter and as generated generic Rust arithmetic. `NativeEvaluator<T>`
supports `f64`, `DoubleFloat` and `Float`; `PrecisionEvaluator` (Rust) and the
Python `prec` keyword provide decimal-precision `Complex<Float>` interfaces.
The direct native route is a transcription of that graph, not an alternative
numerical OneLOop algorithm, a Fortran call, or a high-precision wrapper around
binary64 output. The core still depends directly only on Symbolica.

See [the native evaluator](src/native/mod.rs),
[numeric primitives](src/native/primitives.rs),
[precision adapter](src/precision.rs) and [source generator](examples/generate_native.rs).
`get_expression` exposes the fully expanded exact Laurent coefficients with
explicit resource limits; see [expression inspection](README.md#combining-integrals).

## Contract

- `prec`, `PrecisionEvaluator::new` and `PrecisionEvaluator::with_backend`
  specify positive **decimal digits**.
  Python defaults to 16. The Rust binary-precision constructor instead takes
  bits and does not add guards.
- Decimal requests use a conservative decimal-to-binary conversion plus 32
  fixed guard bits. This is a fixed-precision
  calculation, not an adaptive error estimate or a promise that all requested
  digits survive ill-conditioning. Working bits are exposed by
  `PrecisionEvaluator::binary_precision()`.
- Inputs and exact rational expression coefficients are converted directly to
  the working Float precision. Python Decimal and large integers are parsed
  from exact text, including the squared renormalization scale. Ordinary Python
  floats retain their already-rounded binary64 values.
- Raw `NativeEvaluator<f64>` uses 53 bits; `NativeEvaluator<DoubleFloat>` uses
  106 bits. `NativeEvaluator<Float>` defaults to 128 bits and accepts explicit
  positive working bits. Fixed-type constructors reject a different bit count.
  DoubleFloat elementary arithmetic retains its compensated components; only
  Li2 converts both components directly to Float at 160 bits and back. It never
  narrows the whole computation to f64.
- Output precision is not artificially increased after evaluation. Python
  returns DecimalComplex components with at most the requested decimal digits
  and no more than the computed Float precision supports. Python Decimal
  construction from text does not round to the ambient decimal context.
- Real external squared invariants, nonpositive imaginary squared masses,
  positive real `mu_squared`, and finite input components remain required.
  Domain checks use the supplied high-precision values, not f64 approximations.
- Rust supports scalar and flat row-major batches; Python supports scalar and
  list-of-rows batches. Each row returns finite, simple-pole, and double-pole
  coefficients. Arbitrary-precision batches evaluate scalar rows through the
  selected Native or Expression backend, not SIMD. Shape/domain validation is
  completed before any batch output is written.
- Public master Symbols always use direct native callbacks for `Complex<f64>`,
  `Complex<DoubleFloat>` and `Complex<Float>`, independently of the manual
  default. Their Float caches are keyed by the maximum component bit precision
  supplied to the callback, including `mu_squared`. Transparent function-map compilation
  remains available and avoids making the master opaque to an outer optimizer.

`EvaluationBackend::{Native, SymJit, Expression}` selects the manual Rust route.
`DEFAULT_BACKEND` is Native, selected from measured mixed-family throughput.
`PrecisionEvaluator::with_backend(..., EvaluationBackend::Native)` selects direct
native Float; `Expression` selects the interpreter, and `SymJit` is rejected for
Float. Without an explicit selection, `PrecisionEvaluator` defaults
to Native. `evaluate_with_backend` and `evaluate_batch_with_backend` provide
equivalent explicit choices for binary64 convenience calls.

Python `backend="auto"` likewise uses Native for ordinary
float/complex/small-integer calls at `prec=16`, returning Python complex numbers.
Decimal inputs, large integers, or any nondefault precision select the Float
native implementation and DecimalComplex outputs. `backend="native"` explicitly
selects direct native evaluation in either domain; `"expression"` selects the
interpreter, and `"symjit"` rejects AP inputs. In particular, supplying a Decimal
at the default precision preserves
its value until conversion to the requested working precision; it does not
request unlimited precision implicitly.

## Dependency scope

The local Symbolica patch includes precision-scaled complex dilogarithm
convergence and exact endpoint predicates. The previous generic path limited
its convergence threshold to roughly 900 bits and narrowed some checks to f64;
merely exposing a larger `prec` would not have fixed either problem. The new
path uses arbitrary-precision stopping comparisons. The separate fast binary64
dilogarithm specialization is retained.

A separate precision-provenance audit found problems when complex components
have different relative precisions, as can happen after cancellation. The patch
also supplies Float-only Cartesian complex roots, dominant-axis phase evaluation,
precision-safe exact power identities, and direct half-integer root powers used
by the scalar expressions. Adding or subtracting zero no longer inflates the
nonzero operand's precision. Generic domains keep their existing fallback;
computed small corrections retain their uncertainty, rather than being padded
to a larger requested precision. Independent MPC root comparisons, component
precision matrices, signed-axis tests and exact zero-arithmetic checks exercise
these changes separately from integral parity.

These repairs concern the arithmetic used by these OneLOop expressions, not a
claim that every unrelated Symbolica special function has been validated at
1000 digits. The patch is bundled in `patches/symbolica-dev.patch`; it has not
been accepted upstream. Existing f64 portable caches remain f64 caches and do
not serialize a different precision's evaluator.
The independent [precision audit](ARBITRARY_PRECISION_AUDIT.md) records concrete
before/after probes and every power operation reached by the generated masters.

## Usage

For a State/JIT-free Float workspace with exactly the supplied working bits:

```rust,ignore
use oneloop::{NativeEvaluator, ScalarIntegral};
use symbolica::domains::float::Float;

let mut evaluator = NativeEvaluator::<Float>::with_binary_precision(
    ScalarIntegral::C0, 3456)?;
// Supply Complex<Float> arguments directly; evaluate / evaluate_batch use
// finite, simple-pole, double-pole output order, with mu_squared last in each row.
```

See the decimal-precision [Rust example](README.md#stack-and-precision) and the
[Python examples](python/README.md#arbitrary-precision), including Decimal complex
masses, scales, batches, and per-call precision/backend changes. Symbol access,
`initialize()`, Python import, and higher-level convenience construction eagerly
prepare all five binary64 Native and all five SymJIT backends. Float workspaces
are precision-specific and built on request. That eager symbolic work requires
an adequate calling-thread stack and the installed Symbolica license's thread
limits. Raw `NativeEvaluator` constructs no State, Atom, interpreter or JIT; it
does not remove the dependency's licensing terms or spawn a hidden worker.

The prior performance survey measures binary64 evaluation only. Its observed
slowdowns are not upper bounds for arbitrary precision. Likewise, the retained
Fortran reference tables are double-precision references: matching them cannot
establish 1000-digit accuracy. High-precision identity and convergence tests
must be considered separately from sampled Fortran parity.

## Precision checks (2026-09-08)

The final standalone Python suite passes **27 tests** (83.755 seconds), including
the Native default, explicit Native/SymJIT/Expression selection and rebuilds,
all five families at 32/1000 decimal digits, Decimal elementary vacuum controls,
context28 independence, scales from `1e-1000`, a 1050-digit tiny width,
1024-row machine/AP32 batches, exact large integers, invalid inputs and
helper-free expression inspection. Community integration passes a compile-only
check; runtime inside the community host is not established.

The rebuilt state-free native primitive gate passes **11 tests** against real
Symbolica, covering all three numeric types, signed axes, subnormal widths,
extreme magnitudes, compensated low components, unbalanced Float component
precision, and 1000-digit Li2/root controls. Three independent simple-sector
tests additionally cover 468 parameter combinations for binary64 and Float,
including signs, mass permutations, scales, and exact off-sector rejection.

Separate dependency gates include five precision-provenance regressions at
512/3456 bits, four independent MPC root regressions (including a 324-case
component-precision/magnitude matrix), sixteen existing Numerica Float tests,
and exact zero-arithmetic tests. Four state-free Symbolica Li2 checks at 3456
bits cover Catalan/reflection/inversion identities, refinement, tiny components,
both cut lips, exact predicates, checked convergence limits and reflected-series
coefficients. The binary64 Li2 grid retains a maximum scaled error near
`4.15e-16`. See [the independent audit](ARBITRARY_PRECISION_AUDIT.md) for the
distinction between direct-source, dependency-only and integration evidence.

The final rebuilt all-target Rust suite passes **114 tests**, including all five
complex native masters at 3456 bits against 1120-decimal-digit Expression
evaluation, independently per coefficient component to better than `1e-1000`.
The [release audit](NATIVE_RELEASE_AUDIT.md) records the complete test inventory,
six ignored entries and separately failing binary64 stress test. Cross-backend agreement of the
same formulas is not an independent integral oracle or a guarantee for every
kinematic point.

The repaired high-precision Li2 uses a reflected-coefficient recurrence and
MPFR's integer-zeta specialization. This avoids both premature stopping and the
very expensive general negative-zeta coefficients.
