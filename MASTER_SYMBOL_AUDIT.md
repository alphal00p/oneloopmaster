# Public master Symbols: API audit

Audit date: 2026-09-07. **The requested public master-Symbol interface is not
implemented in this snapshot.** This is separate from the failing numerical
parity audit. No numerical master callback was added as part of this review.

## Current implementation

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

## Recommended implementation, not yet applied

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
