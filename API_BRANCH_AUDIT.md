# Expression API and branch-selection audit

Audit date: 2026-09-09. Scope: public Rust/Python expression interfaces, native
Symbolica assumptions and replacements, shared-kernel Python packaging, and the
new analytic-region selection helper. This is a focused API audit, not a claim
of exhaustive numerical parity over every integral region.

## Required contract

- `get_expression(master)` accepts an actual Symbolica master-call expression.
  All kinematic invariants, squared masses, and the final squared renormalization
  scale may themselves be composite expressions. It returns the three Laurent
  coefficients in `(finite, simple_pole, double_pole)` order.
- Inspection calls omit the Laurent tag. The existing evaluate-able coefficient
  symbols use a leading tag `0`, `-1`, or `-2`; the distinction must be explicit
  in examples and argument errors. The `oneloop::` inspection shorthand must not
  be advertised as an independently registered numerical evaluation hook.
- Symbolic inputs are never parsed from strings by OneLOop. Users construct
  expressions with Symbolica's `E`, `S`, and `N`. Existing Symbolica attributes
  are retained and queried through its own assumption logic. There are no
  OneLOop `real=[...]` or `positive=[...]` declarations.
- `select_branch(expression, replacement_rules)` uses Symbolica `Replacement`
  objects. Rules are applied only to a temporary condition used to choose an
  `if` arm. Chosen formulas retain their parameters; an undecidable condition
  also remains parametric. Nested `if`s are processed in conditions, retained
  arms, ordinary function arguments, sums/products, and power bases/exponents.
- Native `if(c, yes, no)` means **nonzero versus zero**, including complex `c`;
  it is not a direct positivity test. A sample selects a region, not a formula
  valid in every region. Unresolved or undefined numerical predicates must not
  be silently treated as proof of a branch.
- This symbolic-interface change does not revoke the separate numerical API:
  numerical scalars and batches, `Decimal`/`DecimalComplex`, precision settings,
  and backend configuration remain supported. A backend name is configuration,
  not a symbolic expression supplied as text.

## Source findings

The original standalone Python inspection path parsed `Vec<String>` and added
its own assumption declarations. The community path already exchanged genuine
expressions but still selected master families with strings. These are the
interfaces requiring replacement. Rust's scalar/tensor constructors already
accept `Atom` inputs, including composite masses and scales; its enum-based
family API does not have the string-parsing defect. `compile_native` and the
Rust evaluator builders already accept expression-valued parameters and should
not unnecessarily narrow them to symbols.

The pinned Symbolica source provides the necessary native boundaries:

- `src/api/python/expression.rs`: `ConvertibleToExpression` explicitly rejects
  strings, accepts native expressions and supported numeric conversions, and
  `FromPyObject<Symbol>` accepts an expression exactly when it is a variable.
  OneLOop does not need a competing expression-to-symbol conversion policy.
- The same file exposes `PythonExpression::replace_multiple`. The inner Rust
  field of `PythonReplacement` is private, so a Python adapter must use the
  public native replacement method instead of reconstructing rules and losing
  matching conditions, wildcards, or rule settings.
- `src/id.rs`: native `is_real`/`is_positive` already propagate assumptions
  through supported composite expressions. The old local recursive
  `known_real` inference duplicated that logic and has been reduced to the
  native `is_real()` query. Symbolica's positivity includes zero; it cannot
  independently prove a nonzero `if` predicate.
- `src/normalize.rs`: native `if` normalization resolves numeric conditions by
  testing zero. The branch helper must check that sampled conditions are
  defined/finite before relying on numeric truth, rather than assuming every
  numeric atom is a valid analytic-region decision.
- `src/api/python.rs`: `SymbolicaCommunityModule` is the supported mechanism for
  sharing Symbolica state. Registering Python classes/functions is separate
  from symbol/evaluator initialization. A separately linked standalone kernel
  must not accept host-kernel expression objects by pretending they share state.

These paths refer to the local pinned Symbolica dependency described in
`Cargo.toml`, not to an independently installed Python wheel version.

## Verification checklist

Completed source review and runtime results are recorded below. These checks
cover the API change; they do not certify every possible symbolic identity or
analytic region of the underlying integrals.

- [x] Rust and Python expression-only signatures; strings and separate
  assumption lists rejected with clear errors.
- [x] Composite master arguments and native namespaced attributes retained;
  invalid master names, wrong arity, and misplaced Laurent tags rejected.
- [x] Integer, rational, floating-point, and high-precision replacement values.
  No conversion of arbitrary-precision condition inputs through binary64.
- [x] Every surviving nested `if` is visited, including deeply nested
  conditions and conditionals inside powers and generic functions.
- [x] Partial/unknown replacements preserve the original condition and all
  parameters in value branches; no-`if` expressions are unchanged.
- [x] Native replacement order and wildcard/conditional-rule semantics survive
  the Python adapter. Undefined discarded branches are not traversed.
- [x] Equal-mass parametric B0 samples remove every `if` using integer and
  floating-point sample points, preserve `psq`/`m2`, and numerically agree with
  the complete expression in the selected analytic region.
- [x] Exact zero, threshold, below-/above-threshold, spacelike, and complex
  nonzero condition cases; unmatched/partial rules leave conditions unresolved.
- [x] Single-kernel Python import and native `Expression`/`Replacement` identity;
  manual numeric/Decimal/batched evaluation regression tests remain passing.

## Implementation cross-audit

The final iterative traversal has been inspected independently: it visits
conditions before arms, skips discarded arms, preserves the original condition
when a decision is unavailable, and rebuilds ordinary function arguments and
power positions. The Python adapter obtains native replacements through
Symbolica's public `PythonTransformer::replace_multiple` and exposed transformer
chain. It then calls the same Rust helper, retaining native rule settings instead
of rebuilding Python replacements. No private dependency fields or second
Symbolica kernel are used.

An independent precision probe found and verified the repair of a concrete issue:
native Symbolica replacement in `sqrt(1+x)-1`, with
`x = Float::parse("1e-1000", Some(64))` (also at 256 bits), normalizes the sampled
condition to numeric zero. The same sample represented at 4096 bits gives the
expected nonzero value near `5e-1001`. Evaluating the already-normalized zero at
higher precision cannot recover the lost term. This occurs before the branch
helper's predicate evaluator receives its input; merely avoiding binary64
conversion is not sufficient to prevent it. This result was reproduced against
the pinned local Symbolica build using native `replace_multiple` followed by
4096-bit evaluation.

The repair converts finite floating-point replacement values to
their exact represented rational values before native substitution. Static
review confirms this preserves native patterns/settings and protects the
returned parametric formulas. An independent runtime probe also identified a
second obstacle in the initial repair: Symbolica's exact radical normalization
unconditionally factorizes large residual integers. That repaired probe entered
that path for a roughly 3322-bit numerator and used a full CPU for over 80
seconds without finishing; the audit stopped only its own temporary probe.
The relevant path is `simplify_perfect_power_factors` in Symbolica's
`src/coefficient.rs`.

Both findings are resolved without modifying Symbolica: sample replacements
temporarily wrap the exact value `r` in its native `if(pi, r, r)`. This identity
defers arithmetic to the numeric evaluator, avoiding both premature floating
rounding and symbolic integer factorization. It is used only in the private
sample and cannot appear in returned parametric formulas. The helper retains
coefficient precision when choosing arbitrary-precision working bits.

An independent executable compiled the final `branch_selection.rs` directly
against the pinned Symbolica library. It passed the `sqrt(1+x)-1` reproduction
with 64-, 256-, and 4096-bit floating replacement values in under 0.1 seconds,
retaining the correct symbolic arm. Additional independent assertions passed
for partial replacements preserving the original condition and a native
`ReplaceWith::Map` callback returning the tiny floating value.

The implementation lane also ran all eight final branch-selection Rust tests in
an independently compiled harness: all passed in 5.26 seconds. Coverage included
4096 nested value branches, 1024 nested conditions, callback/transformer rules,
undefined-condition handling, and equal-mass B0 across seven samples and nearby
points within the same regions, including exact structural equality to the
below-threshold root/log formula.

The final repository-wide release Cargo run passed 124 tests, with zero failures
and six pre-existing ignored tests across 41 test binaries. It includes the
eight branch-selection tests and both expression-first master-inspection tests.
The new Rust documentation example also passed. Release Clippy passed with
`-D warnings` for the core's full target set and both Python build modes; the
patched SymJIT dependency still emits its existing uncommon-codepoints warning.
Reproduce the Rust checks with:

```sh
cargo test --release --offline --all-targets --no-fail-fast -- --test-threads=1
cargo test --release --offline --doc -- --test-threads=1
cargo clippy --release --offline --all-targets -- -D warnings
```

The final shared-host suite was inspected: 29 tests ran in 56.759 seconds,
with 27 passing and two standalone-only checks explicitly skipped. Of those
27 executed tests, 22 use the actual Symbolica host and five test the performance
harness with fake evaluators. All five new inspection tests passed, including:

- Actual `Expression` outputs and rejection of string/assumption-list inputs.
- Equal-mass B0 with samples `3.23`, `3`, `1/2`, `-2`, `5`, `0`, and `4` at
  squared mass/scale one: no remaining `if`s, symbolic parameters retained, and
  agreement with direct native B0 evaluation.
- Native attributes on composite masses, native family-symbol extraction, and
  expression-valued numeric input `1/3` at 64 decimal digits.
- Nested conditions, 1200 nested value branches, function/power positions,
  partial replacements, and native wildcard/conditional-rule matching.
- The tiny floating-point cancellation regression described above.

The same final host run covered genuine shared-kernel expression compilation,
native/SymJIT fixture agreement, 1024-row machine/Decimal batches, and all five
scalar families at 32/1000 decimal digits. Reproduce it using
`python/host/test.sh` in the documented isolated environment. The standalone
extension's cold-load and environment-guard checks are intentionally not counted
as host checks. These results do not establish global analytic parity or a
certified numerical zero-decision procedure for arbitrary transcendental
identities and singular boundaries.

The rebuilt standalone numeric Python extension also passed its final suite:
29 tests ran in 66.447 seconds, with 24 passing and the five shared-kernel-only
inspection tests explicitly skipped. This separately verifies the preserved
numeric API, eager startup, batching, arbitrary precision, and the rejection of
symbolic entry points in a separately linked kernel.

The captured equal-mass demonstration is stored in
`examples/expressions/b0_equal_mass.txt`. At `psq=3.23, m2=1`, the finite
coefficient is a 758-character parametric root/log expression, including fully
qualified symbol names. It contains no `if`, `abs`, or `conj`; the pole
coefficients are exactly `(1, 0)`. Root aliases used in the README only improve
the mathematical display and are not opaque functions in the returned atom.
