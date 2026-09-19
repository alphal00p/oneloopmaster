# Upstream evaluator fixes (2026-09-19)

All three build roots select unmodified Symbolica and bundled Numerica at
`821b02451256a92039a0665006628bd5d91470cc`. This main revision incorporates all
seven fixes from the [regression report](patches/symbolica-main-regressions.md):
constant lifting indices, exported constant definitions, exact JIT constants,
JIT helper registration, JIT settings restoration, initialization reentry and
Git metadata watches. The original patch files and tarball are historical;
none is applied to a dependency.

The cache discriminator and Python revision metadata follow the new pin.
OneLOop retains its source-evaluator cache format and prepared helper scopes;
the latter still avoid the separate non-inlined compilation-performance issue.
The explicit-pi wrappers are no longer needed for upstream callback-index
correctness, but remain part of that shared-scope preparation. SymJIT 2.25.6
and its numerical callbacks/settings remain as in the previous port.

Validation of this update:

- All 11 focused regressions from the patch bundle pass on upstream main,
  along with all 20 existing Symbolica evaluator tests. Only test files were
  copied into the isolated checkout; its implementation is unmodified.
- All five caches were regenerated, serialized, restored and checked against
  all 324 fixture rows (three Laurent coefficients each). Their documented
  sizes and SHA-256 checksums match the new assets.
- All 33 focused release Rust tests pass (one subprocess helper is ignored):
  cache reloads and mixed batches, nested/tagged definitions, bare master JIT
  callbacks, native primitives, cold startup, arbitrary precision and native
  code generation. The cold-start test now directly calls transcendental
  accessors from a dependent initializer to exercise the upstream reentry fix.
- All ten shared-host Python checks pass: elementary calls, benchmark fixtures,
  expression interop and all seven decimal/arbitrary-precision tests, including
  all five scalar families at 32 and 1000 digits. The shared host used a debug
  build. Its first interop run found missing Nix NumPy runtime libraries;
  setting the GCC/zlib library paths fixed the environment and the rerun passed.
- Both optimized standalone Python startup tests pass, including fresh-process
  calls/evaluators across all five families, heterogeneous batches, at most an
  8 MiB main-thread stack, and environment/file configuration guards.
- Core all-target and standalone Python Cargo checks pass. Clippy with warnings
  denied, all three build roots' formatting checks and `git diff --check` pass.

Tests used Linux x86-64, Rust 1.98.1 and CPython 3.12; Clippy used Rust 1.97.1.
The standalone extension and Rust tests used release builds. SymJIT is unmodified
2.25.6. No dependency source patches are applied.

The reports below preserve validation of the earlier revisions.

---

# Upstream radial precision update (2026-09-19)

The three build roots now pin unmodified Symbolica and bundled Numerica at
`08defb398420ab07a428876d4fd10eaf6016470f`. Relative to the previous port, upstream
changes only `Real::hypot` and `Complex::log`: exact constants use the stronger
component's precision, including the small real logarithm near the unit circle.

OneLOop's temporary `__precision_log` and `__precision_abs` registrations and
lowering have been removed. Interpreted expressions use upstream operations
directly. The existing mixed-precision regression remains, with an additional
near-unit regression at 128, 512 and 3456 bits. Native code generation no longer
needs aliases for the removed callbacks. The cache format and Python revision
metadata follow the new pin; the upstream patch artifacts retain their original
bases as historical references.

Validation of this update:

- All eight independent MPFR/MPC radial precision regressions pass against
  unmodified upstream Numerica.
- All five caches were regenerated, serialized, restored and checked against
  all 324 fixture rows (three Laurent coefficients each).
- `cargo clippy --all-targets -- -D warnings` and formatting checks pass.
- All 34 targeted Rust tests pass (one existing subprocess helper is ignored):
  box limits, portable evaluators, function registration and mixed precision,
  native primitives/startup, arbitrary precision, and native code generation.
  This includes the formerly failing box comparisons and the new near-unit
  logarithm regression at 128, 512 and 3456 bits.
- Ten shared-host Python tests pass: elementary calls, all benchmark fixtures,
  expression interop, and all seven decimal/arbitrary-precision tests, including
  all five scalar families at 32 and 1000 digits.
- Both standalone Python startup tests pass: fresh main-thread direct/evaluator
  calls across all five families and heterogeneous batches, with at most an
  8 MiB stack, plus environment/file configuration guards.

The completed validation below records the previous `3cb73a6` port.

---

# Migration to Symbolica main (2026-09-19)

Target: `3cb73a6d958755eaa86088ff713dd13e6bda217b`, fetched after the force push.
All three Cargo build roots select unmodified upstream Symbolica and its bundled
Numerica at this revision. SymJIT is the unmodified 2.25.6 release.

The source migration uses public interfaces:

- Numerical dilogarithms resolve the registered `polylog(2,z)` callback once per
  numerical domain; DoubleFloat conversion keeps both compensated components.
- OneLOop retains the original exact definitions for inspection and registers
  equivalent wrappers passing pi explicitly to nested functions. This avoids
  Symbolica's lifted-constant/external-function indexing defect.
- Contexts share the prepared definitions. Closed helpers use common internal
  argument slots to avoid recursive capture-pruning rebuilds. Public signatures
  and inspection retain the original arguments. Evaluator construction disables
  the default whole-map Horner pass to retain the original arithmetic structure.
- Portable caches serialize the numerical `ExpressionEvaluator` through bincode,
  validate its dimensions, and recompile every level with explicit JIT settings.
- Native code generation lifts symbolic constants as parameters before using the
  public instruction export, and checks for unresolved constant placeholders.
- JIT construction lowers complex square roots and conditions to registered
  numeric callbacks; exact/interpreted expressions retain their original
  operations. This avoids wrong root branches and imaginary conditions being
  treated as zero by upstream SymJIT.
- Assumption queries and floating sign-bit checks use the current APIs.
- Interpreted logarithms and magnitudes use registered callbacks that build the exact identity
  in the scaled magnitude at the stronger component's precision. Upstream
  `hypot` constructs it from a possibly two-bit ratio, degrading otherwise
  well-determined logarithms. Computed components retain their tracked precision;
  inspection and generated native arithmetic keep their existing formulas.
  A 128-bit `2 + i*2^-120` with a two-bit imaginary component reproduces a
  five-bit real logarithm upstream; the adapter retains at least 120 bits.

The separate [upstream Numerica patch](patches/numerica-hypot-precision.md)
repairs the radial precision defect in the dependency itself. It is supplied in
standalone and Symbolica-tree formats and is not applied by OneLOop.

Validation:

- `cargo check --all-targets` and the shared Python host compile.
- `cargo clippy --all-targets -- -D warnings` passes (Rust/Clippy 1.97.1).
- The complete release unit/integration suite passes: 121 tests, six ignored.
- Eight existing Numerica zero/sqrt/phase/power regressions pass on unmodified main.
- Eleven native primitive regressions pass, including thousand-digit polylog.
- Five focused registration tests pass for tagged calls, inspection, captures,
  mixed-precision magnitudes/logarithms, and scalar/batched complex branches
  before and after serialization.
- Both box-limit tests pass, including the nine initially failing 128-bit
  comparisons, after the logarithm and magnitude adapters were added.
- Eight native code-generator tests pass.
- The Rust documentation test passes.
- The shared Python host passes 27 tests; two standalone-only tests are skipped.
- Both standalone-only tests pass in the standalone extension, covering fresh
  main-thread startup through direct and evaluator calls, all five scalar
  families, heterogeneous batches, and environment/file configuration guards.
  The cold processes use at most an 8 MiB main-thread stack.
- Both Numerica patch formats fail four new regressions before the fix and pass
  all eight afterward. Full standalone/bundled suites pass 194/196 tests; the
  standalone patch combined with the earlier precision fixes passes 207 tests.

All five embedded caches have been regenerated. Their numerical evaluator
payloads pass all 324 fixture rows after serialization and reload.

Direct translation omits fractional-power call targets, so ordinary compilation
is selected. SIMD crashes on zero-input batches and is disabled. Packed complex
arithmetic with fastmath disabled avoids spurious imaginary residues in
conjugate products; the generic complex compiler failed D0 `line_233`. Registered
callbacks preserve complex root branches and predicates. Helper aliases preserve
the original relative symbol order, with a cancellation regression covering it.

The earlier release-only investigation below is historical; its statement that
the dependency switch was reverted describes that earlier attempt.

---

# Symbolica 3.0 migration compatibility check (2026-09-19)

The patch-free migration is not complete. Published Symbolica 3.0.0 does not
provide all the APIs and numerical repairs used by the current OneLOop code.
The attempted dependency switch was reverted; the existing manifests, patches,
Python lockfiles and evaluator assets remain unchanged.

The crates.io release resolves to Symbolica 3.0.0, Numerica 3.0.0, Graphica
3.0.0 and SymJIT 2.25.6. Symbolica's release archive records upstream commit
`1f38ba9e0de48cb9f9e8c616b2ccd729fa45f606`.

## Compilation findings

Temporarily replacing all three Symbolica dependencies with `version = "=3.0.0"`
and removing every `[patch]` section allows dependency resolution, but
`cargo check --all-targets` fails with 18 library errors:

- `dilog_complex_f64` and `dilog_complex_float` are absent from the public
  transcendental module. The generated native evaluator calls these kernels
  without initializing Symbolica state, including for arbitrary precision.
  The public `polylog(2, z)` evaluator can replace this API access, with changes
  to the native evaluator's state-free behavior. The absence of these particular
  function names does not by itself prevent migration.
- `JITCompiledEvaluator` has no `export_portable`, `import_portable`,
  `input_count` or `output_count` methods. The release does implement bincode
  serialization, so a replacement should use that interface and validate
  nested callbacks, strict compilation settings and dimension checks.
- `FunctionMap::get_definition` is absent. The expression inspection code uses
  it to expand tagged definitions with lexical parameter bindings.
- Assumption queries now return `ConditionResult`, `evaluate_with_prec` needs
  an explicit hash-map type, and `Real::hypot` conflicts with OneLOop's method.
  These are ordinary source migration changes, distinct from missing APIs.

Source inspection also finds that `ExportedInstructions` has no
`constant_functions` field, which `examples/generate_native.rs` needs to
distinguish registered constants from rational placeholders. The library
failure prevents Cargo from reaching this example during the attempted build.

## Numerical compatibility

The release's integer polylogarithm series still computes its stopping threshold
as `2f64.powi(-(binary_prec.min(900) as i32))` and converts term magnitudes to
binary64. This is the path addressed by the former arbitrary-precision patch;
a successful API migration alone would not establish the existing precision
contract.

Eight existing dependency-only regression tests were initially copied unchanged
into an isolated crate: one passed and seven failed. A follow-up corrected sign
assertions for Numerica 3.0's API: `is_negative()` now means strictly below zero;
`is_sign_negative()` tests the sign bit, including negative zero. With those
assertions migrated, **two pass and six fail** on the unmodified release.

| Test file | Result |
| --- | --- |
| `tests/complex_sqrt_precision.rs` | 1 failed: nearly imaginary roots lose required precision. |
| `tests/float_zero_precision.rs` | 2 passed after migrating sign-bit assertions. No zero-arithmetic repair is required. |
| `tests/complex_precision_provenance.rs` | 5 failed: mixed-precision phase, integer powers, half powers, power branch lips, and signed-axis phase. |

These tests exercise Numerica arithmetic directly through Symbolica's public
types. They require neither OneLOop startup nor evaluator caches, so the failures
are independent of the missing OneLOop integration APIs. They establish specific
regressions against the current contract, not exhaustive dependency coverage.

The follow-up [Numerica precision patch](patches/numerica-3.0.0-precision.md)
addresses the reproduced Numerica defects and includes the migrated regression
tests. It does not resolve the separate Symbolica API and polylogarithm blockers.
The subsequent [Symbolica polylog patch](patches/symbolica-3.0.0-polylog-precision.md)
addresses the polylog precision defects, including the public expression route,
and passes its regressions with unmodified released Numerica. Neither patch
artifact has been applied to OneLOop's dependencies.

## Reproduction

The compile output is saved locally at `/tmp/oneloop-3-check.log`, and the
reverted trial diff at `/tmp/oneloop-symbolica-3-attempt.diff`.

An isolated dependency-only crate at `/tmp/oneloop-symbolica-3-probe` uses
published Symbolica 3.0.0 with the same explicit features as OneLOop and copies
the existing numerical regression tests unchanged. It does not depend on the
OneLOop library or its local patches. Its output is saved at
`/tmp/oneloop-3-probe.log` and `/tmp/oneloop-3-precision.log`.

Completing the migration requires replacements for the missing interfaces,
validation of the numerical repairs against the released dependencies, and
regeneration and validation of the portable evaluator assets. Moving the old
patches into a vendored dependency would not meet the requested patch-free setup.
