# Upstream status

All seven fixes below are upstream in Symbolica main
`821b02451256a92039a0665006628bd5d91470cc`, selected by OneLOop and both Python
build roots. The patch files and tarball are historical artifacts against
`08defb39`; do not apply them to current main.

| Fix | Upstream commit |
| --- | --- |
| Git metadata watches | `cc5cdaa` |
| Callback indices after constant lifting | `65ff188` |
| Exported constant definitions | `35d424a` |
| Direct JIT constant resolution | `1699ac5` |
| Helpers in the Symbolica namespace | `b6736cb` |
| JIT compilation settings after restoration | `78ee71d` |
| Transcendental initialization reentry | `821b024` |

All 11 focused regressions from this bundle and 20 existing evaluator tests
pass on the unmodified upstream implementation. The standalone tests were
copied into an isolated checkout without applying the source patches.

Current-port validation is recorded in the
[migration notes](../SYMBOLICA_3_MIGRATION.md). The report below records the
original reproductions and proposed patches.

---

# Remaining Symbolica regressions

Verified against unmodified main
`08defb398420ab07a428876d4fd10eaf6016470f` on 2026-09-19. This report excludes
non-inlined-function compilation performance and the SymJIT defects already
submitted to its maintainer. The JIT issues below are in Symbolica's
constant resolution, function registration and serialization. No patch changes
SymJIT code.

These patches are upstream proposals. None is applied to OneLOop's dependencies.
The earlier hypot, complex-logarithm, square-root, phase, power and polylogarithm
precision repairs are already in the selected Symbolica revision.

| Issue | Observable failure | Patch |
| --- | --- | --- |
| Lifted constants leave stale callback indices | Exporting/evaluating a nested non-inlined function can panic or select the wrong callback. | [Callback indices](symbolica-main-constant-function-indices.patch) |
| Exact instruction exports lose constant definitions | `pi` and `polylog(2,1/4)` both export an anonymous rational zero instead of enough information to recover their values. | [Exported constants](symbolica-main-exported-constants.patch) |
| Exact evaluators JIT-compile unresolved constants as zero | Direct JIT compilation of `pi` returns `0`; registered complex constants also become zero. | [JIT constants](symbolica-main-jit-constants.patch) |
| User functions in the Symbolica namespace are skipped | Non-inlined helpers fail complex/vector JIT compilation with an unsupported-operation error. | [JIT namespace](symbolica-main-jit-namespace.patch) |
| Cold transcendental accessors deadlock | A first `polylog()`, `tan()` or `bessel_j()` call hangs when a dependent state initializer calls the accessors. | [Startup reentry](symbolica-main-transcendental-startup.patch) |
| JIT restoration drops child compilation settings | A round-trip silently recompiles non-inlined children with default options and can change their numerical results. | [JIT settings](symbolica-main-jit-settings.patch) |
| Git metadata watches are incorrect | An unchanged linked worktree rebuilds repeatedly; branch commits can leave the version stale in ordinary clones. | [Build metadata](symbolica-main-build-metadata.patch) |

## Callback indices after constant lifting

The compiler hoists a sub-evaluator's constants into the parent, removes their
entries from `external_fns`, and leaves runtime `ExternalFun` instructions with
indices into the old table. The test on current main panics at instruction export:
`index out of bounds: the len is 2 but the index is 3`.

The existing patch remaps the surviving callback indices before filtering.
Its tagged, nested-function regression covers mixed tables, constant-only tables,
ordinary callback tables, numerical evaluation, bincode restoration and JIT
compilation. The patch originally targeted `3cb73a6`; it applies unchanged and
has been revalidated at `08defb39`. This is a correctness fix, independent of
the excluded non-inlined-function performance problem.

## Exact exported constants

On the unmodified base, both expressions below export the same instruction and
constant data:

```text
pi              -> Assign(Out(0), Const(0)), constants=[0+0i]
polylog(2,1/4)   -> Assign(Out(0), Const(0)), constants=[0+0i]
```

Their numerical evaluator values are approximately `3.141592653589793` and
`0.2676526390827326`. The zeros are internal placeholders. The exported structure
does not identify them as placeholders or retain the registered function, tags,
or fixed arguments, despite promising the data needed for external evaluation.

The patch adds `ExportedConstantFunction` and
`ExportedInstructions::constant_functions`. Each entry identifies the constant
slot, registered symbol, canonical tags and exact fixed arguments. Consumers can
resolve those definitions at their own target precision. The ordinary constants
array and evaluator serialization stay unchanged.

Two tests reconstruct the exported constant values at 53, 128 and 512 bits,
covering pi, fixed polylogarithms, differently tagged constants, coefficient
conversion, bincode round-trips and constants lifted from a non-inlined function.
This is an additive public-API change: downstream code constructing an
`ExportedInstructions` literal must supply the new field.

## Constants in direct JIT compilation

`ExpressionEvaluator<Complex<Rational>>::jit_compile` sends the exact evaluator's
constant placeholders directly to SymJIT. It resolves runtime callbacks, but
skips the values of entries marked with `constant_index`. Direct compilation
therefore returns `0` for `pi` and `0+0i` for a registered `2+3i` constant.
Mapping coefficients to a numerical domain first avoids the defect.

The patch resolves registered constants and fixed-argument function calls at
the target domain's precision before compiling. It propagates errors from
constant callbacks and domain conversion through the existing `Result` API.
Unlike the export-metadata patch, this repairs Symbolica's own JIT path and
requires no public API or serialized-format change.

Four tests cover pi, fixed polylogarithms, tagged constants, rational mixtures,
complex constants, rejection of non-real constants in a real domain, callback
errors, and constants lifted from non-inlined helpers with input parameters.
The latter exercise scalar real, scalar complex and four-lane real evaluation,
with direct translation both enabled and disabled.

## User helpers in the Symbolica namespace

Three JIT registration paths skip any symbol classified as built-in without
registered evaluation metadata. That classification includes all names in the
`symbolica` namespace, even user-defined functions that have a compiled
sub-evaluator. A nested helper therefore fails complex compilation with:

```text
op_code cplx_symbolica_jit_namespace_g is not found or is not supported
```

The patch skips such built-ins only when they have no sub-evaluator. Its
regression uses two levels of non-inlined helpers in the `symbolica` namespace
and checks real, complex and four-lane real results with direct translation
both enabled and disabled. This is a function-registration correctness defect,
separate from the excluded compilation-performance issue.

## Transcendental initialization

The lazy symbol caches acquire their own initialization lock before entering
Symbolica State. State then runs dependent initializers, which may call the same
accessors and wait forever on the pending cache lock.

All three fresh-process baseline cases timed out: special functions (`polylog`),
geometric functions (`tan`), and Bessel functions (`bessel_j`). The patch enters
State before locking a cold symbol cache. Warm accesses use the cached value
directly. Its subprocess test checks all three entry routes and successful
completion of the dependent initializer.

## JIT compilation settings during restoration

Symbolica serializes the external-function definitions and top-level SymJIT IR,
but `JITCompiledEvaluator::load` reconstructs non-inlined children with
`JITCompilationSettings::default()`. Explicit options such as `fast_complex`,
`fastmath`, SIMD policy, direct translation and optimization level are lost.

The baseline probe compiles `f(z)=z*conj(z)` with `fast_complex=false`,
`fastmath=false`, SIMD disabled and O2. At `z=0.1+0.2i`, a bincode round-trip
changes its imaginary result from `1.6653345369377347e-18` to zero because the
child switches backends. The discrepancy between those backends is an already
reported SymJIT issue; **Symbolica's separate defect is silently losing the
requested configuration**. The regression checks preservation of results and
policy, not the absolute accuracy of either backend.

The patch retains the requested settings in all four JIT numeric domains and
passes them recursively during restoration. Its regression covers two levels
of non-inlined functions, direct/indirect translation, both complex backends,
scalar and batch evaluation, and two successive bincode round-trips.

**Compatibility:** this patch adds settings to the serde/bincode JIT payload.
Previously serialized `JITCompiledEvaluator` payloads must be regenerated.
`ExpressionEvaluator` serialization and OneLOop's current source-evaluator cache
format are unaffected.

## Build metadata

The build script watches `.git/HEAD`. In linked worktrees `.git` is a file, so
that watched path does not exist. A minimal Cargo project using the unmodified
build script rebuilt its library on two consecutive checks with no changes.
With the patch, the second check reports both the build script and library as
fresh. In an ordinary clone, watching only `HEAD` also misses updates to the
branch reference to which it points.

The patch asks Git for the real HEAD/ref locations and watches only existing
paths. Its test covers ordinary clones, linked worktrees, loose and packed refs,
and source directories without Git metadata.

## Applying and testing

All seven patches apply to `08defb398420ab07a428876d4fd10eaf6016470f`. Apply them
from a clean Symbolica checkout, using the paths in the table above. The callback
index patch is the existing artifact; do not apply it twice. The patches may be
reviewed and applied independently; the three JIT patches modify separate parts of `src/evaluate/backend.rs`.

Run the focused tests after applying the corresponding patches:

```sh
RUST_MIN_STACK=134217728 cargo test --no-default-features \
  --features integer-gmp,float-mpfr,native_code_generation,bincode,serde \
  --test lifted_constant_indices --test exported_constants \
  --test transcendental_startup --test jit_settings_roundtrip \
  --test build_metadata --test jit_constants --test jit_namespace \
  -- --test-threads=1
```

Validation uses Rust 1.98.1 on Linux x86-64 and unmodified SymJIT 2.25.6.
Validation completed:

- All seven patches pass `git apply --check` individually on the pristine base
  and apply together. The resulting sources and tests match the tested checkout.
- All 11 focused regression tests pass. The startup harness additionally runs
  its ignored helper in three fresh subprocesses.
- All 20 existing evaluator tests and all 53 existing transcendental tests pass.
- Both serde and bincode feature paths compile; JIT restoration is tested at
  runtime through bincode.
- Rust formatting and `git diff --check` pass for the changed code.
- The linked-worktree Cargo reproduction rebuilds before the fix and stays
  fresh after it.

The existing suites were run serially with `--lib evaluate::test::` and
`--lib transcendental::tests::` respectively. The transcendental suite validated
the startup fix before the subsequent JIT-only changes. Validation covers Linux
x86-64; no claim is made about other JIT architectures.
