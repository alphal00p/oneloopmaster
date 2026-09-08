# Implementation audit — 2026-09-07

For the subsequent direct-Rust/ArbPrec implementation, see the
[current native release audit](NATIVE_RELEASE_AUDIT.md): 114 release Rust tests
and 27 standalone Python tests pass. The independent precision audit, explicit
14-failure native binary64 diagnostic, expression-inspection limits and matched
new-backend timings are recorded there. The snapshots below preserve the earlier
SymJIT-default implementation and must not be read as current backend status.

For the **current 2026-09-08 release validation**, see the dated sections of
[SCALAR_PARITY_AUDIT.md](SCALAR_PARITY_AUDIT.md) and
[MASTER_SYMBOL_AUDIT.md](MASTER_SYMBOL_AUDIT.md). All 1,026 expanded Fortran
coefficients now pass at 256 bits, and master Symbols are implemented. Full
default release-suite revalidation passes **57 tests**, with zero failures and
four ignored entries (three opt-in diagnostics and the initialization helper),
including all 879 acceptance coefficients through both master hooks and native
maps. The new bare-master outer-JIT/portable-batch regression also passes: all
293 rows in canonical tag order, plus eight controls per family in each of the
other five orders. The final expanded
fixed-f64 rerun retains 14 small-momentum B0/dB0 failures, down from 44. Separate
repeated-root box and vacuum-triangle repairs now pass their focused regressions;
the scalar audit records the remaining 128-bit tiny-width diagnostic explicitly.
Global analytic-region coverage remains unproved, and measured runtime does not
meet the requested 1.5-times-Fortran target. Publication, final dependency
configuration, licensing permission, and actual community-host validation remain
open. These results use dev `fb845d34` plus the documented dependency patches;
they do not validate Python or other CPU architectures by themselves.

The text below preserves the **2026-09-07 audit baseline**. Its old test counts,
missing master interface, and then-unrepaired counterexamples describe that
earlier stage, not the current implementation.

This is a development audit, not certification of complete Fortran equivalence
or a publication approval. The original Fortran wrapper and the independent
dependency-free numerical Rust core were not changed.

**Subsequent scalar-domain audit: full parity fails.** The wider scan found
44 disagreements among 1,026 coefficients from 342 accepted original-Fortran
points. See `SCALAR_PARITY_AUDIT.md` and its retained failing fixtures. The
passing acceptance results below remain valid but are not comprehensive.

The dedicated `oneloopmaster` checkout was revalidated after export: all 31
acceptance tests pass, the new `cargo run --example basic` runs and checks its
output, and all-target Clippy, formatting and warning-free rustdoc pass.
The explicitly enabled 342-point exploratory test was rerun in this checkout
and failed with the same 44 coefficient disagreements (exit status 101).
README now includes setup, usage, conventions and credits to Andreas van Hameren
with the original OneLOop references.

Three follow-up auditors investigated the mismatches by family and cross-checked
the proposed branch/limit corrections. See the follow-up section of
[SCALAR_PARITY_AUDIT.md](SCALAR_PARITY_AUDIT.md). Those are diagnostic findings;
no production formula was changed in this documentation/audit follow-up.
An additional Symbolica-source/API audit found that the requested public master
Symbols with dedicated `EvaluationInfo` hooks are missing. The distinction from
native function-map definitions and the proposed design/tests are recorded in
[MASTER_SYMBOL_AUDIT.md](MASTER_SYMBOL_AUDIT.md). This criterion remains open.

## Four independent auditors and cross-review

| Auditor | Primary scope | Cross-review |
| --- | --- | --- |
| Analytic | A/B, triangle/box sectors, branches, original formulas | Test inventory and infrared coverage |
| Tests | Oracle independence, precision, tolerances, false passes | Analytic regression inputs and native endpoint defect |
| Native design | Transparent function maps, evaluator patches, API and stack | Bubble lips, undefined dB0, vacuum predicates |
| Release | Reproducibility, packaging, artifacts, documentation | Analytic fixes, oracle hardening, dependency patches |

The lead agent coordinated all numerical runs serially to respect the Symbolica
runtime license. Auditors exchanged findings rather than merely issuing four
independent approvals.

## Corrections made during this audit

- Massive real-axis B0/dB0 logarithmic lips above threshold; the corresponding
  derivative endpoint logs are expressed as ratios to avoid signed-zero cuts.
- Mixed-sign massless C0 and D0 keep the separate continued logarithms instead
  of collapsing them into a principal logarithm of a ratio.
- Massive triangle infrared sectors preserve sheet arithmetic. IR2 also had a
  wrong dilogarithm argument: the reference uses Li2(1-x), not Li2(x).
- Exact equal-mass D0 vacuum limit, `1/(6(m²)²)`, and nonfinite rather than zero
  evaluation for the undefined fully scaleless derivative.
- Real-axis reciprocal roots use real square roots explicitly, avoiding tiny
  imaginary leakage from multiplying two complex negative-axis square roots.
- Native exact-zero predicates use Boolean `if` chains instead of squaring
  invariants or complex mass differences.
- Symbolica's native polylog endpoint/real-axis predicates retain the input
  precision instead of recognizing rounded f64 values as exact endpoints.
- Independent Fortran tests now cover A0/B0/dB0 and infrared poles. Inventory,
  shapes, finite values, topology tags, mass masks and pole coverage are enforced.
  Tolerances are dimensionally normalized: the former absolute floor could have
  allowed a 41.7% relative error in one small box fixture.
- Auxiliary sheet fixtures now enforce all 324 rows. A clean-checkout bootstrap
  and the full local Symbolica patch are included; Cargo publication is disabled.

## Validation

Post-audit validation passed:

- `cargo test --offline -- --test-threads=1` (from the crate root):
  **31 tests**, comprising 27 unit tests and four integration tests.
- Independent Fortran oracle: **293 points / 879 coefficients**, including 48
  points with nonzero poles; all comparisons passed. Regeneration was
  byte-identical. The corrected timelike triangle agrees at 128/256 bits to
  better than 1e-28.
- Native 128/256-bit polylog endpoint regressions distinguish 1 +/- 2^-80 and
  -1 + 2^-80; 324 sheet operations also pass reference/precision checks.
- Strict all-target Clippy, formatting checks, warning-free rustdoc and source
  whitespace checks passed. Unified-patch context lines are excluded from the
  whitespace check; the bundled patch matches the sibling dependency.
- The unchanged dependency-free numerical core's nine tests passed separately.

The suite includes mixed C0+D0 function-map composition, tensors through rank
four, infrared poles, complex masses, permutation regressions, exact vacuum
limits, undefined-input behavior and native evaluator regressions. These are
not substitutes for the broader missing checks listed below. See STATUS.md for
the exact oracle provenance and coverage inventory.

## Remaining correctness and release work

1. **Known mathematical coverage gap in degenerate C0:** squared invariants `[-4,-1,-1]`, squared
   masses `[0,0,1]` reach an unsupported quadratic-degeneration guard and return
   zero although the finite Euclidean integral is nonzero. A subsequent original
   Fortran check also emits an error and returns zero at this point, so it is
   not a clean successful-Fortran parity counterexample. Other finite-sector
   `a=0` guards are not proven analytic limits. These need a dedicated limit pass.
2. The full boxc alternate reductions, arbitrary Gram/Cayley degeneracies,
   threshold limits, negative-real-part masses, and the full timelike/complex
   mass domain are not established by this sample. Massive box IR continuations
   need broader timelike coverage. Every mass mask having a dispatch entry does
   not establish complete implementation of every kinematic limit.
3. Independent tensor oracle coverage and broad higher-precision box/threshold
   sweeps remain incomplete. Complex native polylogs above approximately 900
   requested bits retain upstream stopping-criterion limitations.
4. Large generic evaluator compilation still needs a large calling-thread stack.
   Cached definitions do not eliminate compilation cost; realistic mixed-amplitude
   performance and Python calling-thread behavior need integration benchmarks.
5. Local Symbolica patches must be accepted upstream or maintained by agreement;
   the entire community extension must use one compatible engine revision.
   The full upstream Symbolica test suite was not run here.
6. Licensing permission remains with the project owner. The manifest identifier
   is provisional, not a claim of permission. Community Python bindings, package
   metadata and final dependency configuration remain release tasks.

No CI, Python registration, package publication, or pull request was added.
