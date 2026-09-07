# Native OneLOop status

The expression crate is independent of the unchanged Fortran wrapper and
standalone numerical Rust port. Symbolica is its only direct dependency.
It is not yet registered as a community Python extension.

**Full scalar parity currently fails.** The expanded scalar-domain audit found
44 coefficient disagreements over 342 accepted original-Fortran points. Known
failures include the regular dB0 pseudothreshold, mixed-sign finite massless C0,
complex finite D0, and massive infrared box continuations. See
`SCALAR_PARITY_AUDIT.md`; the passing acceptance suite below is narrower.

## Expression coverage

- A0 and tadpole tensors through rank four.
- B0, dB0 and bubble tensor coefficients through rank four.
- C0 and D0 dispatch for every mass mask, including the existing infrared sectors.
- The negative-Källén C0 representation and the all-zero-momentum massive
  triangle limit.
- Exact finite two-mass box boundary branches (including vanishing r12/r13),
  and the equal-mass, zero-external-momentum box `1/(6(m²)²)`.
- Continued root/log/Li2 definitions shared through a native FunctionMap.
  No numerical OneLOop callback or Fortran runtime is used.

C0/D0 return compact mapped series; OneLoopExpressions shares one cached definition
map across an amplitude. Cyclic mass sectors reuse canonical definitions.
The implementation is split into two-point, triangle, triangle-continuation,
box, shared-definition and sheet-arithmetic modules.

## Reproducible checks

The post-audit suite passes **31 tests** (27 unit tests and four integration
tests), including all 879 Fortran coefficients. Validation commands and their
final audit results are recorded in `AUDIT.md`.
Tests arrange their own stack without an external RUST_MIN_STACK setting.
Symbolica can emit normalization
warnings while constructing singular expressions in inactive native if branches;
the evaluated active coefficients are required to remain finite by the tests.

The checked-in Fortran suite contains 293 points / **879 Laurent coefficients**:
4 A0, 13 B0, 12 dB0, 125 C0, and 139 D0 points. It includes 48 points with nonzero
Laurent poles, every C0/D0 mass mask, scales between 1e-4 and 1e4, Euclidean and
selected timelike/scattering invariants, complex masses, cyclic permutations and
finite zero-boundary regressions. The comparison rejects NaN/infinity and applies
2e-10 + 2e-8*|reference| **after mass-dimension normalization** (using mu², mu⁴,
or 1/mu² as appropriate). This prevents small, dimensionful box values from
passing with large relative errors. The test enforces fixture inventory, row
shape, valid integral kinds, and mass-mask/pole coverage.

Regeneration with tests/generate_fixtures.py reproduced the checked-in table
byte-for-byte. Oracle provenance:

- OneLOopBridge revision: 402bf0d70bfa2a9d8ae66b8251cc4395b7b5ab99.
- Test driver: `tests/support/oracle.f90`, linked to the unchanged original library.
- Oracle executable SHA-256:
  5c4fee3ec2a2742fe92da699a2ed2f0e5954ef765e1d98b0c449157bb8863667.
- Original library archive SHA-256:
  8b429a58d6e9404608d447d3a7e0a542fc74ed2adb545cc71d35ebeafb8c854b.
- Fixed random seed 0x5EED; inputs and expected coefficients are stored in the table.
- Ordinary tests do not require Fortran, Python or the numerical Rust port.

Additional tests cover the existing A/B/tensor/IR references, isolated continuation
regressions, mixed C0+D0 expressions sharing one map, and 324 low-level sheet
products/quotients (with mandatory inventory checks). Those auxiliary sheet references come from the numerical Rust
transcription. Its truncated dilogarithm series has about 8e-8 error in selected
unit-circle cases: auxiliary Li2 comparisons therefore allow 2e-7, while raw
values/sheet indices/logs remain strict. Independently, the native results are
compared against 128-bit evaluation with a 2e-12 bound. This does **not** relax
the independent Fortran integral tolerance.

## Fixes and native evaluation

The wrong reciprocal root convention and expanded real/imaginary projections
could select an incorrect continuation sheet. Projections and sign selection
now use explicit native function definitions with single evaluated arguments.
This preserves their exact meaning without a floating tolerance.

The sibling Symbolica dev checkout is based on
0b57776bf911faeea7e28ea133706fb03740ffeb and contains local fixes for:

1. External-function index remapping when constants are lifted from a non-inlined
   evaluator.
2. Native absolute values on the real axis, avoiding a square/square-root round trip.
3. Polylogarithm convergence near the unit circle, including avoiding premature
   termination at the zeta function's structural zeros, and exact endpoint/
   real-axis predicates that do not round arbitrary-precision inputs to f64.

These are native Symbolica changes, not replacement numerical integral helpers.
The manifest's local paths must be replaced with an agreed upstream revision or
maintained patch before distribution. The changes are supplied in
patches/symbolica-dev.patch for review.

Large generic evaluator builds need an adequate stack. Tests arrange a 128 MiB
thread for the complete symbolic workload; see README.md for application and
restricted-license threading constraints. There is no hidden worker thread in
the library, no custom adaptive-precision gate, and no new CI.

## What these checks do not establish

Passing the sampled cases is not proof of global Fortran equivalence. In particular,
arbitrary Gram/Cayley degeneracies, all threshold/real-axis limits, arbitrary
negative-real-part masses, and the complete timelike/complex kinematic domain
are not exhaustively covered.

A concrete known incorrect point remains: C0 with squared invariants
`[-4,-1,-1]` and squared masses `[0,0,1]` reaches an unsupported quadratic
degeneration and returns zero although the Euclidean integral is nonzero.
The original Fortran also reports an error and returns zero for this input:
this is a mathematical coverage gap, not a successful-Fortran parity case.
The generic `a=0` guards in other finite sectors are not established analytic
limits. Do not use these degenerate points as supported inputs.

The fully scaleless derivative dB0(0,0,0) is undefined in Fortran and excluded
from finite-value fixtures. Direct construction marks it indeterminate; generic
evaluation returns a nonfinite value rather than a spurious zero.

High-precision tests cover selected cases at 128/256 bits, not a broad stability
proof. Independent tensor fixtures and broad higher-precision D0/threshold
sweeps remain needed. Symbolica's existing complex-polylog stopping criteria
cap their requested-precision threshold at approximately 900 bits; larger
requested precisions are not verified here.

The native closed forms were checked against selected cases that Fortran routes
through its boxc contour implementation. This is not a line-by-line transcription
of every alternate boxc T/T1/T13/planar reduction. Those alternate representations
and special degeneracies remain an explicit audit/coverage item before claiming
100% path equivalence. Higher tensor triangles/boxes and higher Laurent orders
are not promised by this API.

Community integration, dependency/version alignment, final packaging metadata and
the licensing permission being handled by the project owner remain release work.
See COMMUNITY_INTEGRATION.md for the inspected Spenso-style adapter plan.
