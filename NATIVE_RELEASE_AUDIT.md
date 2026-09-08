# Native and arbitrary-precision release audit — 2026-09-08

This development snapshot adds direct generic Rust evaluation for A0/B0/dB0/C0/D0,
Rust/Python arbitrary precision, and complete expression inspection with explicit
resource limits. Native is the manual default and the sole backend of bare master
hooks for f64, DoubleFloat and Float. Explicit manual SymJIT O2 and Expression
evaluation remain available. Both binary64 backend sets initialize eagerly.
The original Fortran wrapper and separate dependency-free numerical port are
unchanged. No CI, publication, community-host registration or upstream PR was added.

## Completed integration gates

The final clean release command is:

```sh
cargo test --release --offline --all-targets --no-fail-fast -- \
  --test-threads=1 --nocapture
```

It passes **114 tests, zero failures**, with six ignored entries. Two ignored
entries are fresh-process helpers invoked by their parent tests; the other four
are opt-in construction/precision/parity diagnostics. The new expanded native
binary64 diagnostic was also run explicitly and **fails on 14 coefficients**;
it is not included among the passes. Logs are retained in
[the verification directory](performance/2026-09-08-native/verification/).

Coverage includes:

- All 293 acceptance points / 879 coefficients through direct native f64,
  DoubleFloat, the original expression route, master hooks and portable caches.
  Native f64 mixed batches include 1/3/4/5/31/256/1024 rows and partial tails.
- The 342-point expanded original-Fortran table at native Float 256 bits, with
  all 1,026 coefficients passing unchanged dimension-normalized tolerances.
  This double-precision oracle is not a 256-bit accuracy reference.
- Genuine complex inputs in all five native families at 3456 bits against
  Expression evaluation at 1120 decimal digits, independently per component
  to better than `1e-1000`; native Float hook checks, elementary independent
  targets and C0/D0 refinement from 1000 to 1120 digits are separate tests.
- Eleven native primitive tests, three independent simple-sector tests,
  mixed-precision roots/phases/powers, exact zero arithmetic, DoubleFloat low
  components, signed axes, extreme magnitudes and subnormal controls.
- Six inspection tests, including 864 real-mass triangle rows; two additional
  complex-inspection tests cover 4,410 rows and twelve 1000-digit controls.
  Successful expansions contain no opaque OneLOop helpers. The specialized
  complete C0 outputs have 220 real-mass or 316 complex-mass nodes.
- Native-only f64/DoubleFloat/Float hooks, error recovery without poisoned
  caches, all Laurent tags, the squared scale argument, cold startup, portable
  restoration, and eight generator regressions including executable lazy-branch
  and common-expression-reuse tests.

The final standalone release Python extension passes **27 tests** (83.755 s):
all five families, Decimal/large-integer boundaries, 32/1000-digit evaluations,
1024-row machine/AP batches, backend selection and rebuilds, cold starts,
thread handoff, expression inspection, and the benchmark harness. Both release
Clippy gates pass with `-D warnings`; formatting and warning-free rustdoc pass.
The community adapter passes compilation, **not runtime inside the host**.

## Precision audit findings

The additional independent audit found and repaired avoidable precision loss in
mixed-component square roots, phases and powers, an overflowing inverse-power
intermediate, and inflated precision metadata from zero arithmetic and norms.
It traced the actual generated power inventory, Decimal conversion, hook caches,
exact constants and reachable Li2 routines. The old approximately 900-bit Li2
stopping cap and finite-value f64 predicate narrowing are removed. The full
[precision audit](ARBITRARY_PRECISION_AUDIT.md) records concrete reproductions,
independent MPC comparisons and the scope of its conclusions.

The first full release run exposed a test-reference accuracy problem at
`s=-1.25e-31`, `m_squared=.125`, `mu_squared=.03125`: the original Expression
reference lost roughly 196 bits to cancellation. Its 384-bit value differed
from a 768-bit reference by `3.3e-55`, whereas the compact 384-bit expression
differed by `4e-84`. Raising only the reference to 512 bits reduced its error to
`1.9e-93`. The dense test keeps its compact working precision and `1e-60`
tolerance unchanged, and adds a 384/512/768-bit convergence regression. Both
the initial failure log and final passing run are retained.

No remaining finite-value narrowing defect was identified in the audited
ArbPrec scalar path. That is not a universal accuracy certificate: fixed working
precision, cancellation, singular conditioning, finite exponent ranges and
already-rounded inputs still matter. Output precision is not padded. See
[PRECISION.md](PRECISION.md) for the API contract.

## Portable assets and performance

Fresh SymJIT regeneration passes all 324 restored acceptance/benchmark rows.
A0/B0/dB0 blobs are byte-identical. C0/D0 differ solely in the serialized order
of a HashSet of function names; canonicalizing that order makes the complete
blobs byte-identical. All numeric constants, instructions and nested payloads
match. The shipped assets are unchanged; a structural comparison and exact
hashes are retained with the verification logs.

The [matched native performance report](performance/2026-09-08-native/README.md)
separates direct Rust, explicit SymJIT, native-backed bare Symbol hooks and
batched Python. It records per-workload misses instead of hiding them in an
aggregate. The earlier SymJIT-default report remains historical evidence.

## Remaining work before publication

- The 14 fixed-binary64 finite-coefficient failures remain: six small-momentum
  B0 and eight dB0 cases, sometimes with severe cancellation error. Native Float
  256 passes these retained controls; fixed binary64 does not adapt automatically.
- Global threshold/Gram/Cayley intersections and every officially supported
  analytic region have not been exhaustively tested. The 293 and 342 tables
  contain 613 unique ordered inputs. Analytic symmetry intentionally supersedes
  inconsistent legacy vacuum-triangle outputs in the documented exceptional case.
- Generic complete C0/D0 expansion can exceed its explicit resource budget.
  The compact FunctionMap and direct numerical routes do not require expansion.
- The sibling dependency setup and bundled patches need a maintainable final
  dependency configuration; the full upstream engine suite, other CPU
  architectures and actual community-host runtime remain unvalidated.
- Licensing permission remains with the project owner. Publication is disabled;
  the provisional manifest identifier is not a permission grant.

This is a verified feature update with explicit limits, not full Fortran parity
or publication approval.
