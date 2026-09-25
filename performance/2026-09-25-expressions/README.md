# Complete expression construction and late branch selection

These are expression-construction timings, not numerical-evaluation or Fortran
performance comparisons. Measurements use optimized release builds on a shared
AMD EPYC 9754 Linux host, Rust 1.98.1, and unpatched Symbolica revision
`821b02451256a92039a0665006628bd5d91470cc`. Startup, compilation, printing and
numerical validation are excluded from the construction/selection timers.

## Representation and timings

All branches are built before supplying any probe rules. Generic masters use
Symbolica's native `AliasedAtom`: the root and every common-subexpression
definition are included. Definitions contain native operations and references
to other included definitions, not opaque OneLOop helpers or evaluation hooks.
Late selection returns ordinary, still-parametric Symbolica expressions.

Representative release Rust run:

| Complete input | Construction | Late selection, all coefficients |
| --- | ---: | ---: |
| Generic A0, shared | 1.55 ms | 0.03–0.49 ms |
| Generic B0, shared | 0.72 ms | 0.16–0.45 ms |
| Generic dB0, shared | 2.72 ms | 0.23–0.47 ms |
| Generic C0, shared | 175 ms | 4.4–9.2 ms |
| Generic D0, shared | 10.49 s | 0.78–0.89 s |
| `C0(0,-a,a,a,a,b,mu2)`, ordinary expression | 1.98 ms | 0.74 ms |

Selection ranges cover Euclidean, complex-mass and infrared probes. Generic
D0's finite coefficient has approximately 860,000 bindings and 49 MB of Atom
storage (not process RSS). Its selected finite expression has about 676,000
characters in the Euclidean case, 881,000 in the complex-mass case and 323 in
the infrared case; character counts depend on parameter names.

The release Python host independently measured generic D0 construction at
**10.33 seconds** and late selection at **1.12 seconds**. It returned 860,629,
3,219 and 253 definitions for the finite, simple-pole and double-pole
coefficients respectively. The selected Euclidean finite part agreed with the
native backend to the `1e-45` target. No construction-time probe was supplied.

The ordinary C0 case has 1,804 nodes including 52 native `if` nodes. The Python
example constructs it first and selects afterwards, with no `branch_rules`
during construction. At `(a,b,mu2)=(2,1,1)` the result has four dilogarithms,
no remaining conditionals, and gives

```text
finite = -0.313413985805891556410851729326115326931925200907803748216227
simple pole = 0
double pole = 0
```

The under-30-second result for generic masters applies to **shared storage**,
not an unrestricted duplicated tree. Plain `get_expression` remains available
and retains its resource limits. Partial probes, compound patterns and custom
matching restrictions may retain conditions or require expensive unfolding;
they have no universal timing guarantee. Timings are measurements, not hard
deadlines on arbitrary user-supplied argument expressions or hardware.

## Changes that make this possible

- Carry a physical argument and its continued logarithm for inspection sheets,
  eliminating repeated half-turn parity trees. Products/quotients add/subtract
  continued logs. The numerical backends retain their original representations.
- Share argument substitutions as well as completed native subexpressions.
  Sharing only final results was too late to prevent intermediate duplication.
- Memoize expansions without repeatedly charging whole cached subtrees against
  the work budget; preserve lexical scope and predicate context in cache keys.
- Sample branch predicates directly through shared definitions, preserving
  native lazy `if` evaluation and arbitrary precision. Probe values never enter
  the returned parametric bodies. Arithmetic grouping is preserved while deciding
  predicates, avoiding spurious rounding residuals in exact sign tests.
- Treat Python's explicit unconditional replacement condition and RHS-cache
  setting as equivalent to Rust's ordinary defaults for DAG probing. Matching
  restrictions themselves are not discarded.
- Use an analytically integrated four-dilogarithm formula for the reported C0
  family, including its zero-mass dispatch and real-axis cut prescriptions.

## Validation and scope

The **29 focused Rust tests passed** (six unit tests and 23 integration tests).
They cover:

- Every generic master, every Laurent coefficient, three probe regions,
  selection only after complete construction, and no remaining `if` nodes.
  Comparisons to the independent native Rust backend use 384 bits and an
  absolute `1e-45` target. This is 45 coefficient comparisons, not exhaustive
  coverage of all analytic regions.
- C0 full-versus-selected comparisons at both signs, zero masses, and a
  discriminant-zero boundary; the supplied positive reference value to `1e-65`.
- 5,880 old-versus-new sheet component comparisons, covering both real-axis lips,
  complex values, products, quotients, and repeated sheet operations.
- Existing one-scale triangle grids: 864 real and 4,410 complex-mass points;
  a separate 12-point comparison retains a `1e-1000` target.
- Deeply nested conditions, integer/float/fraction probes, partial replacements,
  native pattern restrictions, independent caches, and native compilation of an
  `AliasedAtom` coefficient.

The **15 focused Python tests passed in 18.673 seconds**, covering ordinary C0
and generic shared D0 construction followed by selection, shared-expression
inspection, native replacements, and the existing Decimal/precision API. No
claim of exhaustive scalar parity follows from these tests. In particular, an
additional 120-point exploratory C0 sweep did not validate eight comparisons
against the native reference: it returned NaNs for `b=0` at `a=-1,1,7`, and
precision-dependent values for `a=-1,b=1/8`, each at scales 1 and 7. These are
unresolved reference/continuation diagnostics, not evidence of a Fortran bug.
The native backend was not changed in this work.

The complete repository test suite was not run. Focused suites and scoped Clippy
checks are used here; the latter allow the pre-existing
`chunks_exact_to_as_chunks` lint in unchanged numerical batching code. Tests run
serially because the supplied legacy-format key was rejected by the pinned
Symbolica revision and the restricted fallback permits one active thread.

## Reproduce

From the repository root, with the shared-kernel Python host installed:

```sh
cargo run --release --example expansion_survey -- C0-case
cargo run --release --example expansion_survey -- D0 10000000 shared
python examples/c0_selected_branch.py --full

RUST_MIN_STACK=134217728 cargo test --release --lib inspection -- --test-threads=1
cargo test --release --test full_expansion --test branch_expansion --test branch_selection --test inspection --test inspection_complex --test inspection_master -- --test-threads=1
PYTHONPATH=python/tests ONELOOP_PYTHON_MODULE=symbolica.community.oneloop python -m unittest test_inspection test_precision -v
```

See the [inspection guide](../../EXPRESSION_INSPECTION.md) for the public Rust and
Python APIs and their limits.
