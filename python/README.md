# Python adapter

This development adapter keeps Python dependencies outside the Rust `oneloop`
core. It offers two distinct build modes:

- A standalone `oneloop_native` extension accepting Python numbers.
- An opt-in community binding linked into the host's single `symbolica.core`
  extension, also accepting and returning that host's Symbolica expressions.

Matching Symbolica revisions in two independently loaded extension binaries does
**not** establish shared expression state. The standalone extension deliberately
does not exchange Symbolica `Expression` objects. This is not a published package
or a completed integration into the community repository.

## Build the standalone extension

Use a Python virtual environment with Maturin installed:

```sh
cd python
maturin develop --release
```

Python 3.9 or newer and Rust edition 2024 are required. The manifest pins the same
Symbolica revision as the core, `fb845d34bda8ccf1fedef6544d3aa46dc24944e3`, and
uses PyO3 0.28. Its development patches expect the sibling checkout
`../../symbolica-dev-v3`, including its `numerica` and `graphica` subcrates.
They also use `../../symjit-2.24.1-patched` for the tested local backend fixes.
Those paths must be adjusted at the consuming Cargo build root if the checkout
layout differs. Patches in a dependency's manifest do not propagate to its host.
The supplied lockfile records the standalone build's dependency resolution.

The core and adapter dependencies use `default-features = false` for Symbolica,
explicitly retaining `tracing_max_level_info`, `integer-gmp`, `float-mpfr`,
`native_code_generation`, and `bincode`. Symbolica's optional `faster_alloc`
global mimalloc allocator is excluded; this standalone build uses Rust's system
allocator. The earlier sequential-thread crash was traced to mimalloc during
PyO3 constructor-string extraction, before evaluator construction. With this
supported feature configuration, the unchanged reproducer and full API suite
pass; no JIT-clone fix or extra cache/thread-lifetime restriction was required.

Import eagerly registers every native symbol and the complete transparent
function map, then prepares all five Native and all five SymJIT caches. With the core's default
`prebuilt` feature, this loads portable intermediate code and regenerates host
machine code; it does not deserialize foreign native code. Without that feature,
startup builds all five evaluators from the expressions instead. A missing or
incompatible enabled cache raises `RuntimeError` during import. Subsequent
ordinary calls use the direct Rust backend. `backend="symjit", rebuild=True`
explicitly recompiles the requested SymJIT evaluator and is substantially more
expensive than warmed evaluation. Native rebuilding prepares fresh constants
and a numeric workspace, without JIT compilation.

The core enables its five embedded cache assets by default. The final release
extension with the system allocator and native binary64 Li2 passes all eleven
standalone tests on CPython 3.13.12: seven API tests and four pure performance-
harness tests. These include cold import, mixed batches, rebuilding, and
sequential fresh test threads. Batch performance has also been measured below;
the overall 1.5× Fortran target is not met. Community mode is compile-checked
only; runtime validation inside an integrated host remains pending.
For these portable scalar APIs, unset `SYMJIT_TOML` and use a working directory
without a `symjit.toml` override. The wrapper rejects either override rather than
silently changing the cached compiler configuration.

## Numerical API

```python
import oneloop_native as olo

assert olo.is_initialized()  # Read-only: import has already prepared all families.
finite, simple_pole, double_pole = olo.A0(2.0, mu_squared=4.0)
value = olo.B0(-1.0, 1.0, 1.0, mu_squared=1.0)

evaluator = olo.Evaluator("B0")
value = evaluator.evaluate([-1.0, 1.0, 1.0, 1.0])
values = evaluator.evaluate_batch([
    [-1.0, 1.0, 1.0, 1.0],
    [7.0, 1.0, 2.0, 1.0],
    [3.0, 0.7 - 0.03j, 1.4 - 0.08j, 4.0],
])
```

Every call returns three coefficients ordered as
`(finite, simple_pole, double_pole)`, corresponding to Laurent tags `(0, -1, -2)`.
Ordinary numeric inputs with the default `prec=16` use binary64 components and
return Python `complex` values through direct Rust arithmetic by default.
Select `backend="symjit"` for the prepared SymJIT O2 path, or
`backend="expression"` for the transparent Symbolica expression interpreter.
Higher precision or Decimal inputs use the arbitrary-precision path below.
The direct functions have these argument orders:

| Function | Ordered arguments before optional `mu_squared=1.0` |
| --- | --- |
| `A0` | `mass_squared` |
| `B0`, `dB0` | `momentum_squared, mass_0_squared, mass_1_squared` |
| `C0` | `p1_squared, p2_squared, p3_squared, mass_0_squared, mass_1_squared, mass_2_squared` |
| `D0` | `p1_squared, p2_squared, p3_squared, p4_squared, s12, s23, mass_0_squared, mass_1_squared, mass_2_squared, mass_3_squared` |

Here `s12` and `s23` are squared channel invariants, not their square roots. All
momenta and masses are already squared quantities. External invariants must be
real; complex squared masses must have nonpositive imaginary parts; `mu_squared`
must be positive and real. Nonfinite input components are rejected. These checks
do not establish that every singular kinematic point has a finite integral or
that every analytic region has been validated.

The reusable `Evaluator` accepts the same ordered arguments **with
`mu_squared` always supplied last**, giving arities 2, 4, 4, 7, and 11. An empty
batch returns `[]`; ragged rows or incorrect argument counts raise `ValueError`.
Rows are passed together to the core's batched evaluator. List conversion and
Python result allocation remain part of Python call overhead and are included
in the measured Python/Fortran ratios below. NumPy is not needed in standalone mode.

Direct functions use the core's shared, mutex-protected family caches; lowercase
aliases `a0`, `b0`, `db0`, `c0`, and `d0` are available. Each reusable object
clones a warmed backend into its own evaluator workspace:

```python
evaluator = olo.Evaluator("B0", rebuild=True)  # Explicit fresh construction.
evaluator.rebuild()  # Replace this object's evaluator only after success.
value = olo.B0(-1.0, 1.0, 1.0, rebuild=True)  # One uncached rebuilt call.
assert evaluator.family == "B0"
assert evaluator.arity == 4
```

### Arbitrary precision

All five direct functions accept `prec=16`, `prec=32`, `prec=1000`, or another
positive integer number of **decimal digits**. The reusable evaluator accepts
the same keyword, and its methods can override it for one call:

```python
from decimal import Decimal, localcontext
import oneloop_native as olo

finite, pole, double_pole = olo.A0(Decimal("2"),
                                 mu_squared=Decimal("4"), prec=1000)
assert isinstance(finite.real, Decimal)
assert isinstance(finite.imag, Decimal)
print(finite.real)  # Approximately 3.386294361119890618834464242916...

bubble = olo.Evaluator("B0", prec=32)
row = [Decimal("-1"), Decimal("1"), Decimal("1"), Decimal("4")]
values = bubble.evaluate_batch([row] * 1024)
more_digits = bubble.evaluate(row, prec=1000)  # Does not change bubble.prec.
assert bubble.prec == 32

# Decimal arithmetic after evaluation follows Python's decimal context.
with localcontext() as ctx:
    ctx.prec = 1000
    twice_real_part = finite.real * 2
```

`Decimal("...")` is parsed directly at the requested working precision, including
`mu_squared`: there is no intermediate float conversion. Creating a Decimal
from a string and creating the output Decimal objects do not depend on the
global decimal context. In contrast, `Decimal(0.1)` already contains the Python
float's rounding error; use `Decimal("0.1")` to specify that decimal value.
Python float/complex inputs at higher `prec` retain their original binary64
values exactly, but cannot supply missing input digits.

Python's built-in `complex` cannot hold Decimal components. Use
`DecimalComplex(real, imag)` for arbitrary-precision complex squared masses:

```python
mass = olo.DecimalComplex(Decimal("0.70000000000000000000000000000001"),
                          Decimal("-0.03000000000000000000000000000002"))
finite, _, _ = olo.B0(Decimal("3"), mass, Decimal("1.4"),
                      mu_squared=Decimal("4"), prec=1000)
print(finite.real, finite.imag)
rounded = complex(finite)  # Explicitly lossy conversion to binary64.
```

High-precision calls return three `DecimalComplex` instances with standard-library
Decimal `.real` and `.imag` properties, including for real-valued results.
They are precision-preserving value containers, not a replacement complex
arithmetic package; operate on their Decimal components under a suitable context.
No additional Python numerical dependency is required.

Selection rules are consistent across scalar and batch calls:

- `prec != 16` always selects arbitrary precision.
- At `prec=16`, any `Decimal`, `DecimalComplex`, or integer beyond the exact
  binary64 integer range selects arbitrary precision too.
- A mixed batch uses arbitrary precision for the entire batch if any input
  requires it. Otherwise the default batch uses binary64 Native arithmetic.
- Arbitrary-precision results use DecimalComplex; ordinary default calls retain
  Python complex outputs. Precision must be a positive integer, not a Boolean.

The core uses the generated generic Rust implementation with Symbolica
`Complex<Float>` and fixed guard bits. `backend="expression"` evaluates the same
exact formulas with Symbolica's interpreter; `"symjit"` is binary64-only.
The reusable object keeps a precision-specific workspace for reuse; changing
precision prepares its constants again. Both sets of five binary64 backends are
loaded eagerly on import, but arbitrary precisions are prepared on request.
Arbitrary-precision batches amortize Python overhead; they do not use binary64
SymJIT SIMD. `rebuild=True` rebuilds the selected backend, and no custom
automatic precision-escalation policy is added.

Requested precision is a working-precision choice, not a certified error bound.
Cancellation, exact singularities, or unverified analytic limits can still affect
results. Returned decimal strings are limited by the computed Float precision,
not padded to imply that lost digits were recovered. The historical binary64
performance ratios do not apply to this backend. See the core
[precision contract and checks](../PRECISION.md).

### Calling thread, stack, and license

The adapter performs all work on its calling thread. It never silently spawns a
worker, changes numerical precision, or changes Symbolica's license policy.
An `Evaluator` must be used on the thread where it was created.

Eager import includes expression-map construction and must run on an adequately
sized stack, as must explicit rebuilding. Most regression checks use 128 MiB;
the two ordinary-main-thread startup regressions passed with an 8 MiB limit in
the final system-allocator release run. For explicit rebuilding in a new
standalone application, arrange an adequate stack **before importing the extension or
otherwise starting Symbolica work**; for example:

```python
import threading

def work():
    import oneloop_native as olo
    print(olo.Evaluator("B0", rebuild=True).evaluate([-1, 1, 1, 1]))

threading.stack_size(128 * 1024 * 1024)
worker = threading.Thread(target=work)
worker.start()
worker.join()
```

Applications already using Symbolica should arrange an adequate calling-thread
stack and observe its licensing constraints. The installed restricted license
permits only one active Symbolica
thread per user; do not run this concurrently with other Symbolica tests or
applications under that license. A valid Symbolica setup is still required.

## Shared-kernel community integration

This mode implements `SymbolicaCommunityModule`; it does not edit or register
anything in an existing community checkout. A compatible host can depend on the
adapter with:

```toml
oneloop-python = { path = "../oneloopmaster/python", default-features = false, features = ["community"] }
```

The host must resolve a single compatible Symbolica/PyO3 dependency graph and
apply the matching Symbolica, numerica, graphica, and SymJIT patches at its own Cargo
build root. Link this Rust library into that host; do not also load a standalone
OneLOop wheel as an expression bridge. Cargo unifies dependency features:
another host dependency enabling Symbolica defaults or `faster_alloc` can
re-enable its global mimalloc allocator despite this adapter's manifest. The
host must audit its resolved feature graph, preserve the intended allocator
configuration, and rerun the thread-handoff/API checks in that actual host.
In the host's existing registration code:

```rust
register_module!(m, oneloop_native::CommunityModule);
```

The corresponding Python facade follows the host's existing pattern:

```python
from symbolica.community.oneloop_native import *
initialize_module()
```

Registration adds Python classes and functions without constructing Symbolica
symbols. `initialize_module()` calls the core's eager initializer, preparing
all symbols, the complete native expression map, and all five scalar caches.
The core also registers this work with Symbolica's `initialize!` inventory, so
the host's first Symbolica state access may have already completed it. This is
idempotent, and `is_initialized()` reports readiness without doing more work.

In this single-kernel configuration, compact master expressions can be combined
with other host expressions and then compiled with their transparent definitions:

```python
from symbolica import S
from symbolica.community import oneloop as olo

x = S("x")
finite, pole, double_pole = olo.master_coefficients("A0", [x, 1])
evaluator = olo.compile_native([finite + pole], [x])
result = evaluator.evaluate_complex([2 + 0j])
```

`master_coefficients(family, arguments)` includes the explicit squared scale last
and returns three genuine host `Expression` objects. Use `compile_native` for
combined expressions so OneLOop's native function definitions remain visible to
the evaluator compiler; simply returning a compact master call is not a substitute
for supplying that map. `compile_native` returns the host's `Evaluator`, initially
using its complex interpreter, with exact rational constants retained. Its
optional JIT settings start with the core's branch-preserving O2 configuration;
enabling JIT is explicit through the host evaluator API. The host's array output
and numerical precision interfaces are separate from the standalone list API.

## Validation status and tests

The adapter tests cover elementary values, all 31 shared benchmark fixtures,
single/batch agreement for 3/4/5/8-row B0 layouts, heterogeneous family batches
of size 4/5/8 including partial tails, explicit rebuilding, and argument errors.
The shared-kernel expression test is available for an integrated host, but the
standalone run checks rejection of foreign expression interoperation instead.
These tests are not evidence of exhaustive analytic-region coverage or of a
Python performance target.

On 2026-09-08, all **eleven standalone tests passed in 19.740 seconds** with
CPython 3.13.12, the final release extension, the system allocator, native
binary64 Li2, and the original fresh-worker-per-test mode
(`ONELOOP_PYTHON_THREAD_HANDOFF=1`). Seven are actual API tests; the other four
are fake-evaluator performance-harness tests that do not run Symbolica. This
supersedes the earlier 107.040-second seven-test debug result.

The standalone suite first tests eager import in two fresh subprocesses, one
entering through a direct function and one through a reusable evaluator after
import. Each checks readiness before its first master call, all 31 fixtures
across all five families, and nine mixed batch rows per family, with no rebuild
fallback or worker thread. On POSIX each child caps its inherited main-thread
stack at 8 MiB without increasing a smaller limit. Separate children verify rejection of
`SYMJIT_TOML` and a working-directory `symjit.toml`, before cache loading or
source construction. They use private temporary directories and child-only
environment changes. Both cold children and both environment guards passed in
the full system-allocator run. This establishes those specific startup cases,
not a general stack-size bound for arbitrary expressions or rebuilding.

After building and installing the extension, run from the Rust crate root
(`..` relative to this directory):

```sh
python -m unittest discover -s python/tests -v
```

To force explicit rebuilding in the numerical checks, use:

```sh
ONELOOP_PYTHON_REBUILD=1 python -m unittest discover -s python/tests -v
```

This test mode passes explicit `rebuild=True` to numerical constructors and
direct calls. It does not disable eager startup or change the core's enabled
features: import still loads embedded caches in the default build. Batch sizes
make the same tests applicable when SIMD is enabled, but do not assert that it
is enabled. Explicit-rebuild mode skips the two positive cold subprocesses;
portable-environment rejection is still checked at import in separate children.
Community mode skips the standalone subprocess tests entirely.

To test a prepared community host, set
`ONELOOP_PYTHON_MODULE=symbolica.community.oneloop` instead.
After the sequential subprocess checks, the remaining test harness explicitly
creates one persistent, adequately stacked test thread and does not require any
third-party Python testing dependency. The unchanged separate
`tests/thread_handoff.py` regression also passes: after the initializing worker
exits, a new worker constructs and drops all five family evaluators. The complete
API suite also passes with a fresh worker per test. The native allocator trace
and configuration repair are retained in `tests/thread_handoff_backtrace.txt`.
These results do not impose a lifetime restriction on the initializing thread;
individual Python `Evaluator` objects still use their documented creating-thread
policy. Validation inside an actual community host remains pending.

To rerun the confirmed sequential-thread checks separately:

```sh
python -X faulthandler python/tests/thread_handoff.py
ONELOOP_PYTHON_THREAD_HANDOFF=1 python -X faulthandler -m unittest \
  discover -s python/tests -p test_api.py -v
```

## Python batch performance survey

The separate standard-library survey defaults to 1024 rows per batch, 20,000
scalar points per repetition, seven repetitions, and a 256-point untimed scalar
warmup per workload. It uses the same 31 named fixtures and explicit squared
scale arguments as the Rust and original-Fortran surveys:

```sh
python python/performance_survey.py \
  --batch 1024 --iterations 20000 --repetitions 7 \
  --build-label release --output /tmp/oneloop-python.json \
  --tsv /tmp/oneloop-python.tsv
```

Run from the Rust crate root with the desired extension already installed.
`--build-label` is a caller-supplied description, not proof of compiler flags;
the report also records the module path and binary hash. Import and all numeric
work stay on the main thread. No source-rebuild fallback is used.

`SAME` repeats each fixture separately; `HETERO` cycles a family's fixtures in
file order. The report records exact point counts, batch counts and partial
tails, raw repetitions and medians, Python/CPU metadata, and fixture hashes.
All three coefficients are validated both before measurement and for every
measured output. Python input lists are built before timing. The timed interval
includes calls to `Evaluator.evaluate_batch`, adapter list-to-Rust conversions,
result allocation, and retention of each returned batch. Checksums, validation,
and final output destruction are outside timing. In contrast, the native Rust
and Fortran measurements include their finite-coefficient checksum accumulation
inside timing; reports explicitly preserve that difference.

Supply `--fortran /tmp/oneloop-fortran.json` for per-workload ratios with exact
matching fixture inputs, call counts, repetitions and warmup. Exit status 2 then
means at least one selected Python/Fortran median ratio exceeds
`--maximum-ratio` (default 1.5); it is not a correctness failure. Its four
fake-evaluator unit tests check counts, tails, output validation and report generation without
loading Symbolica:

```sh
python -m unittest discover -s python/tests -p test_performance_survey.py -v
```

### Measured release results (2026-09-08)

The final survey used 4,097 points per repetition, seven repetitions, and
1,024-row batches: four complete batches plus a one-row tail. All three Laurent
outputs passed validation for all 31 uniform-fixture workloads and five mixed
family workloads. The measured build uses the system allocator, native
binary64 Li2, regenerated strict-v3 assets, and the four-entry coefficient
grouping cache, with warmed O2 evaluators. Each timed process was
pinned to CPU 160 on a shared AMD EPYC 9754 Linux host; other users' work was
not stopped.

The overall ≤1.5× Fortran target **was not met**: six of the 36 workloads met
it and 30 missed it. All five mixed-family workloads missed that target.
Ratios below compare median time per scalar point, with all three outputs:

| Mixed family | Python / original Fortran | Python / prepared Rust, same batch |
| --- | ---: | ---: |
| A0 | 7.0092× | 5.4296× |
| B0 | 4.4102× | 2.5419× |
| dB0 | 5.3069× | 1.8572× |
| C0 | 3.5930× | 1.1003× |
| D0 | 10.1065× | 1.0387× |

The prepared-Rust column measures the same cached scalar evaluator route and
batch size, not automatic symbolic callbacks or expression-map compilation.
Small Python/Rust timing differences do not mean either is near the original
Fortran speed, and shared-host variation limits fine-grained comparisons.
All ratios use the final, matching Fortran baseline, not the earlier snapshot.
Startup, explicit rebuilding, arbitrary precision, and an
integrated community host are outside these measurements. The timing-boundary
differences described above remain relevant; these are sampled throughput
results, not universal guarantees. See the
[full report and retained evidence](../performance/2026-09-08/README.md).
