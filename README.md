# Native OneLOop expressions

An independent Rust library whose only direct dependency is Symbolica. It does
not link the Fortran wrapper or the standalone numerical Rust implementation.
All integral algorithms are native symbolic expressions; continuation definitions
are transparent `FunctionMap` entries. Public master Symbols also have dedicated
evaluation hooks that evaluate those same expressions, never Fortran or a second
numerical OneLOop algorithm.

Development status: the expanded audit passes at 256 bits, but retains 14
fixed-binary64 B0/dB0 disagreements. The measured 1.5×-Fortran performance target
is also not met across all workloads. See `SCALAR_PARITY_AUDIT.md` and the
[performance report](performance/2026-09-08/README.md). Full scalar parity is not
established; this is not publication-ready.

This work is based on **Andreas van Hameren's OneLOop package**. The original
algorithms, analytic continuation conventions and Fortran reference results are
his work; see [Credits and references](#credits-and-references).

Inputs are squared invariants, squared masses with `Im(m²) <= 0`, and a positive
real squared renormalization scale. External invariants are real. Coefficients
are ordered `[finite, 1/epsilon, 1/epsilon²]`, with OneLOop normalization. This is
a Laurent expansion through the finite term, not an unexpanded function of epsilon.

## Development dependency

The manifest currently uses the sibling `../symbolica-dev-v3` checkout for local
fixes: non-inlined evaluator external-function indices, exact native
absolute values on the real axis, and native polylogarithm convergence near the
unit circle, plus reuse of repeated non-inlined definitions before recursive
evaluator compilation. Symbolica is based on dev revision
`fb845d34bda8ccf1fedef6544d3aa46dc24944e3` (latest `dev` inspected on 2026-09-08).
That checkout still declares version 2.2.0; this is preparation for the future
3.0, not a claim that 3.0 has been released. Portable JIT restoration and complex
function registration also require the local patch.
See [patches/README.md](patches/README.md) for a fresh-checkout
bootstrap. Keep Symbolica, numerica and graphica compatible. This temporary path
setup is **not a publishable dependency configuration**. A consuming project's
root manifest must supply these overrides too; library-level patches do not
propagate to consumers.

Both manifests disable Symbolica's default features and explicitly retain
`tracing_max_level_info`, `integer-gmp`, `float-mpfr`,
`native_code_generation`, and `bincode`. This excludes `faster_alloc`, its optional
global mimalloc allocator; the standalone extension uses Rust's system allocator.
The earlier Python sequential-thread crash was in mimalloc during PyO3 argument
extraction, before evaluator construction. The unchanged reproducer passes with
this configuration, without a cache-lifetime or mathematical change.
Cargo features are additive: another dependency enabling Symbolica's defaults
or `faster_alloc` can re-enable mimalloc in an embedding host. Audit the resolved
feature graph at the consuming build root, not just this library's manifest.

## Quick start

Clone this repository as `oneloopmaster`, then follow the sibling Symbolica
[bootstrap instructions](patches/README.md). A recent Rust toolchain supporting
edition 2024 and a working Symbolica license/setup are required (the audited
toolchain is Rust 1.91.1). No Fortran compiler is needed for ordinary use or tests.

```sh
git clone https://github.com/alphal00p/oneloopmaster.git
# Set up the patched sibling symbolica-dev-v3 as described above, then:
cd oneloopmaster
cargo run --example basic
```

[examples/basic.rs](examples/basic.rs) is a complete runnable example. Its core is:

```rust,ignore
use symbolica::prelude::*;

let p = parse!("p");
let mass_squared = Atom::num(1);
let mu_squared = Atom::num(1);
let bubble = oneloop::b0(&p, &mass_squared, &mass_squared, &mu_squared);
let views = bubble.coefficients().each_ref().map(Atom::as_view);
let mut evaluator = Atom::evaluator_multiple(&views, &[p])
    .build()
    .expect("compile B0")
    .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
let mut output = [Complex::new(0.0, 0.0); 3];
evaluator.evaluate(&[Complex::new(-1.0, 0.0)], &mut output);
println!("[finite, simple pole, double pole] = {output:?}");
```

At `p²=-1`, `m0²=m1²=mu²=1`, this gives approximately
`[-0.1520447048200202, 1, 0]`, with zero imaginary parts. The example also prints
the exact A0 expressions. To evaluate complex masses, make the mass arguments
symbolic parameters and supply `Complex::new(real_part, nonpositive_imag_part)`.
Use complex-valued evaluation even for real inputs: timelike integrals can have
an imaginary part. See the precision and stack requirements below.

The scalar entry points are `a0(m², mu²)`, `b0(p², m0², m1², mu²)`,
`db0(p², m0², m1², mu²)`, `c0([p1²,p2²,p3²], [m0²,m1²,m2²], mu²)`,
and `d0([p1²,p2²,p3²,p4²,s12,s23], [m0²,m1²,m2²,m3²], mu²)`.
Here `s12=(p1+p2)²`, `s23=(p2+p3)²`; ordering follows the original OneLOop API.
`db0` differentiates B0 with respect to `p²`, not the renormalization scale.
Arguments are references to Symbolica `Atom`s. The routines construct the three
Laurent coefficients; they do not numerically integrate Feynman parameters.

## Combining integrals

`a0`, `an`, `b0`, `db0` and bubble tensor coefficients return plain native Laurent
series. `c0` and `d0` return `MappedLaurentSeries`, including the definitions needed
to evaluate them. Their `_mapped` names are aliases. For amplitudes, construct a
single `OneLoopExpressions` context and share its definitions:

```rust,ignore
use oneloop::OneLoopExpressions;
use symbolica::prelude::*;

let context = OneLoopExpressions::new();
let p = [parse!("p1"), parse!("p2"), parse!("p3")];
let m = [parse!("m1"), parse!("m2"), parse!("m3")];
let mu2 = parse!("mu2");
let c = context.c0([&p[0], &p[1], &p[2]], [&m[0], &m[1], &m[2]], &mu2);
let combined = c.coefficients()[0].clone() + oneloop::a0(&m[0], &mu2).coefficients()[0].clone();
let parameters = [p[0].clone(), p[1].clone(), p[2].clone(),
                  m[0].clone(), m[1].clone(), m[2].clone(), mu2];
let exact_evaluator = context.evaluator(&[combined], &parameters)?;
// Map exact coefficients with map_coeff for f64, or map_coeff_with_prec for Float.
```

Keep the function map when exporting expressions or building an evaluator through
Symbolica directly. Compact function calls alone are not self-contained formulas.
The context caches the native definitions once per process; cyclic mass sectors
reuse canonical formulas instead of constructing copies.

## Compact master Symbols

`oneloop::{A0(), B0(), dB0(), C0(), D0()}` return the Symbolica
Symbols `oneloopmaster::{A0,B0,dB0,C0,D0}` with dedicated `EvaluationInfo` hooks.
The crate registers a Symbolica `initialize!` hook, ordered after Symbolica's
special-function initializer. At global-state startup it registers all master and
native-helper symbols, prepares their shared definitions, and loads **all five**
SymJIT backends together. Parsing a master does not require a prior accessor call.
Rust crates have no Python-style import event: call `oneloop::initialize()?` at
application startup to force this work and receive startup errors immediately;
otherwise Symbolica's first state access triggers the hook. `is_initialized()` is
a read-only readiness query. Python module import calls `initialize()` explicitly.

Each master call takes a leading Laurent-power tag `0`, `-1` or `-2`, then
the same numeric arguments as the lowercase constructor, flattened into one list.

```rust,ignore
use symbolica::prelude::*;

let p = parse!("p");
let finite = oneloop::B0().call((0, &p, 1, 1, 1));
let pole = oneloop::B0().call((-1, &p, 1, 1, 1));
let context = oneloop::OneLoopExpressions::new();
let exact_evaluator = context.evaluator(&[finite, pole], &[p])?;
```

The context supplies transparent tagged native definitions; Symbolica gives them
precedence over the hooks. Bare calls can instead evaluate through cached hooks
for `Complex<f64>` and `Complex<Float>`. Those hooks evaluate the same native
bodies but are opaque to the *outer* evaluator's optimizer, so prefer the context
for combined amplitudes. They do not supply every Symbolica numeric domain or
automatic error-controlled precision. Invalid tags/arities produce descriptive
errors or panics; numeric callbacks have no `Result` return channel. Stack and
license requirements below apply to either route. See
[MASTER_SYMBOL_AUDIT.md](MASTER_SYMBOL_AUDIT.md) for verification status.

To compile the bare Symbols while retaining their cached evaluation hooks, omit
the function map. This example keeps the squared scale as a runtime parameter:

```rust,ignore
use symbolica::prelude::*;

let parameters = [parse!("p"), parse!("m0"), parse!("m1"), parse!("mu_squared")];
let calls = [0, -1, -2].map(|tag| {
    oneloop::B0().call((tag, &parameters[0], &parameters[1],
                      &parameters[2], &parameters[3]))
});
let exact = Atom::evaluator_multiple(&calls, &parameters)
    .direct_translation(true).build()?;
let mut evaluator = exact.jit_compile::<Complex<f64>>(oneloop::jit_settings())?;
let input: Vec<_> = (0..1024).flat_map(|_| [-1., 1., 1., 4.])
    .map(|value| Complex::new(value, 0.)).collect();
let mut output = vec![Complex::new(0., 0.); 1024 * 3];
evaluator.batch_evaluate(&input, &mut output, 1024);
```

The outer evaluator is compiled here; the five family backends are already ready
from eager initialization. Passing the native map instead allows optimization
across transparent definitions. Both routes are measured separately below.

## Numerical evaluators and batches

`ScalarEvaluator` represents one generic A0/B0/dB0/C0/D0 family and returns all
three Laurent coefficients. Its parameters are runtime inputs, including
`mu_squared` as the **last argument of every row**; the scale is not baked in.
The API preserves manual expression construction and adds
`OneLoopExpressions::jit_evaluator` / `MappedLaurentSeries::jit_evaluator` for
compiling combinations to SymJIT O2.

```rust,ignore
use oneloop::{ScalarEvaluator, ScalarIntegral};
use symbolica::domains::float::Complex;

oneloop::initialize()?; // All symbols and five backends are now ready.
// Clone the prepared O2 backend; no additional IR loading or JIT compilation.
let mut evaluator = ScalarEvaluator::prebuilt(ScalarIntegral::B0)?;
let rows = [-1., 1., 1., 1.,  // p², m0², m1², mu²
            -1., 1., 1., 4.]; // same kinematics, another squared scale
let input: Vec<_> = rows.into_iter().map(|x| Complex::new(x, 0.)).collect();
let mut output = vec![Complex::new(0., 0.); 2 * 3];
evaluator.evaluate_batch(&input, &mut output, 2)?;
// output[0..3] belongs to the first row; output[3..6] to the second.
```

`evaluate` handles one point; `evaluate_batch` accepts flat row-major arrays and
checks their exact dimensions, including empty batches and partial SIMD tails.
For default manual evaluation without constructing an evaluator object,
`oneloop::evaluate(family, input, output)` and
`oneloop::evaluate_batch(family, input, output, rows)` use the same eagerly
prepared cache as the master-Symbol hooks. `ScalarEvaluator::cached(family)`
clones a prepared backend into an independently reusable workspace; an explicitly
compiled custom expression keeps using its own evaluator.
`to_bytes` / `from_bytes` serialize portable SymJIT intermediate code **and nested
function definitions**, not machine code or process pointers. Loading regenerates
host machine code. Cache compatibility is tied to the patched Symbolica/SymJIT
versions; load only trusted artifacts. Float/arbitrary-precision evaluators remain
available through the original exact-expression API and do not use these f64 blobs.
Those evaluators are constructed for the requested precision; eager startup
prepares the five binary64 SymJIT backends, not every possible precision.
Keep `SYMJIT_TOML` unset and do not place `symjit.toml` in the working directory
when building/loading portable caches: these overrides can pin SymJIT to a
particular machine architecture, so the portable API rejects them. Register
any custom numerical Symbol callbacks before loading a manually composed cache.

The `prebuilt` feature is enabled by default and embeds all five evaluators using
`include_bytes!`. Asset provenance, sizes and hashes are in
[assets/evaluators/README.md](assets/evaluators/README.md). To rebuild explicitly,
call `ScalarEvaluator::rebuild` or run:

```sh
cargo run --release --no-default-features --example rebuild_evaluators -- assets/evaluators
```

Without `prebuilt`, startup eagerly builds all five backends from the native
expressions. An optional final family argument writes only that rebuilt asset.
Applications can also
explicitly replace the f64 master-Symbol cache with
`rebuild_cached_evaluator(ScalarIntegral::B0)`. Those hooks use the same numerical
evaluator API, whereas a supplied native function map remains visible to the
outer optimizer. Restored scalar and mixed-batch regression checks pass; the
requested performance target is not met across all sampled workloads. A batch
method does not guarantee a speedup for divergent conditional branches.
The current manifest also patches the published SymJIT 2.24.1 crate in a sibling
checkout for the confirmed direct-power, complex-IF and SIMD defects; see the dependency
[bootstrap notes](patches/README.md). The original Cargo registry is not modified.

## Python (PyO3)

The separate [Python adapter](python/README.md) supplies scalar and batched
complex-number evaluation through PyO3. It adds no Python/PyO3 dependency to the
core Rust crate. It also has an opt-in shared-kernel Symbolica community mode;
that is distinct from a standalone numeric extension and must be built into the
same Symbolica host. See its README for build commands, argument ordering,
rebuilding evaluators, and working examples. The standalone release extension
passes all seven API tests on CPython 3.13.12, including the original sequential
fresh-worker mode, cold eager import, mixed batches, and explicit rebuilding.
The unchanged two-worker allocator regression also passes. These runtime
results do not establish shared-host community integration. Python throughput is
measured separately in the performance survey below.

Build with `cd python && maturin develop --release`. For example, on an
adequately stacked calling thread:

```python
import oneloop_native as olo

bubble = olo.Evaluator("B0")  # Reuses an embedded O2 evaluator.
values = bubble.evaluate_batch([
    [-1.0, 1.0, 1.0, 1.0],
    [3.0, 0.7 - 0.03j, 1.4 - 0.08j, 4.0],
])
# Each row returns (finite, simple_pole, double_pole); final input is mu_squared.
```

Convenience functions `olo.A0/B0/dB0/C0/D0(..., mu_squared=1.0)` use reusable
caches by default and also accept `rebuild=True`. Python inputs/results are
binary64 complex numbers; this list API makes no arbitrary-precision or
Python-throughput guarantee.
The optional [Python performance survey](python/performance_survey.py) defaults
to 1024-row batches and records Python conversion/allocation overhead separately
from initialization. The completed 1024-row survey includes the one-row tail,
all three outputs, and matched Rust/Fortran measurements; see the
[performance report](performance/2026-09-08/README.md).

## Stack and precision

Eager startup and generic C0/D0 compilation traverse a large conditional expression graph.
The Rust regression tests run the entire Symbolica workload on a thread with a 128 MiB
stack. For a standalone program, create that thread **before constructing any
Symbolica values**, and run construction/evaluation inside it. A restricted
Symbolica license permits only one active Symbolica thread per user; do not
overlap Symbolica work across threads or processes under that license.
Applications already using Symbolica must arrange an adequate calling-thread
stack rather than silently spawning another Symbolica worker. Both Python cold
startup checks pass on main threads capped at 8 MiB, but this is not a general
stack bound for arbitrary expressions or explicit rebuilding. Sequential thread
handoff is separately tested; the initializing thread need not remain alive for
the verified cached-evaluator use case.

Use exact rational inputs when testing precision convergence. Mapping coefficients
and supplying values as `Complex<f64>` is fixed-precision evaluation; it cannot
recover digits already rounded out of the input. `Complex<Float>` evaluators can
be built at a chosen precision using `map_coeff_with_prec`. No custom adaptive
precision or near-degeneracy gate is implemented here.

## Verification

Run local checks from this directory:

```sh
cargo test -- --test-threads=1
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

The original acceptance table has 293 points / 879 Laurent coefficients. The
broader retained audit has 342 accepted Fortran points / 1,026 coefficients;
run its fixed-f64 comparison explicitly with:

```sh
ONELOOP_AUDIT_FIXTURES=tests/data/scalar_audit.txt cargo test --test parity \
  audit_extra_fortran_fixtures -- --ignored --test-threads=1 --nocapture
```

The current release suite passes 57 tests, with four ignored test entries
(three opt-in diagnostics and a subprocess helper invoked by the startup test).
The broader table passes at 256 bits after the branch/limit repairs; the fixed-f64
command above still fails on 14
small-momentum B0/dB0 coefficients. No comparison tolerance has been relaxed
and no disagreement is converted into a passing assertion. Passing the sampled
tables or default suite is not a full-parity certificate. The acceptance and
expanded tables have 613 unique ordered inputs, not 635 independent points.
See [the exceptional-vacuum audit](VACUUM_SYMMETRY_AUDIT.md) for the exact disputed
permutations, independent analytic targets, and gaps in region coverage. Analytic
symmetry takes precedence over inconsistent reference outputs there. See the dated audit
for the latest completed runs and additional exact-boundary findings.
The restored strict-v3 scalar caches separately pass all three outputs for the
293 acceptance rows and 31 benchmark rows. The release Python suite passes all
11 tests (seven API/startup/thread tests and four benchmark-harness tests), using
the original fresh-worker mode. Arbitrary-precision, callback, mapped-expression,
portable-cache and batched routes are tested separately, not inferred from one
another. Bare master calls compiled without a FunctionMap pass every acceptance
row in canonical tag order, scalar and five mixed/repeated batch layouts, before
and after portable restoration. Five additional tag orders cover controls from
each family. The exact repeated-root box regression passes at 128/256 bits; the
nearby width-1e-30 diagnostic at 128 bits asserts finiteness only, with its target
comparison performed at 256 bits.

The performance survey compares the same 31 numeric fixtures, all three Laurent
outputs, and 256 warmup calls against the unchanged original Fortran source.
It separates eager startup from warmed evaluation, checks uniform and mixed rows,
and reports batches 1/4/32/256/1024 individually:

```sh
python3 tests/benchmark_fortran.py /path/to/OneLOop --iterations 2000 \
  --warmup 256 --repetitions 7 --output /tmp/oneloop-fortran.json
cargo run --release --example performance_survey -- 2000 7 all jit \
  > /tmp/oneloop-native.tsv
python3 tests/compare_performance.py /tmp/oneloop-fortran.json /tmp/oneloop-native.tsv
```

Use `all symbols` instead of `all jit` to compile public master calls with their
transparent native FunctionMap, using the same batch sizes. Use `all hooks-jit`
to JIT-compile bare master calls without that map, invoking the automatically
registered cached callbacks. `all hooks` keeps that bare outer evaluator
interpreted and is scalar-only (compare with `--batch 1`). These are distinct
measurements; mapped compilation is not the automatic callback path.
Distinct Laurent tags at bitwise-identical inputs reuse one backend result.
A bounded four-entry FIFO accommodates tag-major SIMD lanes, including identical
points, but each tag is consumed only once per group: repeated same-tag calls
cannot be memoized away. Seven counting-backend tests cover every tag order,
successive 1/2/4-lane groups, signed zeros, changed inputs, eviction, explicit
invalidation and failures during evaluation. More heavily interleaved points or
threads can reduce reuse. The outer JIT callbacks still enter the scalar family
backend per point; they do not combine those callbacks into a family SIMD batch.

The optional fifth survey argument selects batches, for example
`4097 7 all jit 1,1024`. This evaluates four full 1024-row batches plus a one-row
tail per repetition. Compare each selected batch separately with `--batch`.
The Python survey in `python/performance_survey.py` defaults to batches of 1024;
its timings include Python input conversion and output-object allocation, unlike
the Rust-only route. Compare runs on the
same machine without competing benchmark/build jobs, with matching call counts
and repetitions. The comparator rejects incomplete or invalid measurements and
returns exit 2 if any selected median exceeds 1.5 times the Fortran time.
The 1.5× target is not met across all workloads; full measurements and limitations
are in the [performance report](performance/2026-09-08/README.md).
The comparator's independent protocol tests
run with `python3 -m unittest discover -s tests -p test_compare_performance.py -v`.

The permanent Fortran fixture generator is standard-library-only Python. Compile
the small test driver against an existing, unchanged original OneLOop build:

```sh
gfortran -I /path/to/OneLOop tests/support/oracle.f90 \
  /path/to/OneLOop/libavh_olo.a -o /tmp/oneloop-oracle
python3 tests/generate_fixtures.py /tmp/oneloop-oracle
```

It prints the fixture table to standard output and rejects rejected/incomplete or
nonfinite Fortran results. The original example program does not print the full
A0/B0/dB0 outputs needed by this table, so use this test driver. The original
library is not a production dependency. Set its compiler runtime library path if required
by your local toolchain. The checked-in fixture table makes normal tests independent
of Fortran, Python, and the standalone numerical port.

For a whole-table high-precision exploratory run, also set
`ONELOOP_AUDIT_BITS=256` (or another precision) in the audit command above.
This selects one numeric domain for every point, not an error-triggered fallback.
`branch_limits`, `triangle_limits` and `box_limits` separately retain focused
precision/branch regressions; see the dated audit for which checks pass.

See `AUDIT.md` and `STATUS.md` for verified coverage and known incorrect or
unsupported cases. This is a development implementation, not a claim of global
Fortran equivalence. Cargo publishing remains disabled. See
`COMMUNITY_INTEGRATION.md` for the inspected Spenso-style integration plan.

## Credits and references

The mathematical and algorithmic basis of this transcription is
[OneLOop, by Andreas van Hameren](https://helac-phegas.web.cern.ch/OneLOop.html).
The original Fortran implementation is the independent reference used by the
scalar parity tests. This repository is a separate, still-incomplete Rust /
Symbolica transcription, not an official replacement or a claim of endorsement
by the original author. The original
[OneLOopBridge wrapper](https://github.com/SecretGmG/OneLOopBridge) remains unchanged.

Please credit the original work when using or building on these expressions:

- A. van Hameren, **OneLOop: for the evaluation of one-loop scalar functions**,
  *Computer Physics Communications* **182** (2011) 2427–2438,
  [arXiv:1007.4716](https://arxiv.org/abs/1007.4716),
  [doi:10.1016/j.cpc.2011.06.011](https://doi.org/10.1016/j.cpc.2011.06.011).
- A. van Hameren, C. G. Papadopoulos and R. Pittau,
  **Automated one-loop calculations: a proof of concept**,
  *JHEP* **09** (2009) 106,
  [arXiv:0903.4665](https://arxiv.org/abs/0903.4665),
  [doi:10.1088/1126-6708/2009/09/106](https://doi.org/10.1088/1126-6708/2009/09/106).

[Symbolica](https://symbolica.io/) supplies the symbolic expressions, function
maps and numerical evaluation infrastructure. It has its own licensing terms.
Permission for the intended distribution/relicensing of this transcription is
being handled separately; it has not been recorded as granted here. The manifest's
license identifier is provisional, not a grant of rights to the original code.
Cargo publication is disabled pending permission and technical readiness.
