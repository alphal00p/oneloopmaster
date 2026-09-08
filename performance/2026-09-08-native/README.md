# Direct Rust performance survey — 2026-09-08

Native Rust is the new default. On the five mixed-input workloads it is faster
than manually selected SymJIT, with native/Fortran ratios of approximately
**0.85 / 0.84 / 1.08 / 1.22 / 1.51** for A0/B0/dB0/C0/D0. This is not a uniform
win over SIMD SymJIT and **does not meet the 1.5×-Fortran target in every sector**.
The earlier [SymJIT-default survey](../2026-09-08/README.md) remains historical;
all ratios here use a freshly measured matched Fortran baseline.

## Protocol

Every route uses the same 31 uniform fixtures and five mixed-family workloads,
all three Laurent outputs, 256 warmup calls, 4,097 evaluations per repetition,
and seven repetitions. Rust routes measure batch sizes 1 and 1,024 separately.
Python uses 1,024: four full batches plus a one-row tail, not a truncated stream.
Fortran performs scalar calls on exactly the same point stream.

All timed processes ran sequentially, pinned to logical CPU 160 of a shared
AMD EPYC 9754 Linux host, after this audit's builds and numerical tests stopped.
Other users' work was not stopped, and the host was not isolated. The original
Fortran library is unchanged, built with GNU Fortran 15.2.0, `-O2 -fno-fast-math`.
Rust uses release builds; manual SymJIT uses strict O2, SIMD with scalar retry
for divergent branches, no fastmath/fast-complex, AVX-512 or internal threading.
The direct Rust family loop is scalar, not a SymJIT evaluator.

Routes are deliberately distinct:

- `native`: reusable direct generic Rust workspace, here instantiated with f64.
- `jit`: manually selected, eagerly prepared portable SymJIT family evaluator.
- `hooks-jit`: bare master Symbols compiled without a FunctionMap; their hooks
  use the direct native backend and share each result across Laurent tags.
- `python`: standalone release PyO3 `Evaluator.evaluate_batch`, explicitly
  selecting `backend="native"` and including list conversion/result allocation.

Mapped transparent Symbols and interpreted outer hooks are not remeasured in
this snapshot. Their old timings must not be relabeled as current native-hook
measurements. Arbitrary-precision throughput is also outside this binary64 survey.

## Mixed-input results

Medians at batch 1,024; lower ratios are better. Times are nanoseconds per point.

| Master | Fortran ns | Native ns | Native / Fortran | Bare Symbol / native | Python / native |
| --- | ---: | ---: | ---: | ---: | ---: |
| A0 | 46.415 | 39.366 | 0.848× | 4.392× | 8.476× |
| B0 | 109.234 | 92.019 | 0.842× | 2.585× | 5.025× |
| dB0 | 124.626 | 134.782 | 1.081× | 2.025× | 3.661× |
| C0 | 2,665.234 | 3,246.944 | 1.218× | 1.131× | 1.257× |
| D0 | 21,170.537 | 31,888.618 | 1.506× | 1.068× | 1.070× |

Mixed D0 takes **21.17 µs in Fortran, 31.89 µs in direct Rust, 34.07 µs through
bare Symbol hooks, 34.11 µs in Python, and 231.44 µs in manual SymJIT**.
Direct Rust is about **7.26× faster than SymJIT** there. The expensive native-hook
route is close to direct evaluation on this sample; that is not true for cheap
A0/B0/dB0, where callback/group-cache overhead matters. Python's mixed
Fortran ratios are 7.189/4.233/3.959/1.532/1.611 respectively: batching avoids
one Python call per point, but does not eliminate conversion/allocation costs.

Native is selected for overall mixed-family throughput, not because it beats
SymJIT on every uniform workload. Uniform complex batches can benefit from
SymJIT SIMD; retain the explicit backend selector when those dominate an application.

## Remaining misses

For direct Rust at batch 1,024, the six misses are:

| Workload | Native / Fortran |
| --- | ---: |
| Timelike real C0 | 1.613× |
| Timelike real D0 | 3.800× |
| Complex D0, nominal scale | 1.518× |
| Complex D0, scale .01 | 1.521× |
| Complex D0, scale 10000 | 1.526× |
| Mixed D0 | 1.506× |

At batch 1, scaleless A0/B0 and one-massless on-shell B0 also miss; timelike D0
is 3.828×. Borderline ratios such as 1.506 are not statistically isolated from
1.5 on this shared host, but they are retained as measured misses, not rounded
into passes. Per-workload counts at or below 1.5 are:

| Route | Batch 1 | Batch 1,024 |
| --- | ---: | ---: |
| Direct native Rust | 27/36 | 30/36 |
| Manual SymJIT | 6/36 | 25/36 |
| Bare Symbol JIT, native hooks | 7/36 | 7/36 |
| Python, native backend | — | 5/36 |

These are workload counts, not global domain guarantees. No aggregate average
is used to hide a slow sector. Correct benchmark outputs also do not resolve
the separate 14 fixed-binary64 stress failures documented in the
[release audit](../../NATIVE_RELEASE_AUDIT.md).

## Evidence and reproduction

Raw reports are [fortran.json](fortran.json), [native.tsv](native.tsv),
[jit.tsv](jit.tsv), [hooks-jit.tsv](hooks-jit.tsv), [python.json](python.json) and
[python.tsv](python.tsv). Exact executable, dependency-patch, fixture and asset
hashes are in [provenance.json](provenance.json); production source hashes are
in [source-sha256sums.txt](source-sha256sums.txt). Compiler provenance is retained
in `fortran-build.json` and `fortran-compiler-trace.txt`.

Rust and Fortran check all fixture coefficients and batch layouts outside timing
and accumulate a finite-coefficient checksum inside timing. Python checks all
three outputs after every timed batch and records all three checksums outside
timing; it excludes prepared input lists and final result destruction. Reports
retain every repetition and exact tail count. These timing-boundary differences
are explicit, not hidden by a single cross-language speed claim.

Eager startup is excluded, not free: direct-Rust survey initialization took
5.33 s, bare-Symbol initialization 5.62 s, and Python import 5.65 s. This prepares
both binary64 backend sets and symbolic definitions. Native per-family constant
setup and outer-JIT construction are recorded in the `*-startup.txt` files.

After the dependency bootstrap and release extension build, pin each process to
the same available CPU and avoid overlapping your own builds/tests/timings:

```sh
python3 tests/benchmark_fortran.py --executable /path/to/verified/bench_oracle \
  --iterations 4097 --warmup 256 --repetitions 7 --output /tmp/fortran.json
cargo run --release --example performance_survey -- 4097 7 all native 1,1024 > /tmp/native.tsv
cargo run --release --example performance_survey -- 4097 7 all jit 1,1024 > /tmp/jit.tsv
cargo run --release --example performance_survey -- 4097 7 all hooks-jit 1,1024 > /tmp/hooks-jit.tsv
python3 python/performance_survey.py --backend native --iterations 4097 \
  --repetitions 7 --batch 1024 --build-label release --fortran /tmp/fortran.json \
  --output /tmp/python.json --tsv /tmp/python.tsv
python3 tests/compare_performance.py /tmp/fortran.json /tmp/native.tsv --batch 1024
```

Ensure the calling thread has sufficient stack for eager symbolic startup;
the timed Python process here used a 128 MiB main-thread stack. Every comparator
route/batch returns **2**: valid data with target misses. Return 1 instead means
malformed/incomplete measurements. Select batch 1 separately for scalar results.

The [verification logs](verification/) retain the clean 114-pass Rust suite,
27-pass Python suite, Clippy/rustdoc/community compilation, the explicitly failing
native binary64 stress diagnostic, reference-precision convergence evidence and
portable-asset regeneration proof. Community-host runtime and non-x86 execution
are not established by this survey.

## If further optimization becomes worthwhile

The implemented inexpensive improvements include scoped arithmetic reuse,
dominating lazy-conditional reuse, branch factoring, scaled primitives and direct
simple-sector formulas. D0's generated source shrank from about 21.9 MB to
1.9 MB, with repeated sorting call sites reduced from 2,034 to two.

The following are **prospects, not demonstrated speedups**:

1. Profile timelike C0/D0 and complex D0 separately, measuring primitive call
   counts and instruction/branch costs before choosing a rewrite.
2. Exploit provably real intermediate quantities in real-mass sectors, avoiding
   unnecessary complex arithmetic; preserve exact axis/lip behavior.
3. Hand-simplify only the hot formulas and share remaining roots, logarithms and
   continuation quantities. If Li2 dominates, optimize its numerical kernel
   separately against independent branch/precision controls.
4. Consider sector-grouped or vectorized native batches only if grouping costs
   and real workloads justify the added complexity.
5. For cheap masters, examine exact pole-tag hook shortcuts and packed Python
   buffers. These address interface overhead, not the raw timelike D0 formula.

No further optimization is assumed necessary merely because a synthetic workload
misses a target. Profile an actual consuming application before expanding this work.
