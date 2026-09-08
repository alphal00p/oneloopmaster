# OneLOop performance survey, 2026-09-08

These are sampled warmed-throughput measurements, not universal numerical or
performance guarantees. The original Fortran wrapper/source was not modified.
The benchmark machine was a shared AMD EPYC 9754 Linux host; every timed process
was pinned to logical CPU 160. Other users' work was not stopped. Compiler and
archive provenance are retained in `fortran-build.json` and
`fortran-compiler-trace.txt`. The original generated Fortran source was compiled
with GNU Fortran 15.2.0, `-O2 -fno-fast-math`.

## Final development snapshot

**The requested 1.5×-Fortran target is not met across all workloads.** Correctness
and report-integrity gates pass for these benchmark inputs; that does not remove
the 14 fixed-binary64 B0/dB0 stress-test failures or prove full analytic-region
coverage. See [the scalar audit](../../SCALAR_PARITY_AUDIT.md).

All final routes use 4,097 points per repetition, seven repetitions, 256 warmup
calls and all three Laurent outputs. There are 31 uniform workloads plus five
mixed-family workloads. Batch 1,024 means four full batches **plus a one-row
tail**, not rounding away the remainder. Rust uses release builds; the embedded
evaluators use strict SymJIT O2 with SIMD, scalar retry on divergent branches,
and no fastmath, fast-complex, AVX-512 or internal threading. Exact binary,
dependency-patch, fixture and asset hashes are in [provenance.json](provenance.json).

The routes are deliberately distinct:

- `jit`: a reusable clone of the eagerly prepared scalar-family evaluator.
- `symbols`: public master calls compiled with their transparent native map.
- `hooks-jit`: bare master calls compiled without that map; automatic callbacks
  invoke the already prepared family backends. A four-entry FIFO shares one
  result across distinct Laurent tags, including tag-major SIMD lanes. Each tag
  is consumed once per group, so repeated points are not memoized away.
- `hooks`: the same bare callbacks with an interpreted scalar outer evaluator.
- `python`: the standalone PyO3 list API, using `evaluate_batch` with 1,024 rows.

### Mixed-input results

Ratios below are medians at batch size 1,024. Each denominator is explicit;
lower is better. The Fortran reference performs scalar calls on the identical
point stream. These five mixed workloads contain each family's uniform fixtures
in cyclic order and exercise branch divergence.

| Family | Prepared Rust / Fortran | Mapped Symbols / prepared | Bare Symbol JIT / prepared | Python / prepared | Python / Fortran |
| --- | ---: | ---: | ---: | ---: | ---: |
| A0 | 1.291× | 1.259× | 2.779× | 5.430× | 7.009× |
| B0 | 1.735× | 1.103× | 1.455× | 2.542× | 4.410× |
| dB0 | 2.857× | 1.020× | 1.113× | 1.857× | 5.307× |
| C0 | 3.265× | 0.993× | 1.006× | 1.100× | 3.593× |
| D0 | 9.730× | 1.004× | 1.080× | 1.039× | 10.106× |

For scale, mixed D0 takes 22.413 µs per point in Fortran, 218.084 µs in prepared
Rust, 219.052 µs through mapped Symbols, 235.589 µs through bare Symbol JIT, and
226.517 µs in Python. Mixed A0 takes 47.950 ns in Fortran and 61.899 ns in prepared
Rust, making wrapper/conversion overhead much more visible.

The compiled Symbol routes are close to prepared evaluation for expensive C0/D0
on this sample, **not uniformly for cheap integrals**. Bare callbacks still enter
a scalar family backend per lane; the outer JIT does not turn those callbacks
into one family SIMD batch. At batch size 1, bare-JIT/prepared ratios are
1.882, 1.280, 1.175, 1.004 and 1.067 for A0/B0/dB0/C0/D0. Interpreted-hook/prepared
ratios are 2.061, 1.357, 1.198, 0.995 and 1.011 respectively. Uniform/mixed and
scalar/batched results must not be substituted for one another.

### Workloads meeting the requested bound

Each entry counts workloads with median time at most 1.5× the matched Fortran
median; the denominator is all 36 workloads, not just the fast cases.

| Route | Batch 1 | Batch 1,024 |
| --- | ---: | ---: |
| Prepared Rust | 6/36 | 25/36 |
| Mapped Symbols | 4/36 | 21/36 |
| Bare Symbol JIT | 1/36 | 1/36 |
| Interpreted Symbol hooks | 1/36 | — |
| Python | — | 6/36 |

All five mixed Python workloads miss the bound. Passing a uniform workload does
not establish a family-wide performance guarantee. No aggregate average is used
to hide the slow sectors.

### Timing contract, evidence and reproduction

Python timings include input conversion and result-object allocation, but exclude
input-list preparation, correctness checks, checksums and final result destruction.
Rust and Fortran include a native finite-coefficient checksum in their timed
loops. Every route separately verifies all three fixture coefficients outside
timing; Python additionally checks every timed output afterward. The independent
audit validates all 756 Python Laurent checksums (252 samples × three outputs),
exact JSON/TSV agreement, inventories, repetitions and provenance hashes.

Eager startup and expression/JIT construction are excluded from warmed timings,
not assumed free. Rust startup took about five seconds on this machine, including
all five backends. The separate `*-startup.txt` files record initialization,
prepared-backend cloning and outer compilation times. Python import/startup is
recorded separately in its JSON report.

All benchmark processes from this audit ran sequentially. The shared host was
not isolated: the final Fortran mixed medians were 4.4–10.1% slower than the
earlier same-day run; across all 36 workloads final/earlier medians ranged from
0.839 to 1.114. All final comparisons use the **final** baseline. Few-percent
route rankings and near-1.5× classifications should therefore be treated as
sampled observations, not statistically isolated wins.

Raw final evidence is retained in [fortran.json](fortran.json), [jit.tsv](jit.tsv),
[symbols.tsv](symbols.tsv), [hooks-jit.tsv](hooks-jit.tsv), [hooks.tsv](hooks.tsv),
[python.json](python.json) and [python.tsv](python.tsv). From the repository root:

```sh
python3 tests/compare_performance.py \
  performance/2026-09-08/fortran.json \
  performance/2026-09-08/hooks-jit.tsv --batch 1024
```

The comparator returns **2**, meaning valid measurements that miss the target;
malformed/incomplete measurements return 1. Select batch 1 for `hooks.tsv`.
The final Rust files contain batches 1 and 1,024; compare them individually.
The harness also supports 4/32/256, but these sizes are not part of this final
timing snapshot. The separate correctness suite exercises additional batch sizes.

To repeat the matched runs after the dependency bootstrap, use a verified O2
original Fortran build and an installed release Python extension:

```sh
python3 tests/benchmark_fortran.py /path/to/verified/OneLOop \
  --iterations 4097 --warmup 256 --repetitions 7 --output /tmp/fortran.json
cargo run --release --example performance_survey -- 4097 7 all jit 1,1024 \
  > /tmp/jit.tsv
cargo run --release --example performance_survey -- 4097 7 all symbols 1,1024 \
  > /tmp/symbols.tsv
cargo run --release --example performance_survey -- 4097 7 all hooks-jit 1,1024 \
  > /tmp/hooks-jit.tsv
cargo run --release --example performance_survey -- 4097 7 all hooks 1 \
  > /tmp/hooks.tsv
python3 python/performance_survey.py --iterations 4097 --repetitions 7 \
  --batch 1024 --build-label release --fortran /tmp/fortran.json \
  --output /tmp/python.json --tsv /tmp/python.tsv
```

Pin the timed processes to the same available CPU and avoid overlapping your
own benchmark/build jobs. The Fortran harness compiles a driver, not the existing
library; its flags alone do not establish how that library was built.

The [verification logs](verification/) retain the final **57-pass Rust suite**,
**11-pass Python suite**, Clippy checks, community-mode compile check and asset
regeneration gate. The expanded 256-bit/f64 logs were taken earlier with the same
native formulas/dependency fixes; the later coefficient FIFO and regenerated IR
do not enter that direct-expression audit. Community-host runtime and non-x86
execution are not validated here.

### Remaining optimization work

The native binary64 Li2 specialization removed the previous per-call MPFR cost
without adding a second numerical OneLOop algorithm. The remaining C0/D0 cost
needs profiling and further transparent-expression work. Source review identifies
three bounded experiments, **not implemented or credited with a speedup here**:

1. Bind repeated raw sheet products/quotients and normalization signs once before
   passing them to the three component definitions.
2. Share the D0 root pair's discriminant square root and stable-root selection,
   retaining exact linear/repeated-root and branch-label behavior.
3. Eliminate calls to exactly zero pole definitions in fully massive sectors.

Measure each independently, including compilation memory/time and portable size,
then rerun scalar, sheet, exact-limit, restored-cache and mixed-batch regressions.
Public-wrapper inlining alone cannot accelerate the prepared route, which already
bypasses those wrappers. A packed-buffer Python interface is another possible
follow-up for cheap integrals; the current list API retains conversion overhead.

## Before the native binary64 dilogarithm specialization

`before-li2/` retains a complete, deliberately unsuccessful baseline. It uses
the strict-v3 embedded evaluators, the system allocator, and the previous
Symbolica binary64 polylogarithm callback, which delegates each call to MPFR.
The later Python per-row temporary-vector removal was also not present.

Every workload used the same input stream, all three Laurent outputs, 256
warmup calls, 1,025 points per repetition and three repetitions. Python used
one 1,024-row batch plus a one-row tail. Both implementations passed their
untimed coefficient checks; Python also checked every timed output afterward.
The report includes 31 uniform workloads and five mixed-family workloads.

| Mixed family | Python / Fortran median time |
| --- | ---: |
| A0 | 6.7303× |
| B0 | 4.5791× |
| dB0 | 5.4513× |
| C0 | 199.1422× |
| D0 | 223.8638× |

None of the 36 workloads met the 1.5× target. The slowest uniform workload was
`D0_complex_mu_0.01`: Python 8.776 ms versus Fortran 39.629 µs.
These numbers are **not measurements of the optimized implementation**.

Python timings include input conversion and result-object allocation. Rust
and Fortran timings include a native finite-coefficient checksum; Python
checksums and output validation are outside timing. The reports record this
distinction explicitly. Input-list preparation, initialization and JIT
construction are not included in warmed timings.

The checked-in comparator validates fixtures, inventories, repetitions and
checksums, and correctly returns exit status 2 for this baseline:

```sh
python3 tests/compare_performance.py \
  performance/2026-09-08/before-li2/fortran.json \
  performance/2026-09-08/before-li2/python.tsv --batch 1024
```

This earlier baseline uses 1,025 points × three repetitions, unlike the final
4,097 × seven survey above. The large Li2-related improvement is observable, but
these are not strictly paired before/after runs. Neither this baseline nor
isolated dilogarithm accuracy tests establish integral-wide performance or parity.
