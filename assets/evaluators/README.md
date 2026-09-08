# Portable scalar evaluators

Generated from the native expressions by `examples/rebuild_evaluators.rs` on
2026-09-08, using Symbolica dev `fb845d34bda8ccf1fedef6544d3aa46dc24944e3`
and the two local dependency patches documented in `../../patches/README.md`.
These are versioned portable intermediate-code/data files, not native CPU code.
Loading recompiles code for the host and restores nested native definitions.
The later zero-argument SIMD input-transpose repair changes code generation only,
not this serialized IR; these scalar-family assets do not require regeneration
for that repair.

Final regeneration after the arbitrary-precision repairs again passed all
293 acceptance and 31 benchmark rows after restoration; the shipped assets
were kept unchanged. A0/B0/dB0 matched byte-for-byte. C0/D0 differed only in
the serialization order of 12/18 function names from SymJIT's `HashSet` table:
their complete name sets match, and sorting only those entries makes each whole
blob byte-identical. Thus numeric constants, instructions, nested evaluator
payloads and all other metadata are unchanged; export is not byte-deterministic
in this table. See the [regeneration log](../../performance/2026-09-08-native/verification/asset-regeneration.log),
[hashes](../../performance/2026-09-08-native/verification/asset-sha256sums.txt),
and [structural comparison result](../../performance/2026-09-08-native/verification/asset-semantic-comparison.log).
The [read-only comparison source](../../performance/2026-09-08-native/verification/asset-compare-portable.rs)
records the exact binary-layout proof and the temporary paths of that run.

SymJIT settings: O2, direct translation, SIMD enabled, scalar retry on divergent
branches, threads/fastmath/fast-complex/AVX-512 disabled. Native x86-64 execution
is tested; generated code on other architectures has not been exercised here.
The original registry, Fortran library, and numerical Rust port are unchanged.

| Family | Bytes | SHA-256 |
| --- | ---: | --- |
| A0 | 1,437 | `d9e13db197c6684931efb84ad82b83b4e541cfe8744468ed2c880a8cc300e2fb` |
| B0 | 8,968 | `77b9ea18bdf4ccb8a2eb50c7c9b8dc61a70eaeaa61068712d193bd2da1ac467f` |
| dB0 | 39,225 | `1289b1d1e8e9fd2d43d30d3b75f87facf69b4b1f20da5f963d42c8dbb6233874` |
| C0 | 545,750 | `b579a59dbd74e0c8484f51c798f40e505d0eaee56a253f2cbb72f4e0651afde9` |
| D0 | 10,441,150 | `2587812c4f1ea2ec1610b64b4403a3e324bbbc46a9fb28152bc9ddac904579e3` |

The format discriminator identifies the strict `symjit-2.24.1-patched-v3` backend.
The generator reloads each serialized family and verifies all of its rows in
`tests/data/parity.txt` and `tests/data/benchmark.txt` before writing it (293 + 31
rows, all three Laurent coefficients). All five families passed this gate. A
validation failure leaves that family's existing asset untouched. Fresh-process
loads, heterogeneous batch layouts, broader original
fixtures, and Python are separate tests; see the current audit status for their
results. These checks are not proof of exhaustive analytic-region coverage.

From the crate root, rebuild all five files with:

```sh
cargo run --release --no-default-features --example rebuild_evaluators -- assets/evaluators
```

An optional final family name regenerates just one file. Commit matching source,
dependency patches and regenerated assets together; update this manifest's sizes
and hashes too. A backend option, patch or expression change may require a format
revision and rebuilding all assets. Do not load untrusted serialized evaluators.
