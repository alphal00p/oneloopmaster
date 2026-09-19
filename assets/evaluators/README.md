# Portable scalar evaluators

Generated on 2026-09-19 with unmodified Symbolica and Numerica at
`821b02451256a92039a0665006628bd5d91470cc` and released SymJIT 2.25.6.
No dependency source patches are applied.

These files contain bincode-serialized numerical `ExpressionEvaluator` graphs,
including nested definitions and registered callback metadata. They contain no
native machine code. Loading recompiles every level for the host with OneLOop's
explicit settings: O2, ordinary translation, packed complex arithmetic enabled,
and fastmath, SIMD batches, threads and AVX-512 disabled. Registered numerical
callbacks preserve complex square-root branches and complex conditions. Packed
complex arithmetic avoids unwanted FMA contraction in the generic compiler.

The generator restored each file and checked all three Laurent coefficients
against all 293 acceptance and 31 benchmark rows before writing it. All five
families passed. Fresh-process, batch and Python validation is tracked in the
[migration notes](../../SYMBOLICA_3_MIGRATION.md). Native x86-64 execution is tested;
other architectures have not been exercised here. These fixtures do not establish
exhaustive analytic-region coverage.

| Family | Bytes | SHA-256 |
| --- | ---: | --- |
| A0 | 1,891 | `7cfe260396001602a94b9b768899a092793c3584dfb50fdd1ed8302264aabd26` |
| B0 | 14,880 | `e22b4a035d092f9ac9d34e7587e81c8c5bdb388713fa558ae67ac74e72ba8ca9` |
| dB0 | 71,901 | `ff7a37d25d2150a26d6d933996f014b664067584944630151516cd2e6957c98c` |
| C0 | 590,631 | `9b6b597e2de2f2cc9873ce0caa1624b533c867baba6eedab1bc78718217542b2` |
| D0 | 13,160,627 | `6d32bfa9a7d3a9032303dd732e9e1e872fb889c8e5920b85a92eb6d19933579e` |

The `oneloop-evaluator-v2` discriminator binds these files to the selected
Symbolica/SymJIT versions and settings. Earlier patched-backend caches are
incompatible. Helper aliases retain their original relative symbol order so
hash-map iteration cannot change the order of floating-point operations.

Regenerate with:

```sh
cargo run --release --no-default-features --example rebuild_evaluators -- assets/evaluators
```

An optional final family name writes only that family. Startup still prepares
all five evaluators. Validation failure leaves the affected existing file
untouched. Commit matching source and regenerated assets together, and update
this manifest's sizes and hashes. Load only trusted serialized evaluators.

The previous generation and its dependency patches are documented in the
[2026-09-08 audit](../../performance/2026-09-08-native/README.md).
