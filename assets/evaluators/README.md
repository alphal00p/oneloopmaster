# Portable scalar evaluators

Regenerated on 2026-09-29 with released Symbolica and Numerica 3.0.1 and
SymJIT 2.26.0. No dependency source patches are applied.

The C0 asset was regenerated on 2026-09-22 after restoring the zero-root
ordering of the finite three-mass triangle. Its round-trip validation passed
all existing C0 acceptance and benchmark rows. An additional regression checks
45 equal-mass kinematic points, including both sides of the timelike threshold,
through the native, expression and portable SymJIT backends; native, expression
and public master callbacks also agree with independent analytic values to
110 decimal digits.

These files contain bincode-serialized numerical `ExpressionEvaluator` graphs,
including nested definitions and registered callback metadata. They contain no
native machine code. Loading recompiles every level for the host with OneLOop's
explicit settings: O2, ordinary translation, packed complex arithmetic enabled,
and fastmath, SIMD batches, threads and AVX-512 disabled. Registered numerical
callbacks preserve complex square-root range/branches and exact predicates.
Packed complex arithmetic retains the validated production policy; the earlier
strict-FMA defect is fixed in SymJIT 2.26.0.

The generator restored each file and checked all three Laurent coefficients
against all 293 acceptance and 31 benchmark rows before writing it. All five
families passed. Fresh-process, batch and Python validation is tracked in the
[migration notes](../../SYMBOLICA_3_MIGRATION.md). Native x86-64 execution is tested;
other architectures have not been exercised here. These fixtures do not establish
exhaustive analytic-region coverage.

| Family | Bytes | SHA-256 |
| --- | ---: | --- |
| A0 | 1,891 | `8a3a160b9e0909a60461179f2327f4bddd85153d79679b4ed1e5b8a73e3edf71` |
| B0 | 14,880 | `cc8611afca806d3d28f16e566da9601963dcf1f4340c99c1f00f6968bf3563c2` |
| dB0 | 71,901 | `5aa2786250a8e131d0b751d7e05b6c89431b55e624dc3e5e07cee9a7e2779f77` |
| C0 | 179,716 | `506ec67a01cc34eaa378b1c3fea3b5092a9a8cd738f3ddac1e0595ad502bb09a` |
| D0 | 4,333,475 | `af637c421c71922dbda4e34c0bb4ee22c72d52c52eba4f438b9963fa5d325869` |

The embedded files retain the `oneloop-evaluator-v2` discriminator. Its old
SymJIT version field describes their producer; loading always recompiles the
Symbolica IR with the current backend. New saves use `oneloop-evaluator-v3`,
which binds the Symbolica IR schema and strict compilation policy without
claiming backend binary compatibility. Only the known v2 source formats are
accepted for migration. Original helper symbols and signatures are preserved. Upstream shares
non-inlined bodies, so explicit-pi aliases and padded argument lists are no
longer needed. Symbol ordering still controls floating-point operation order.

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
