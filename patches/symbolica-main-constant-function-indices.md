Upstream status: fixed by `65ff188`, included in OneLOop's selected main
`821b02451256a92039a0665006628bd5d91470cc`. The patch below is historical.

# Symbolica: preserve callback indices when lifting constants

[symbolica-main-constant-function-indices.patch](symbolica-main-constant-function-indices.patch)
applies at the Symbolica repository root, against main commit
`3cb73a6d958755eaa86088ff713dd13e6bda217b`.

Revalidated unchanged against main `08defb398420ab07a428876d4fd10eaf6016470f`:
the regression still fails before the fix and passes afterward. See the
[remaining-regressions report](symbolica-main-regressions.md) for the other
Symbolica patches and current validation.

## Problem and fix

When compiling a non-inlined function, Symbolica hoists its coefficients and
registered constants into the caller and supplies their values as additional
arguments. This keeps coefficient conversion and arbitrary-precision evaluation
consistent with inlined functions.

The sub-evaluator then removes the hoisted constants from `external_fns`, but
its existing `Instr::ExternalFun` instructions still contain indices into the
original table. Those indices can select the wrong callback or go out of bounds.
The regression reproduces a panic in instruction export: the filtered table has
two entries, but a call still refers to index three.

The fix computes an old-to-new index map for surviving callbacks and rewrites
the call instructions before filtering the table. Constant hoisting, symbolic
tags, fixed arguments, and numeric argument order retain their existing behavior.
No public API or serialized format changes are needed.

## Regression coverage

The added test exercises nested non-inlined functions, a tagged helper, two
distinct tagged constants, and two distinct tagged callback implementations.
It covers mixed tables, tables containing only constants, and tables without
constants. Each case checks instruction export and numerical evaluation, plus
bincode round-trips and JIT compilation when those features are enabled.

The regression fails on the unmodified base with the out-of-bounds panic and
passes with the patch, including serialization and JIT checks. All 20 existing
`evaluate::test::` tests also pass, including nested sub-evaluators, coefficient
conversion, C++ export, and JIT compilation. Forward and reverse
`git apply --check`, Rust formatting, and `git diff --check` pass.

## Apply and reproduce

In a clean Symbolica checkout at the commit above:

```sh
git apply --check /path/to/symbolica-main-constant-function-indices.patch
git apply /path/to/symbolica-main-constant-function-indices.patch
RUST_MIN_STACK=134217728 cargo test --no-default-features \
  --features integer-gmp,float-mpfr,native_code_generation,bincode \
  --test lifted_constant_indices -- --test-threads=1
RUST_MIN_STACK=134217728 cargo test --no-default-features \
  --features integer-gmp,float-mpfr,native_code_generation,bincode \
  --lib evaluate::test:: -- --test-threads=1
```

Validation uses Rust 1.98.1 on x86-64 Linux and unmodified SymJIT 2.25.6.
This is a standalone upstream patch artifact; it does not change OneLOop's
dependency selection or apply the patch to its dependencies.
