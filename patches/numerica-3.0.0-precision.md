# Numerica 3.0.0 precision patch

[numerica-3.0.0-precision.patch](numerica-3.0.0-precision.patch) applies at the
root of the standalone Numerica repository, against release commit
`e125f2da40f55d5e602d46c9fb7734d8833656c8` (3.0.0). The published crate records
the same commit in `.cargo_vcs_info.json`; application was also checked against
its unmodified source tree.

SHA-256: `af38ad1a20138cc5495a4aa15a6bd9ced72bb81d1ac5be4dde4ec512f1bbf32b`.

## Changes

A tiny complex component can have only two relative precision bits after
cancellation while the other component retains hundreds or thousands. The
release can let that weak component reduce the precision of a well-determined
root, phase or power. This patch:

- Adds a default `Real::complex_sqrt` hook containing the release's existing
  generic algorithm. Float overrides it with scaled Cartesian arithmetic,
  signed-axis handling and explicit nonfinite limits. Computed quantities keep
  their tracked precision; only exact constants use the stronger precision.
- Computes mixed-precision `Float::atan2` as a small correction around an
  accurately constructed axis angle. Small angles retain their available
  relative precision. Axis and nonfinite inputs use the backend conventions.
- Starts nonzero integer powers with the first actual factor, avoiding an
  unnecessary multiplication by a weak-precision identity. Zero powers construct
  the identity at the stronger component precision.
- Evaluates exact powers ±1/2 and ±3/2 through Cartesian square roots and
  products. Negative powers invert before cubing to avoid intermediate overflow.

The patch changes three source files and adds four test files with 13 tests.
Tests cover 512/3456-bit inputs, two-bit cancellation residuals, signed zeros,
all quadrants, zero arithmetic, extreme magnitudes, and nonfinite roots.
Independent MPC square-root references check a component-precision matrix,
axes and inputs as large or small as `1e±10000`. They check that small uncertain
components are not padded with fictitious precision.

Numerica 3.0 already preserves nonzero precision when adding or subtracting
zero. Its `Float::is_negative()` means strictly below zero; sign-bit tests must
use `is_sign_negative()`. The included zero-arithmetic tests pass without source
changes to those operations. After adapting these assertions, six of the eight
original OneLOop precision regressions fail before this patch; two of the four
additional MPC tests also fail before it.

## Apply and test

```sh
git clone https://github.com/symbolica-dev/numerica.git numerica-3.0.0-precision
cd numerica-3.0.0-precision
git checkout --detach e125f2da40f55d5e602d46c9fb7734d8833656c8
git apply --check /path/to/oneloopmaster/patches/numerica-3.0.0-precision.patch
git apply /path/to/oneloopmaster/patches/numerica-3.0.0-precision.patch
cargo test --all-targets
cargo check --no-default-features --features integer-malachite,float-astro
```

The new reference tests require the default MPFR feature; they are gated with
`float-mpfr`. The Astro/Malachite configuration is compile-checked separately.
Validation uses Rust 1.98.1 on x86-64 Linux. The standalone patch passes both
forward and reverse `git apply --check` against a fresh release checkout.
The exact patch applied to that checkout passes **199 tests**: 172 library
tests, 14 existing API regressions and 13 added precision regressions. The
Astro/Malachite compile check and formatting check also pass.

This is a Numerica-only patch artifact. It does not alter OneLOop's dependency
selection or include the separate Symbolica polylogarithm and evaluator API
changes needed for the full Symbolica 3.0 migration.
