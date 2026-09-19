# Numerica radial precision patch

Upstream status (2026-09-19): Symbolica includes the `hypot` repair in
`1140960b` and a further complex-logarithm constant-precision fix in `08defb39`.
OneLOop now selects `821b02451256a92039a0665006628bd5d91470cc` and uses these
operations directly. The artifacts below retain their original base revisions;
do not apply them to the updated Symbolica checkout.

These are two formats of the **same fix**; choose the one matching the checkout:

| Patch | Base |
| --- | --- |
| [numerica-hypot-precision.patch](numerica-hypot-precision.patch) | Standalone Numerica `e125f2da40f55d5e602d46c9fb7734d8833656c8` (3.0.0 and current standalone main when checked). |
| [symbolica-main-numerica-hypot-precision.patch](symbolica-main-numerica-hypot-precision.patch) | Symbolica `3cb73a6d958755eaa86088ff713dd13e6bda217b`, modifying only `lib/numerica/`. |

Both patches change `Real::hypot` and add eight regression tests. No OneLOop
manifest applies either patch; its dependencies remain unmodified upstream.

## Problem and fix

The scaled formula is `a * sqrt(1 + (b/a)^2)`, with `a >= b >= 0`. Previously,
`1` was constructed using `(b/a).one()`. A tiny component can have only two
relative precision bits after cancellation. Giving that precision to the exact
identity degrades the whole magnitude, even when the dominant component remains
accurately known. This also affects complex magnitude and the radial part of
complex logarithms.

For example, take a 128-bit real component `2` and an imaginary component
`2^-120` stored at two-bit precision. Upstream computes the real logarithm at
only five bits (approximately `0.688`), instead of retaining an accurate
`log(2) = 0.6931471805599453...`.

The patch constructs only the exact identity from the input with more precision.
It does not increase the precision of either input, their ratio, or any computed
result. A weak dominant component, or a comparably large weak component, still
limits the output precision. The existing scaled arithmetic and zero handling
remain in place.

## Regression coverage

The eight tests cover:

- 128-, 512- and 3456-bit dominant components paired with two-bit tiny components,
  including swapped inputs and both signs.
- Weak dominant and comparable components, which must not acquire fictitious
  precision.
- Independent MPFR magnitude references over a component-precision matrix.
- Exponents as large and small as `1e+10000` and `1e-10000`, signed-zero axes,
  and fixed-precision scaling controls.
- Complex magnitude and logarithm accuracy using independent MPFR/MPC references.
- A tiny uncertain phase and the nonzero real part of `log(1+i*tiny)`, which must
  retain their actual uncertainty rather than being padded or rounded to zero.

On both base revisions, four of these eight tests fail before the fix.
The regression file is enabled with `float-mpfr` for its MPFR/MPC references.
The standalone Numerica patch passes **194 tests**: 172 library tests, 14 existing
API regressions, and the eight new regressions. The Astro/Malachite build also
passes `cargo check`. The Symbolica-tree version passes **196 tests** against
its bundled Numerica. Applying both standalone precision patches together passes
**207 tests**. Validation uses Rust 1.98.1 on x86-64 Linux.

## Applying

For standalone Numerica, apply `numerica-hypot-precision.patch` at the repository
root, then run:

```sh
cargo test --all-targets
cargo check --no-default-features --features integer-malachite,float-astro
```

For Symbolica main, apply `symbolica-main-numerica-hypot-precision.patch` at the
Symbolica repository root, then run:

```sh
cargo test --manifest-path lib/numerica/Cargo.toml --all-targets
```

This is an addendum to the earlier
[Numerica sqrt, phase and power patch](numerica-3.0.0-precision.md). Those repairs
are already present in Symbolica's current bundled Numerica. On standalone
Numerica 3.0.0, the earlier patch and this new patch can be applied together;
this patch does not duplicate or replace those earlier repairs.

SHA-256:

- Numerica format: `3b14b9983dd16784fbc999cf1318d2b24b09967c423b0d6d37d94dcf79c2a12f`
- Symbolica format: `608f27de405ca39075e4d025820fd069cc16161eac4d9d32fcc2b3d4c8f49243`
