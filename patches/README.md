# Local Symbolica and SymJIT fixes

`symbolica-dev.patch` is the complete local delta against Symbolica dev
`fb845d34bda8ccf1fedef6544d3aa46dc24944e3`. It includes the evaluator index repair,
early reuse of repeated non-inlined definitions before recursive compilation,
real-axis absolute values, polylogarithm convergence changes, a native binary64
dilogarithm specialization, and upstream-style
regression tests. It also adds portable JIT export/import with nested callback
definitions, the complex conjugation export alias, exact external-constant
resolution during direct JIT compilation, and registration of non-inlined user
functions even in the default Symbolica namespace. It has not been submitted or
accepted upstream.

The build script also resolves Git metadata paths with `git rev-parse --git-path`
before registering existing HEAD/ref paths with Cargo. This supports linked
worktrees, where `.git` is a file, and avoids watching a nonexistent `.git/HEAD`
that would trigger a full library rebuild on every invocation.

Review/apply it to a clean checkout of that revision with `git apply --check`
followed by `git apply`; do not apply it a second time to the already patched
sibling checkout. The OneLOop manifest currently selects that sibling checkout.
Library-manifest dependency overrides do not propagate to a consuming build root.

The numerical batch path additionally uses `symjit-2.24.1.patch` against the
**published SymJIT 2.24.1 crate snapshot**, in `../symjit-2.24.1-patched`. The
published `.cargo_vcs_info.json` records revision
`c9029c4f299556003edd655be04a548d72e80b1e`, which was unavailable from the public
GitHub repository when this patch was prepared. The release archive, not an
unrelated Git checkout, is the reproducible base:

- Archive: <https://static.crates.io/crates/symjit/symjit-2.24.1.crate>
- SHA-256: `9c1ead707a367a48f864f40b8c42a4882c0a9f0a2f183611e86a6e8cfe353da9`
- Local pristine baseline commit: `7bd06dd13ff8c25ecd8dbcaff133d29ce4f7b243`

The patch contains all nine modified source files, including their regressions,
and no build outputs. It fixes direct-translation power-function registration,
scalar complex callback SIMD lane packing, and propagation of nested conditional
bailouts to the existing scalar retry. A fourth fix makes complex IF truth exact
in both branch and join instructions: a condition is zero only when both its
real and imaginary components equal zero. It uses component comparisons, not a
norm or tolerance, and leaves general comparison semantics unchanged.
Two further repairs preserve the principal complex square-root lip in the
lower half-plane (including negative imaginary zero), and route integer powers
`-2`/`-3` through the existing reciprocal-power implementation instead of broken
specializations. Root signs are taken from the original imaginary component;
no numerical regulator or custom OneLOop function is introduced.
Two additional repairs map known-real absolute values to the correct component
registers and honor `fastmath=false` in complex multiplication/division and
generic root/reciprocal norm construction. Separate products must round before
addition: a one-sided FMA otherwise gives `z * conj(z)` a spurious nonzero
imaginary part, which can change an exact-axis logarithm or polylogarithm branch.
The `fastmath=true` paths and explicitly requested fused operations are unchanged.
A ninth repair skips the x86 SIMD input-transpose loop for zero-argument
expressions. Its former do-while loop dereferenced an empty input slice before
decrementing a zero counter. Both four- and eight-lane generators are repaired;
output transposition is retained. The AArch64 zero-count paths already use
bounded loops and require no change. This code-generation-only repair does not
change serialized scalar-family IR or require regenerating its assets.
The patch does not change OneLOop mathematics.
The original Cargo registry and public Git checkout remain untouched.

To prepare this sibling, run from the parent of `oneloopmaster`. The commands
refuse to reuse an existing archive or destination:

```sh
set -eu
test ! -e symjit-2.24.1.crate
test ! -e symjit-2.24.1
test ! -e symjit-2.24.1-patched
curl --fail --location https://static.crates.io/crates/symjit/symjit-2.24.1.crate --output symjit-2.24.1.crate
printf '%s  %s\n' 9c1ead707a367a48f864f40b8c42a4882c0a9f0a2f183611e86a6e8cfe353da9 symjit-2.24.1.crate | sha256sum --check
tar -xzf symjit-2.24.1.crate
mv symjit-2.24.1 symjit-2.24.1-patched
git -C symjit-2.24.1-patched apply --check ../oneloopmaster/patches/symjit-2.24.1.patch
git -C symjit-2.24.1-patched apply ../oneloopmaster/patches/symjit-2.24.1.patch
cargo test --manifest-path symjit-2.24.1-patched/Cargo.toml --lib
```

The consuming build root must include the following override (already present
in OneLOop's manifest):

```toml
[patch.crates-io]
symjit = { path = "../symjit-2.24.1-patched" }
```

All twelve standalone dependency tests pass on x86-64: real and complex powers
with portable save/load, two- and four-lane callback packing, and real/complex
batches through two nested conditional helpers. The additional complex-IF test
checks scalar/SIMD, strict/fast-complex, and original/restored evaluators with
exact truth assertions, signed zeros and smallest-subnormal conditions.
Two tests cover negative integer powers and principal complex roots
through all quadrants, signed axes and the origin, including scalar/SIMD,
strict/fast-complex, and portable restoration.
Three further tests cover known-real component absolute values, exact-axis
conjugate products, and separately rounded root/reciprocal norms. The latter
two include original/restored evaluators and scalar/SIMD batches with
`fast_complex=false`, matching the production setting.
The zero-argument regression covers real/complex constant outputs, scalar/SIMD,
AVX-512 requests off/on, original/restored evaluators, and batches 0/1/3/4/5/8/9.
This constant-only test does not establish AVX-512 callback support.
The integrated probe also passes all eight arithmetic, exact-truth,
complex-polylogarithm and conditional/helper cases
in scalar, SIMD-fallback and SIMD-branch-mask modes, with individual points,
batch sizes 3/4/5/8/12/16 and three outputs. Conditions cover both pure-imaginary
signs, both real signs, every nonzero complex quadrant and all signed zeros;
the exact-truth case also includes the smallest subnormal values and batch16.
It reports zero mismatches. The original seven-case probe, which used only real
conditions, reported 360 mismatches before the first SIMD repairs and zero
afterward; the expanded conditions additionally exposed the complex-IF defect.
After preparing both dependencies, reproduce the expanded probe with:

```sh
cargo run --manifest-path oneloopmaster/Cargo.toml --example simd_probe
```

A separate integrated root/power probe passes native primitive comparisons,
`sqrt(z) * z^(-1/2) = 1`, and the dB0 regression, with scalar/SIMD,
strict/fast-complex, original/restored evaluators and batches 1/3/4/5/8/12:

```sh
cargo run --manifest-path oneloopmaster/Cargo.toml --example power_branch_probe
```

The focused triangle branch probe also passes both formerly failing complex
C0 fixtures, intermediate sheet quantities, and portable restoration after the
strict-contraction repair. The previously spurious imaginary parts of conjugate
products and their square-root magnitudes are now exactly zero:

```sh
cargo run --manifest-path oneloopmaster/Cargo.toml --example triangle_branch_probe
```

The integrated debug `manual_evaluators` test passes after the zero-argument fix
(1/1, 74.30 s). It checks rejected dimension metadata, clone lifetime after
dropping the original, restored batches, and a zero-input `[2, pi]` expression
at batch size 5:

```sh
cargo test --manifest-path oneloopmaster/Cargo.toml --test manual_evaluators -- --test-threads=1
```

The checked production settings are O2, direct translation, `use_simd=true`,
`simd_branch=false`, `enable_simd512=false`, `fastmath=false`,
`fast_complex=false`, and `use_threads=false`. Four-lane AVX SIMD is enabled;
mixed conditional lanes use the existing scalar retry. The AArch64 two-lane
callback layout is unit-tested, and its corresponding bailout and fast complex
root generators are patched, but generated AArch64 machine code has not been
executed here.
AVX-512 remains disabled because the published callback ABI has only one
native-width trampoline, not a separate eight-lane callback. The AVX-512 bailout
generator is patched too; that alone does not establish eight-lane callback
support. This validation is focused, not an exhaustive backend test suite.
Run Symbolica runtime checks serially when using its restricted license.

The bundled SymJIT patch passes reverse `git apply --check` against the live
modified checkout and forward `git apply --check` against a fresh extraction of
the verified release archive.

For a fresh development clone, use this directory layout (run the commands from
the parent of `oneloopmaster`; `symbolica-dev-v3` must not already exist):

```sh
git clone https://github.com/symbolica-dev/symbolica.git symbolica-dev-v3
git -C symbolica-dev-v3 checkout --detach fb845d34bda8ccf1fedef6544d3aa46dc24944e3
git -C symbolica-dev-v3 apply --check ../oneloopmaster/patches/symbolica-dev.patch
git -C symbolica-dev-v3 apply ../oneloopmaster/patches/symbolica-dev.patch
cargo test --manifest-path oneloopmaster/Cargo.toml -- --test-threads=1
```

An already existing checkout must be inspected before changing its revision or
applying the patch. The bootstrap does not modify the original Fortran wrapper.
Cargo publishing is disabled until the dependency and permission questions are
resolved; the manifest's license identifier is provisional, not a permission grant.

The polylogarithm change uses the logarithmic expansion within its convergence
disk, ignores structural zeta zeros when checking convergence, and preserves
arbitrary-precision inputs in endpoint and real-axis predicates instead of
rounding them to `f64`. References:
[DLMF 25.12.12](https://dlmf.nist.gov/25.12#E12) and the unit-circle regression's
[DLMF 25.12.8](https://dlmf.nist.gov/25.12#E8).

The later binary64 specialization selects the exact integer order-2 tag once
when registering real/complex `f64` callbacks. Its finite-argument hot path uses
only binary64 arithmetic, with no per-call MPFR conversion or captured-order
clone. Inversion and reflection reduce the argument to `|u| <= 1`,
`Re(u) <= 1/2`; the Bernoulli series in `w = -log(1-u)` retains terms through
`B28`. Its omitted-series bound is below `1.5e-22` on that reduced domain. This
is not an arithmetic-error bound or a correct-rounding guarantee. Scaled
reciprocals/logarithms avoid squared-norm overflow, and `ln_1p` preserves tiny
arguments. Exact real `x > 1` retains Symbolica's negative imaginary value for
both signed zeros; nonzero imaginary inputs follow their actual side of the
cut. Nonfinite inputs, other tags/orders, and arbitrary-precision callbacks
retain the existing generic implementation. In particular, the specialization
does not change the earlier arbitrary-precision repairs described above.
The identities and coefficient construction follow
[DLMF 25.12.4](https://dlmf.nist.gov/25.12#E4),
[DLMF 25.12.6](https://dlmf.nist.gov/25.12#E6), and
[DLMF 24.2.1](https://dlmf.nist.gov/24.2#E1).

Three state-free numerical tests pass: exact-rational reconstruction of all
14 coefficients and the tail bound; 662 selected finite inputs at each of
128/256 reference bits (1,324 comparisons); and endpoints, cut sides,
callback/fallback behavior, plus 24 analytic tiny-width derivative checks.
The sampled inputs include signed axes, unit-circle/reduction boundaries,
subnormal widths, and scales from `1e-300` to finite binary64 extremes. The
largest observed `|error|/(1+|Li2|)` was `4.1499993958426e-16`; tiny arguments
and tiny-width components have separate relative checks allowing two final
subnormal ulps. These checks do not establish correctly rounded values or
exhaustive binary64 coverage. The tests call numerical functions directly,
without Symbolica state or license initialization. Their filter is
`dilog_f64_` in the Symbolica library test target, using the same patched
dependency graph and explicit feature configuration as this project.

With the new binary64 path, the twelve focused engine tests also pass, as does
the additional native special-function registration test (thirteen unique
state-using tests). These cover real callback registration and complex nested
runtime/fixed-Li2 JIT callbacks, including portable restoration. They are not
an independent high-precision numerical oracle; that role belongs to the
separate numerical checks above. Integral-level regression and timing results
must be checked separately before attributing an application-level speedup.

The refreshed six-file Symbolica patch has SHA-256
`1b8ab1c9a1f8a808a1a7c408187d1c5586181b4686e19d929d2683e2cde22c62`.
It exactly matches `git diff --binary` against the pinned `fb845d34` revision,
passes reverse `git apply --check` against the live checkout, and passes forward
`git apply --check` against a fresh archive of that exact base revision.

OneLOop's tests exercise the native index/absolute-value regressions and compare
324 sheet operations at machine and 128-bit precision. The full Symbolica test
suite has not been run as part of this change. An upstream-compatible revision
containing these fixes, or an explicitly agreed maintained patch, is still
required before packaging for distribution.

The follow-up reuse repair preserves current call arguments and captured-variable
bindings, reusing only the already registered function and its lifted constants.
Focused tests cover shadowed captures, distinct tags, constants mixed with runtime
polylogs, f64/Float remapping and a nested repeated-call graph. The latest pin
also strips unused captures; private helper identities retain the full original
capture scope while reusing the compact compiled body. Public master tags are
unchanged. Twelve unique focused engine tests pass with `--test-threads=1`,
including portable complex JIT restoration, conjugation, exact pi, and nested
function callbacks. Both strengthened portable tests also pass with a fixed-Li2
constant, actual executable dimensions, and the patched SymJIT dependency. This
is not the full upstream test suite. The bundled patch also passes reverse
`git apply --check` against the actual modified checkout.
