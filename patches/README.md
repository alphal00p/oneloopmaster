# Local Symbolica fixes

`symbolica-dev.patch` is the complete local delta against Symbolica dev
`0b57776bf911faeea7e28ea133706fb03740ffeb`. It includes the evaluator index repair,
real-axis absolute values, polylogarithm convergence changes, and upstream-style
regression tests. It has not been submitted or accepted upstream.

Review/apply it to a clean checkout of that revision with `git apply --check`
followed by `git apply`; do not apply it a second time to the already patched
sibling checkout. The OneLOop manifest currently selects that sibling checkout.
Library-manifest dependency overrides do not propagate to a consuming build root.

For a fresh development clone, use this directory layout (run the commands from
the parent of `oneloopmaster`; `symbolica-dev` must not already exist):

```sh
git clone https://github.com/symbolica-dev/symbolica.git symbolica-dev
git -C symbolica-dev checkout --detach 0b57776bf911faeea7e28ea133706fb03740ffeb
git -C symbolica-dev apply --check ../oneloopmaster/patches/symbolica-dev.patch
git -C symbolica-dev apply ../oneloopmaster/patches/symbolica-dev.patch
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

OneLOop's tests exercise the native index/absolute-value regressions and compare
324 sheet operations at machine and 128-bit precision. The full Symbolica test
suite has not been run as part of this change. An upstream-compatible revision
containing these fixes, or an explicitly agreed maintained patch, is still
required before packaging for distribution.
