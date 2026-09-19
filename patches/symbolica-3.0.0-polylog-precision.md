# Symbolica 3.0.0 polylog precision patch

[symbolica-3.0.0-polylog-precision.patch](symbolica-3.0.0-polylog-precision.patch)
applies at the Symbolica repository/crate root, against released 3.0.0 commit
`1f38ba9e0de48cb9f9e8c616b2ccd729fa45f606`. It changes the existing numerical
polylog implementation and adds regression tests. No new public dilogarithm
function or evaluator API is introduced.

SHA-256: `d73065ecd2420e56d37feb8e37b67e58ce0d992a807c41590f4a2c91a72dca70`.

## Behavior

The release rounds argument, order and convergence checks through binary64.
Consequently, `polylog(2, 1e-1000 + 2e-1000*i)` becomes zero, a tiny nonzero
imaginary order is treated as zero, and requesting thousands of bits does not
remove the series' 900-bit stopping cap. The logarithmic expansion can also
stop at a structural zero of zeta before its remaining terms converge.

The patch:

- Uses exact zero/integer/endpoint predicates and Float-valued radius and
  convergence comparisons.
- Scales stopping thresholds and checked iteration limits with working
  precision; exhaustion returns `None` instead of an unconverged partial sum.
- Includes geometric tail estimates and skips structural zeta zeros.
- Uses the logarithmic expansion on `|log(z)| < 4`, within its convergence
  disk `|log(z)| < 2*pi`, including the unit circle.
- Adds internal guard precision for integer-order summation and continuation,
  including the exponent gap between nonzero components. Results are rounded
  down to the requested precision, never padded after a computed precision loss.
- Evaluates the dilogarithm's logarithmic series with recurrent powers and
  reflected positive-integer zeta coefficients. MPFR's integer-order zeta
  routine supplies those coefficients and integer polylog endpoint constants.
- Preserves Symbolica's negative-imaginary convention on the exact real cut
  while distinguishing arbitrarily small nonzero imaginary parts on either side.

The formulas follow [DLMF 25.12.12](https://dlmf.nist.gov/25.12.E12), with
coefficient reflection from [DLMF 25.6](https://dlmf.nist.gov/25.6).
The continuation checks use [DLMF 25.12.4](https://dlmf.nist.gov/25.12.E4)
and the reflection identity in [DLMF 25.12.6](https://dlmf.nist.gov/25.12.E6).

General noninteger orders retain the existing `|z| < 0.95` numerical domain.
The patch does not provide a new arbitrary-order analytic continuation or
repair unrelated gamma, zeta, or Bessel algorithms. Convergence checks and
guard precision are not a correctly-rounded arithmetic guarantee.

## Tests

Nine added tests cover:

- Li2(i) against Catalan's constant, complex reflection and inversion, and
  precision refinement at 3456 bits (over 1000 decimal digits).
- Relative accuracy of tiny arguments and imaginary components near the cut,
  including signed-zero exact-axis controls.
- Near-endpoint real arguments against independent MPFR `li2` values.
- Exact order predicates and high-order endpoint constants.
- Reflected coefficients, a trilogarithm unit-circle identity, and a duplication
  identity for noninteger complex orders above the former 900-bit cap.
- Checked iteration limits and refusal to return an unconverged sum.
- The public `polylog(2, z)` expression-evaluation API at 3456 bits.

Three of the initial four numerical regression groups fail on the unmodified
release. The coefficient identity is a passing control. The tests use released
Numerica 3.0.0 and Graphica 3.0.0 from crates.io, with no dependency patches.
The separate Numerica precision patch is not required for these polylog checks.
All nine new tests pass, and the broader `transcendental::` run reports
**56 passed, 0 failed**. Existing optional GiNaC comparisons return early because
`ginsh` is unavailable; the new analytic-identity and MPFR checks all execute.

## Apply and reproduce

The published archive gives a build root without the Git repository's local
development overrides:

```sh
curl --fail --location https://static.crates.io/crates/symbolica/symbolica-3.0.0.crate \
  --output symbolica-3.0.0.crate
tar -xzf symbolica-3.0.0.crate
cd symbolica-3.0.0
git apply --check /path/to/oneloopmaster/patches/symbolica-3.0.0-polylog-precision.patch
git apply /path/to/oneloopmaster/patches/symbolica-3.0.0-polylog-precision.patch
cargo test --no-default-features \
  --features tracing_max_level_info,integer-gmp,float-mpfr,native_code_generation,bincode \
  --lib transcendental:: -- --test-threads=1
cargo check --no-default-features --features integer-malachite,float-astro
```

The new tests require `float-mpfr`. Native code generation and bincode are
enabled because the release's existing library tests also compile evaluator
tests. Forward application was checked against both the pristine release
checkout and the published crate. A freshly applied patch matches the tested
source files byte-for-byte and passes reverse `git apply --check`.

The Astro/Malachite compile check and formatting check also pass.
Validation environment: Rust 1.98.1, x86-64 Linux. OneLOop's dependency manifests
and existing evaluator assets are unchanged by this patch artifact.
