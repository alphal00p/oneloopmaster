# Symbolica / SymJIT audit, 2026-09-21

Target: unmodified Symbolica main `a19c760dd567c239f30d87e4e924ca2f8b8457ab`,
with bundled Numerica and released SymJIT 2.26.0 (release VCS revision
`530304a07d1be6d5abc80291aa5e92cf28ac5546`). Tests run on Linux x86-64.
No dependency source patches are applied.

## Resolved during this update

The first fetched main, `3f588e06`, failed to compile because
`src/transcendental.rs` declared an absent `numeric` module. The force-pushed
main `a19c760d` removes that missing-module dependency and compiles.

Symbolica's new shared non-inlined function bodies remove the need for OneLOop's
explicit-pi wrappers and common argument slots. Native instruction exports now
use `Arc` for shared bodies; the generator consumes the public constant metadata
instead of lifting pi into a synthetic argument.

All twelve standalone reproducer/control invocations for the six previous
SymJIT defects pass on 2.26.0: principal roots and signed lips, imaginary
conditions in both branch and join, strict conjugate products, direct
fractional powers, direct negative powers, and zero-input SIMD batches.
Sources and the complete transcript are in
[`reproducers/symjit-2.26.0`](../reproducers/symjit-2.26.0/README.md).
These results supersede the previous 2.25.6 status, not the historical evidence.

## Remaining: SymJIT complex square-root scaling

Both `fast_complex=false` and `true` fail five of eight elementary finite-input
cases, with O2, indirect translation, scalar evaluation and `fastmath=false`:

| Input | Expected principal root | Observed |
| --- | --- | --- |
| `1e300 + 0i` | `1e150 + 0i` | `inf + 0i` |
| `-1e300 - 0i` | `0 - 1e150 i` | `0 - inf i` |
| `1e-300 + 0i` | `1e-150 + 0i` | `7.071067811865476e-151 + 0i` |
| `-2e300 i` | `1e150 - 1e150 i` | one infinite component and one zero component |
| `-2e-300 i` | `1e-150 - 1e-150 i` | one infinite component and one zero component |

The root and its components are representable. Zero and small-width
`1 +/- 1e-20 i` controls pass. The MRE emits only `symbolica_sqrt`; there is no
Symbolica, OneLOop, custom callback or unsafe consumer code. Componentwise
relative comparisons prevent a tiny squared error from underflowing to zero.
The backend still forms an unscaled squared norm. It needs a scaled magnitude
and root construction, preserving the repaired signed-zero lips.

```sh
cd reproducers/symjit-2.26.0
cargo run --locked --bin complex_sqrt_scales
cargo run --locked --bin complex_sqrt_scales -- generic
```

Both commands are expected to exit 101 on this release. Complete outputs:
[specialized backend](../reproducers/symjit-2.26.0/sqrt-scales-packed.txt),
[generic backend](../reproducers/symjit-2.26.0/sqrt-scales-generic.txt).
OneLOop retains its numerical square-root callback.

## Remaining: SymJIT complex reciprocal scaling

Both complex backends fail six of seven finite-input reciprocal cases under
the same settings. `1/(1e300+0i)` becomes zero instead of `1e-300`;
`1/(1e-300+0i)` becomes `inf + NaN i` instead of `1e300 + 0i`.
Pure-imaginary and equal-component inputs fail as well. The ordinary `1/2`
control passes. The MRE emits one integer power with exponent `-1`.

The generic `complexify::recip` and x86 packed `amd::complex::recip` paths
form `re*re + im*im`, overflowing or underflowing before division. These paths
need scaled complex division/reciprocal arithmetic. This is separate from the
now-fixed sign of the direct `-2`/`-3` specializations.

```sh
cargo run --locked --bin complex_reciprocal_scales
cargo run --locked --bin complex_reciprocal_scales -- generic
```

Both commands are expected to exit 101. Complete outputs:
[specialized backend](../reproducers/symjit-2.26.0/reciprocal-scales-packed.txt),
[generic backend](../reproducers/symjit-2.26.0/reciprocal-scales-generic.txt).
OneLOop's JIT adapter does not replace general reciprocal arithmetic, so the
passing ordinary integral fixtures do not establish correctness at these scales.
The default native backend and explicit arbitrary-precision backend are distinct
paths; this audit does not claim that all extreme integral inputs are supported.

## Remaining: scalar complex callbacks corrupt SIMD lanes

A standalone `Defuns::add_sliced_func` callback implements
`(2-3i)*a + (-5+7i)*b`. All individual scalar calls are correct. A four-row
SIMD batch produces different numbers; e.g. the first result is `-251+38i`
instead of `-287+11i`. Disabling SIMD makes the same batch correct.

`Defuns::add_sliced_func` still selects the generic scalar trampoline for
`Complex<NativeSimd>`. That trampoline treats interleaved scalar complex pairs
as if they had the native SIMD layout (all real lanes followed by imaginary
lanes). The earlier 2.24.1 patch contains a component-aware trampoline, but
this release does not contain that repair.

```sh
cargo run --locked --bin complex_simd_callback          # fails, exit 101
cargo run --locked --bin complex_simd_callback -- scalar # passes
```

[Failure transcript](../reproducers/symjit-2.26.0/complex-simd-callback.txt),
[control](../reproducers/symjit-2.26.0/complex-simd-callback-control.txt).
No Symbolica or special-function implementation is involved.

## Remaining: nested SIMD conditional fallback is lost

A compiled child computes `if(condition, z, 2*z)`. Its parent returns
`child(condition,z) + z`. All scalar calls pass. For a four-row batch with
mixed true/false lanes and `simd_branch=false`, evaluation returns only `z`,
rather than `2*z` or `3*z`. With SIMD disabled, all four rows pass.

This isolates the propagation of an incomplete child kernel to the existing
scalar retry. The prior 2.24.1 patch addressed that path; it remains broken in
2.26.0. The MRE uses `Defuns::add_applet`, with no scalar user callback.

```sh
cargo run --locked --bin nested_simd_fallback          # fails, exit 101
cargo run --locked --bin nested_simd_fallback -- scalar # passes
```

[Failure transcript](../reproducers/symjit-2.26.0/nested-simd-fallback.txt),
[control](../reproducers/symjit-2.26.0/nested-simd-fallback-control.txt).
Both SIMD MREs require an actual compiled SIMD kernel and fail explicitly when
one is unavailable, rather than silently testing only scalar execution.

## Integrated evidence and scope

`cargo run --release --example simd_probe` reports **603 mismatches** across
callback and nested conditional cases. All eight scalar cases and all
arithmetic, exact-truth and top-level-conditional modes pass. The isolated
SymJIT-only reproducers above establish that these are backend defects, not
OneLOop expressions or Symbolica's polylogarithm implementation. Production
OneLOop continues to disable SIMD; the diagnostic deliberately enables it.

The ten [standalone Symbolica evaluator/startup regressions](../reproducers/symbolica-main/README.md)
pass on the final main revision. Function bodies now intentionally see globals
and their own parameters even when inlined; caller-local expansion is available
through `FunctionMap::add_aliases`. The former OneLOop capture test has been
migrated to that documented contract. No additional Symbolica correctness
defect was established by these focused checks.

Follow-up on 2026-09-22: the full community integration exposed an
[allocator-path crash](community-host-allocator-2026-09-22.md) during Python
thread transitions. The local host now selects the system allocator; this
separate finding does not invalidate the evaluator regression results above.

The exact-root callback also passes a new OneLOop regression at scales `1e300`
and `1e-300`, including a negative imaginary axis. The root/power and complex
triangle probes pass, including portable restoration. Eight native code-generator
tests pass, and all five native families export successfully with registered
constants recovered from upstream metadata. Newly exported Rust files were
written to a temporary directory; the previously validated native formula files
in the repository were not replaced.
