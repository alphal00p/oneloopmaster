# SymJIT 2.26.0 regression check

The six standalone sources from the existing 2.25.6 reproductions, unchanged,
with only the dependency upgraded to released 2.26.0. Run `python3 run.py`.
The runner now fails if any of the twelve subprocess cases fails and disables
core dumps. No Symbolica, OneLOop or custom numeric callbacks are involved.

All twelve cases pass on Linux x86-64 (2026-09-21). See `observed.txt` for inputs,
outputs, compiler and CPU. This covers complex root branches and signed-zero
lips, imaginary conditions in branch and join, strict conjugate multiplication,
direct fractional and negative integer powers, and zero-parameter SIMD batches.
It is not exhaustive SIMD, callback, nonfinite or cross-architecture coverage.

Two additional binaries reproduce remaining failures: `complex_sqrt_scales`
and `complex_reciprocal_scales`. Each accepts `generic` to select the generic
complex backend; the default uses packed complex arithmetic. They are excluded
from the twelve passing regressions in `run.py` and deliberately exit 101 on
2.26.0. See the [audit](../../patches/upstream-status-2026-09-21.md) and the four
saved scale-probe transcripts in this directory.

`complex_simd_callback` and `nested_simd_fallback` isolate two remaining SIMD
bugs. Run each with no arguments to reproduce failure (exit 101 on an AVX-capable
x86-64 host), and with `scalar` for a passing control. The default mode asserts
that a SIMD kernel was actually compiled. These are also separate from `run.py`.
