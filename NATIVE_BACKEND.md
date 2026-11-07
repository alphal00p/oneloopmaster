# Direct Rust backend

This lane adds ahead-of-time Rust functions for the same scalar masters. It is
separate from the existing native-expression/SymJIT backend. Direct Rust is now
the manual default and the sole numerical backend of bare master-Symbol hooks.
Manual `backend="symjit"` and `backend="expression"` remain available; native
compilation is not a promise of a speedup for every input or batch.

## One formula source

`examples/generate_native.rs` exports the exact repaired expressions and their
non-inlined definitions, then emits generic Rust functions under `src/native/`.
Generated functions use ordinary arithmetic, named locals, function calls and
lazy Rust `if` blocks. There is no instruction interpreter, expression builder,
function-map lookup, or JIT compiler in a warmed native scalar evaluation.

The generator tracks branch-local definitions and live outputs explicitly. It
shares repeated arithmetic and identical lazy conditional results only where
the earlier result dominates the later use. Conditional bodies are not hoisted
or evaluated speculatively. Large branches are factored into ordinary Rust
functions; there are no per-point hash tables or runtime expression caches.
Executable generator regressions cover captures, tuple ordering, branch-local
aliasing and an inactive division by zero.

Small exact sectors are written directly in [simple.rs](src/native/simple.rs):
scaleless cases, massless one/two-scale triangles, and the all-on-shell massless
box. They retain separate continued logarithms. For example, with
`L_s = log((-s-i0)/mu_squared)`, the normalized one-scale massless triangle is
`[L_s²/(2s), -L_s/s, 1/s]`; the massless on-shell box is
`[(2 L_s L_t - pi²)/(s t), -2(L_s+L_t)/(s t), 4/(s t)]`.
Exact zero dispatch is used, not a small-parameter cutoff.

## Hand-written binary64 sector kernels

For `f64`, the dominant massive sectors use hand-written kernels in
[`src/native/handwritten/`](src/native/handwritten/) instead of the generated
code: the three-mass triangle (`triangle.rs`'s `trif3` representation, nonzero
momenta with a nonnegative Källén function), and for three or four massive
propagators the contour box of `box_complex.rs` (OneLOop `boxc`) and the
four-mass Denner-Nierste-Scharf box (`boxf4`), with the same channel exchange
and contour selection as `box_series`. Each is a branch-for-branch transcription
of the expressions: the same exact zero tests (`x - |x|` sign tests included),
infinitesimal directions, root labels and continued sheets (`sheet.rs` mirrors
`sheet_exact.rs`). They evaluate one sector dispatch per call and never an
inactive branch, which makes them two to five times faster than the generated
code. Like the expression, the triangle takes the quadratic root that would
cancel in `(b ± d)/(2a)` from the product of the roots (OneLOop's `solabc`), so
that an exact zero produced by roundoff cannot reach a logarithm. Every other sector,
and every sector for `DoubleFloat` and `Float` (whose tracked precision follows
the expression's own grouping), uses the generated code. Regenerating the
expressions does not update these kernels: an expression change in these
sectors must be transcribed as well (as were the continuation fixes checked by
`tests/permutation_regressions.rs`: exact r-function roots at zero momenta and
thresholds, the roots, chart choice and root lips of `trif3`, the root lips of
`boxf4` (and its Denner-Dittmaier chart of `tests/d0_four_mass_orderings.rs`:
the relabelling that puts an invariant below its pseudo-threshold in the `p23`
slot, principal `x r13` and the explicit `eta(x, r13) Log R` term), and the
contour box's double-root split and real residue logarithms;
and the conditioning fixes of `tests/c0_coincident_roots.rs` and
`tests/d0_tiny_invariant.rs`: the term-wise dilogarithm divided difference of
`sheet_exact::near_difference`, `log(x)/(1-x)` next to one through
`log_one_plus_body`, and the contour box's exact coefficient sums; the
contour representation for every real-mass three-mass box it supports, of
`tests/d0_real_mass_orderings.rs`; and the contour box's Källén
discriminants, massless-end roots, exact line sums, expansion about
beta = 1 and far-root `r1` of `tests/d0_contour_thresholds.rs`; and the line
edges' pole offsets from Vieta's product of `tests/d0_vanishing_gram.rs`, the
product formed only where an offset cancels; and of
`tests/d0_contour_regressions.rs`, the product only for a line parallel to
the contour direction, the zero-external labelling with a massless corner at
the origin unless its chart has a double contour direction, the dilogarithms'
complements from the poles, formed only where a dilogarithm takes
`log(1 - x)` on a winding sheet or has an `x` that rounds to one, a far pole
2^52 beyond the near one where `K` rounds to zero, the other contour's chart
where the selected one's `a1` keeps fewer than 24 bits, and the tests that
return NaN where no chart is regular and the one used moves under its line's
or its `K`'s rounding, evaluated in a cold path only there and where `K`
rounds to zero).
`tests/native_handwritten.rs` checks them
against the generated expressions at 256 bits on randomized kinematics and on
every distinct C0/D0 call of the recorded Ninja benchmark streams
(`tests/data/ninja_master_rows.txt`); `src/native/handwritten/box_contour.rs`
checks each contour-box definition against its registered expression. Timings:
[performance/2026-10-17-claude-c0d0](performance/2026-10-17-claude-c0d0/README.md).

The massless finite triangle (sector 0 with three nonzero invariants, OneLOop
`trif0`) has a hand-written kernel too, for binary64, `DoubleFloat` and `Float`
(`triangle::massless`, the default wherever the hand-written kernels run and,
alone of the hand-written kernels, for `Float`). It
transcribes `triangle::triangle_finite_massless` with the conditioning fixes
of `tests/c0_massless_orders.rs`, which next to light-like legs and next to
nearly equal invariants made the value depend on the argument order: both
roots from Vieta's product, the dilogarithms' complements `1 + p X` from the
kinematics, the Källén function and the sums `b`, `u`, `v` with the small
invariant last (`massless_kallen`, `ordered_sum`), and the offset of the
roots' quotient from their gap `d/a`. Its residual error grows like `1/g` at a
relative gap `g` between the two largest invariants (a few `1e-15/g` in
binary64), as the generated code's did in the other orders.
`generated_c0.rs` and the prebuilt SymJIT asset `assets/evaluators/c0.bin`
were NOT regenerated with this expression (under the current Symbolica pin
the generator's output changes every C0 sector), contrary to
`assets/evaluators/README.md`'s rule to commit them together: until the next
full regeneration, `EvaluationBackend::SymJit` (binary64) and
`set_handwritten(false)` still run the direct forms of sector 0. Next to
light-like legs those are off by up to the whole value in binary64, SymJIT
included; native `Float` ran them too until lit-olm-ap, and there they kept up
to 60 bits fewer than they tracked in the orders with a light-like leg first
(at 64 bits up to the whole value with every bit tracked), so native `Float`
takes the hand-written kernel since. Its residual next to nearly equal
invariants is about `1/g²` at 128 and 256 bits (28 bits at `g = 1e-4`), as in
the interpreted expressions, and the tracked precision sees most of it (up to
9 bits beyond numerica's 8-bit slack in `tests/c0_massless_orders.rs`, at most
`1/(256 g²)`; the generated direct forms missed 36 to 40 bits there). The interpreted expressions
(`EvaluationBackend::Expression`, binary64 and `Float`) and evaluators built
from the expressions at run time have the fix.

## C0 at a vanishing Källén function (`native/gram.rs`)

Where the Källén function `λ(p1, p2, p3)` of a triangle's invariants (a
multiple of the Gram determinant) vanishes, the roots of every sector's
quadratic coincide. The one-mass sector forms its discriminant
`b² - 4ac = λ/(|m|² m)` from rounded coefficients (the two- and three-mass
sectors take it from `λ`), so there its roots split by rounding noise and its
dilogarithm divided differences are formed from that noise. Binary64 and
`DoubleFloat` stay near their precision through the near-coincidence
expansion of `sheet_exact::near_difference` (`DoubleFloat` keeps only binary64
accuracy, about 1e-16), but native `Float`, whose wider precision keeps the
direct differences, was wrong by up to the whole value at every precision:
`C0(0, p, p; m, 0, 0)` with `p = -0.2128`, `m = 0.3607` gave -5.871 with all
512 bits tracked for -4.928 (numerica keeps the tracked precision of a
computed factor of exactly `±1`, here a `log(r)/(1 - r)` of two noisy roots,
and of `log` next to one). Zero-invariant two-mass triangles were wrong by up
to 1e154 in the orders with the massless line between the two masses (on
every route), and the massless triangle was NaN in two of its six orders.

Native `Float` therefore reduces such a triangle exactly to bubbles. With the
modified Cayley matrix `Y_ij = m_i + m_j - p_ij`, any `b` with `Y b = e` gives
`C0 = -Σ_i b_i B0^(i) + (N - d - 1) (Σ_i b_i) C0^(d+2)` (`B0^(i)` without
propagator `i`), and `Σ_i b_i = e^T Y^-1 e` is proportional to the Gram
determinant, so at `λ = 0` the triangle is `-Σ_i b_i B0^(i)`, in every mass
sector and argument order. Two cases are reduced, decided by exact predicates
on the represented inputs (`NativeFloat::exact_rational`, so not by a rounded
`det Y`): `det Y ≠ 0`, `b = adj(Y) e / det Y` (an infrared-divergent triangle
has `det Y = 0`, so this case never meets a scaleless bubble); and three
vanishing invariants with at most one massless line, `b = (e_a - e_c)/(m_a -
m_c)` for the two masses that differ most. Otherwise (`det Y = 0` at a nonzero
invariant, e.g. the equal on-shell masses of STATUS.md's open C0 limitations)
the sector code runs. Binary64 and `DoubleFloat` keep their kernels; they
reach their precision with `set_gram_reduction(true)` (validation). Next to
`λ = 0` the one-mass sector code still loses about `log2(s²/|λ|)` bits (native
`Float` tracks them; the interpreted expressions not always, STATUS.md).

The expressions take the same reduction (`triangle::with_gram_reduction` in
`triangle_series`), decided by exact zero tests in the evaluator's own
arithmetic like the sectors' Gram tests; at three vanishing invariants they
weight the bubbles symmetrically, `b_i = (3 m_i - Σ m)/(3 Σ m² - (Σ m)²)`. The
bubbles are calls of the registered B0 master, so the shared definitions gain
no function or symbol: registering new ones moved an interpreted D0 at 128
bits by 1e-30 (its evaluator's instruction order follows the definitions).
This serves the interpreted expressions (`EvaluationBackend::Expression`, in
binary64 and `Float`) and every evaluator built from them at run time, among
them NinjaRS's `SymbolicaMasters<Float>` (its arbitrary-precision scalar
masters). The generated code and the SymJIT asset were not regenerated:
regenerating them would give binary64 and `DoubleFloat` the reduction too.
`tests/c0_gram_singular.rs` checks both routes against closed forms and
mpmath Feynman-parameter integrals in every order (native `Float` at 128 to
1024 bits, the expressions at 128 and 256).

This keeps the expression implementation as the mathematical source. The older
dependency-free numerical port is not copied into production: its finite-width
regulators and then-unrepaired sectors would undo repairs already made here.
The original Fortran wrapper and that separate numerical port remain untouched.

Exact rational constants, pi and fixed-argument dilogarithms are described in
generated initialization data and converted once at working precision.
Symbolica's export passes the constants a shared body uses as arguments,
through every enclosing body; the emitter reads a helper parameter that is the
same constant at every call from the context instead, and leaves it out of the
helper's signature. The Symbolica instruction export includes symbolic
constant metadata: its numeric zero placeholders must never be mistaken for
exact zeros by another backend.

## Numeric domains

The generated functions are generic over `NativeFloat`, implemented for `f64`,
Symbolica `DoubleFloat`, and arbitrary-precision `Float`. The core still has only
Symbolica as a direct dependency. Complex masses use `Complex<T>`.

The native primitives preserve exact complex zero tests and principal square-root
and logarithm lips. For binary64 off the axes and inside `[2^-500, 2^500]`, complex
reciprocals use one division and norms `sqrt(x²+y²)` (within a few ulps; tested
against 256 bits); axes, non-finite values and extreme scales keep the scaled
algorithms. They avoid converting arbitrary-precision predicates through
f64. Dilogarithms call the same state-free native Symbolica numerical kernels:
binary64 calls `symbolica::transcendental::dilog_complex_f64` directly (a
Bernoulli-series kernel with Symbolica's lips and signed zeros, Symbolica PR89),
without the registered-callback lookup or MPFR. Its normwise relative error is at
most 1.2e-15 (about 10 units of 2^-53, reached next to Re z = 1/2 where the
reflection cancels; 99.86% of points within 4 units) on a 200,745-point review
grid against mpmath, where the previous MPFR route was correctly rounded. DoubleFloat's
dilogarithm uses a precision-preserving MPFR conversion of both
compensated components, while ordinary generated arithmetic remains DoubleFloat.
This is not a runtime scalar-master evaluator callback.

## Tensor coefficients

[`tensor.rs`](src/native/tensor.rs) is a hand-written generic port of OneLOop's
bubble tensor routines (`bub0`-`bub1111`, `bnlog`, `solabc`) and `tadpn`, exposed
as `NativeTensorEvaluator<T>` with `TensorIntegral::{An, B11, B111, B1111}`. It is not
generated from the symbolic expressions: their Passarino-Veltman forms cancel
like `(m²/p²)^rank` at small momenta, while OneLOop's root-pair representation
with the logarithmic `bnlog` series stays accurate there. It keeps OneLOop's
mass ordering and reconstruction, replaces its `10*EPSN` degeneracy tolerances by
exact zero tests and its `IEPS` regulator by explicit lips, and shares
`log(1-1/x)` between ranks. A massless propagator's unit root `x = 1` is formed
exactly (`(p²-m²)/p²` is divided, not scaled by `1/p²`). `f64` evaluators run
a binary64 specialization of the same algorithm
([`tensor/binary64.rs`](src/native/tensor/binary64.rs)): real kinematics in
real arithmetic, one `bnlog` (`2 Re`) for the complex-conjugate roots of real
masses below threshold, and the `bnlog` series by Estrin's scheme on a fixed
8, 16 or 24 coefficients instead of Horner's rule; `DoubleFloat` and `Float`
keep the generic code. `tests/native_tensor_binary64.rs` compares it with the
generic binary64 code and 320 bits on 8,000 random rows over every sector, with
the symbolic `b0`-`b1111`/`an` expressions at 256 bits, and with OneLOop on every
bubble/tadpole call of Ninja's recorded streams. Besides OneLOop's
values and Ninja's B111 formulas, `tests/native_tensor.rs` checks all five
bubble coefficients against direct Feynman-parameter integrals (728 rows,
`tests/generate_tensor_feynman_fixtures.py`, mpmath): 320-bit results agree
to 3e-27. Validation and timings are in
[the 2026-10-17 report](performance/2026-10-17-claude/README.md).

### Rank four (`B1111`)

`TensorIntegral::B1111` returns `[B0, B1, B11, B111, B1111]` for reductions with
numerator rank `n+2`, with `B1111 = 1/(5 eps) - int_0^1 x^4 log(D(x)/mu²) dx`.
It is `bub1111` with one deliberate difference: after the mass swap (taken for
`p² = 0` with `|m1| >= |m0|`, and for `p² != 0` with `|m1| < |m0|`, including the
massless branch) the result is `B_n = (-1)^n sum_k C(n,k) B_k'`, i.e.
`B1111 = B0' + 4 B1' + 6 B11' + 4 B111' + B1111'`. OneLOop 3.6 and 3.4 negate
this sum in all three swap branches (`avh_olo.f90` of OneLOop-3.6 in
MadGraph5_aMC@NLO's CutTools, lines 3400, 3424 and 3439; `avh_olo_bub.h90` of
OneLOop-3.4 in IREGI, lines 89, 115 and 132), so their `B1111` has the wrong
sign wherever the swap is taken, not only at zero momentum (e.g. the corrected
`B1111(0; 2, 5) = -0.2998475029970859`, where OneLOop returns `+0.2998...`).
Neither OneLOop nor a library built on it (such as Ninja's `avholo`) is a valid
`B1111` reference;
the tests use Feynman-parameter integrals, the exact `b1111` expression and the
swap identity instead.

At 53 bits rank four needs 9 series coefficients at `|y|_1 = 0.01` and 17 at
`0.2`, one more than the binary64 tables hold. Rank four therefore keeps the 8,
16 and 24 coefficients with its own radii, `0.0085` and `0.17` (8 and 16
coefficients reach `0.0090` and `0.1736`; all 24 reach `0.534`), checked by
`Binary64::new`. The rank-four length is chosen apart from the others only in
the `B1111` instantiation (a compile-time condition), so ranks zero to three are
evaluated exactly as before: `B11` and `B111` are bitwise unchanged, and the
first four coefficients of `B1111` are bitwise those of `B111` on every route
(binary64 kernels, generic binary64 code, `DoubleFloat`, `Float`).

`tests/native_tensor_b1111.rs` checks every route against the exact
`b0`-`b1111` expressions (2,576 rows; 320 bits, or 1280/5120 where the
Passarino-Veltman form cancels beyond 320; Float(256) within 2e-61, DoubleFloat
2e-29, binary64 within 5e-14 or the generic code's own accuracy), against the
Feynman-parameter integrals per sector, on a dyadic ladder `p² = ±2^-4..2^-24`,
and the swap identity in every branch and both mass orders (binary64 4.4e-15,
DoubleFloat 8e-26 at `|p²| = 1e300` where its logarithm keeps about 90 bits,
320 bits 1.5e-83, relative to the size of the terms). Binary64 errors of
`B1111` are measured relative to the largest coefficient of the row: the swap
reconstruction cancels to a value up to about ten times smaller than its terms.

The exact `b1111` expression is a Passarino-Veltman recursion that divides by
`2p²` at each rank above one and has no small-momentum branch. In fixed
precision it loses about `4 log10(m²/p²)` digits: in binary64 at masses
`(1, 1.5)` its relative error is 4e-6 at `p² = 2^-8`, 0.6 at `2^-12` and 1e15
at `2^-24` (also the `1/eps` coefficient, which it obtains through the same
recursion); in `DoubleFloat` it is 4e-10 at `2^-16` and 2.9 at `2^-24`
(masses of order one). The native family keeps
2e-15 relative to `max_k |B_k|` over the whole ladder. Fixed-precision callers
(`f64`, `DoubleFloat`, binary64 fast paths) should take bubble tensors from
`NativeTensorEvaluator`; `Float` can use either at a working precision that
covers the cancellation.

## Verification and performance gates

The completed gates compare all three Laurent coefficients independently:

- Original acceptance fixtures for all five families, scalar and mixed batches
  including partial tails; explicit diagnostics for known fixed-f64 failures.
- DoubleFloat against original controls and higher-precision primitives.
- Genuine complex C0/D0 at 1000 digits against higher-precision expression
  evaluation, plus independent analytic limits and tiny-width controls.
- Public Symbol hooks, Python Decimal boundaries, and explicit manual SymJIT
  evaluation checked as separate routes.
- Matched warmed Rust/native, Rust/SymJIT, hook, batched Python, and unchanged
  Fortran measurements using the existing survey protocol. The target remains
  at most 1.5 times Fortran; measured misses must be reported explicitly.

The prior binary64 report remains a historical measurement of the SymJIT
implementation at its recorded commit. It is not evidence for this new backend.

The final release suite passes **114 Rust tests**; the standalone extension
passes **27 Python tests**. The expanded 342-point table passes at native Float
256 bits, while the explicit binary64 diagnostic retains 14 B0/dB0 failures.
See the [release audit](NATIVE_RELEASE_AUDIT.md),
[precision audit](ARBITRARY_PRECISION_AUDIT.md), and
[matched timing report](performance/2026-09-08-native/README.md) for the full
inventory, retained failures and per-workload performance limits.

### Baseline correctness checks (2026-09-08)

The first generated implementation passed the 293-point acceptance set for
all five families with `f64` (scalar and mixed batches 1/3/4/5/31/256/1024) and
`DoubleFloat`. A separate test passed genuine complex inputs in all five
families at 3456 working bits against 1120-decimal-digit expression evaluation,
requiring relative error below `1e-1000` independently for real and imaginary
components. These are sampled checks, not exhaustive analytic-region coverage.

The unoptimized test times were 0.76 seconds, 2.44 seconds and 152.32 seconds
respectively, including the tests' different setup costs. They are **not**
matched Fortran performance measurements. The final factored implementation's
matched measurements are recorded separately in the report linked above.
