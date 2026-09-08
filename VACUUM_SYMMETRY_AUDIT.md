# Exceptional vacuum triangle: analytic symmetry versus the reference

The intended result is the analytic, permutation-symmetric vacuum limit. We do
not copy inconsistent exceptional-point outputs from the reference Fortran.
This is a documented exception to bitwise or numerical reference parity, not a
relaxed comparison tolerance.

## Reproduction

The inspected reference is the unchanged OneLOop source bundled in
`SecretGmG/OneLOopBridge` revision
`402bf0d70bfa2a9d8ae66b8251cc4395b7b5ab99`. Its combined Fortran source has a
2024-11-23 banner. At `p1²=p2²=p3²=0`, `mu²=1`, the original DP library returned:

| Ordered squared masses | Finite C0 |
| --- | --- |
| `[-1, 2, 3]` | `-0.361861096127785498 - 0.261799387799149741 i` |
| `[2, -1, 3]` | `-0.361861096127785387 + 0.785398163397448390 i` |

Thus the specific disputed permutation is **213**, exchanging the first two
propagator masses while all three external invariants stay zero. This is not a
claim that every invocation of a particular permutation routine is faulty.

There is also an unambiguous positive-mass example, without negative-axis branch
choices. All six orders of `[1,2,3]` should give
`2 log(2) - 3 log(3)/2 = -0.2616240718822739...`. The bundled DP library gave:

| Permutation | Squared masses | Finite C0 |
| --- | --- | --- |
| 123 | `[1,2,3]` | `-0.261624071882274156` |
| 132 | `[1,3,2]` | `-0.214288222551819341` |
| 213 | `[2,1,3]` | `-0.261624071882273879` |
| 231 | `[2,3,1]` | `+1.78726335506856970` |
| 312 | `[3,1,2]` | `-0.26162407188227349` |
| 321 | `[3,2,1]` | `+0.163299316185545218 + 0.463477515359967507 i` |

A fresh `-O0` compilation of the original combined source reproduced the DP
results, excluding a stale archive and an optimization-only explanation. QP
fixes the specific negative-mass 213 discrepancy, but not all other exceptional
orders. At nonzero Euclidean invariants `[-1,-1,-1]`, all six positive-mass
permutations instead agree to about `1e-15`.

## What the source supports

The bundled `src/avh_olo_c0.h90` sends negative Källén discriminants to
`trif3HV`, but zero to `trif3`. At the exact vacuum point, the latter's selected
root chart can have all three quadratic coefficients zero. Floating residuals
then determine a mathematically indeterminate root calculation. The alternate
`trif3HV` implementation contains a vacuum limit and agrees with the independent
partial-fraction expression in the examined cases.

The evidence therefore supports an **exceptional zero-momentum / degenerate-chart
handling limitation in this bundled source**, not a general permutation bug or
a claim about an untested current upstream release. The original paper itself
discusses limitations at exceptional points; see
[OneLOop, arXiv:1007.4716](https://arxiv.org/abs/1007.4716).

## Independent target

For distinct squared masses `a,b,c`, partial fractions give

```text
C0(0,0,0;a,b,c) = - sum_i [ m_i² log(m_i² - i0)
                            / product_{j != i}(m_i² - m_j²) ].
```

Here `m_i²` denotes each complete squared-mass argument, not its real part.
For finite negative widths use the analytic logarithm in the lower half-plane.
For `[-1,2,3]` the result is `2 log(2)/3 - 3 log(3)/4 - i*pi/12`.
The scale cancels and both Laurent poles vanish. The implementation preserves
the inherited logarithmic sheets of mass **ratios**, rather than replacing them
with unrelated principal logarithms. Repeated masses use separate exact limits.

`tests/vacuum_symmetry.rs` is the regression against this independent formula:
all six mass orders, three sign configurations, zero and two nonzero widths,
three scales, and 128/256-bit evaluation (324 evaluations, each checking all
three coefficients). It passes against the patched `fb845d34` Symbolica engine
on 2026-09-08. This focused result does not certify other exceptional surfaces.

## Scope of the broader audit

The two main tables contain **613 unique ordered scalar inputs**, not 635
independent points: A0 15, B0 56, dB0 56, C0 218, D0 268. Tests of hooks, maps,
other precisions, and batch layouts often reuse these same inputs. They check
different execution paths but do not establish additional analytic-region
coverage.

The previous audit was **not exhaustive**. In particular there is no complete
analytic branch/leaf inventory; mass-mask coverage is not equivalent to region
coverage. Important remaining dimensions include all consistent momentum/mass
permutations, independent width limits, broad threshold and Gram/Cayley
intersections, large mass hierarchies, and all channel-zero patterns. A
"physical" random-input label in older diagnostics denotes a sign pattern, not
validated external four-momenta. High-precision convergence alone also does not
prove that an analytic continuation is on the correct sheet.
