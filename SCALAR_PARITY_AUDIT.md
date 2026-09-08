# Scalar-integral parity audit

## Native/ArbPrec follow-up — 2026-09-08

The [new native release audit](NATIVE_RELEASE_AUDIT.md) records 114 passing Rust
tests, 27 passing standalone Python tests, native-only master hooks, arbitrary
precision and complete-expression inspection. The direct native backend passes
the 342-point table at 256 bits; its separately selected binary64 stress audit
still fails on 14 finite coefficients (six B0, eight dB0), with unchanged
tolerances. Global scalar parity remains unestablished. The following section
preserves the earlier SymJIT-default snapshot and its historical test counts.

## Preserved SymJIT-default release validation — 2026-09-08

**Full scalar parity is not established, and the requested performance target is
not met.** The current release checks use Symbolica dev
`fb845d34bda8ccf1fedef6544d3aa46dc24944e3` plus the recorded Symbolica/SymJIT
patches, including native binary64 Li2. The original Fortran wrapper and the
separate dependency-free numerical Rust core remain unchanged.

- The default release suite passes **57 tests**, with zero failures and four
  ignored test entries (three opt-in diagnostics and the initialization helper).
  This includes public master hooks/native maps, portable restored evaluators,
  mixed batches and tails, manual evaluator shape/clone checks, and fresh-process
  eager initialization. The ignored entries are not all unresolved failures:
  one is the child-process helper executed by the initialization parent test.
- The added bare-master outer-JIT regression passes without a native function
  map, checking all three coefficients against fixtures and prepared backends.
  Original/restored evaluators cover scalar and 1/3/4/5/31-row batches, repeated
  inputs and tails. All 293 acceptance rows use canonical tag order; the other
  five orders use eight controls per family. This does not constitute testing
  every accepted row in every tag order or adding new analytic regions.
- The explicitly selected expanded table passes at **256 bits: 342 points /
  1,026 coefficients**, with unchanged dimension-normalized tolerances. The
  current fixed-f64 rerun still fails on **14 finite coefficients: six B0 and
  eight dB0 small-momentum cases**. Higher-precision success is not automatic
  recovery by the fixed-f64/SymJIT API. The retained failing inputs and earlier
  convergence analysis are preserved below.
- The repeated-root real-box regression is now enabled and passes. Its exact
  `D0([0,0,0,0,4,-3],[1,1,1,1])` target is checked at both 128 and 256 bits.
  The broader threshold test retains one explicit limitation: at
  `s=2, t=-4, width=1e-30`, the 128-bit sample checks finiteness only; its
  independent numerical target is asserted at 256 bits. Passing that test does
  not certify target accuracy for every sampled width at both precisions.
- The corrected vacuum triangle passes **324 evaluations** against the
  independent symmetric lower-lip formula: six mass orders, three real-sign
  sets, three widths, three scales, and 128/256 bits. For squared masses
  `[-1,2,3]`, the target has imaginary part `-pi/12` in every order. The user
  has selected analytic symmetry over inconsistent exceptional-point legacy
  values. See [VACUUM_SYMMETRY_AUDIT.md](VACUUM_SYMMETRY_AUDIT.md) for exact
  permutations, original-source evidence, and the scope of that conclusion.

The 293 acceptance rows and 342 expanded rows represent **613 unique ordered
inputs**, not an exhaustive inventory of analytic regions. Repeating them at
other precisions or batch layouts validates execution paths, not additional
kinematic regions. Arbitrary threshold/Gram/Cayley intersections, input
hierarchies, limiting widths, and the full original domain remain open coverage
work. The dependency's full upstream test suite and runtime on other CPU
architectures have not been established by these checks. The most recent Python
and benchmark status is recorded separately in [STATUS.md](STATUS.md) and the
[performance report](performance/2026-09-08/README.md); this native release run
must not be presented as a Python/community-host runtime test.

## Preserved prior-engine repair follow-up — 2026-09-08

The following follow-up records the earlier patched `0b57776b` engine. Its
44-test count, ignored repeated-root failure, unresolved vacuum sheet, and
pending analytic-versus-legacy decision describe that historical stage; the
current validation above supersedes those status statements. The detailed
counterexamples and independent targets remain useful evidence and are retained.

Full parity is still **not established**. Three agents are investigating branch
conventions/master hooks, complex boxes, and original-source coverage, with
cross-review of the proposed formulas. The original wrapper remains unchanged.

Verified during this follow-up:

- The complete default suite passes **44 tests**, with four opt-in diagnostics
  ignored. One of those diagnostics is the confirmed repeated-root box failure
  described below; a passing default suite does not hide or resolve that gap.
- The complete expanded corpus now passes at **256 bits: 342 accepted Fortran
  points / 1,026 coefficients, zero mismatches**, with unchanged comparison
  tolerances. Family counts are `[13,49,50,100,130]` for A0/B0/dB0/C0/D0.
  This is a chosen-precision evaluation of every point, not automatic recovery
  by a fixed-f64 evaluator and not proof of unrestricted domain parity.
- The final fixed-f64 rerun completed with **14 failures / 1,026 coefficients**
  (exit 101, 127.13 seconds), down from the baseline's 44. All remaining failures
  are finite small-momentum bubble coefficients: six B0 and eight dB0. Their
  fixture lines are `21,23,25,27,29,31,33,53,55,57,59,61,63,65` in
  `tests/data/scalar_audit.txt`. All sampled A0/C0/D0 coefficients and all pole
  coefficients pass. The 256-/384-bit bubble convergence checks below establish
  recoverable fixed-precision cancellation for these retained points; the test
  still fails at f64 and no tolerance or expected value was changed.
- The five tests in `tests/branch_limits.rs` pass. All 49 B0 and 50 dB0 points
  from the 342-point corpus match Fortran at **both 256 and 384 bits**, using
  unchanged dimension-normalized tolerances. The two precisions agree to a
  normalized absolute difference below 1e-32. This includes small-momentum
  cancellation, not a claim that fixed f64 evaluation is stable there.
- A native real-axis grouping identity `if(z,z,0)` prevents conjugation from
  distributing a composite argument into differently rounded subexpressions.
  Without it, one timelike B0 passed at 256 bits but selected the wrong lip at
  384 bits because its nominally zero axis predicate had a real residual.
  The focused test now verifies exact-zero predicates and both logarithm lips.
- Negative-mass dB0 zero-momentum and on-shell logarithms, mass-exchanged
  coincident roots, and the original normal-threshold prescription pass focused
  regressions. The normal-threshold prescription is an original API convention,
  not a finite one-sided derivative limit.
- The formerly failing finite complex four-mass box passes all four cyclic
  rotations at f64 and 128 bits using the new transparent `box_complex` map.
  This is a focused eight-comparison check, not full boxc certification.
- `tests/triangle_limits.rs` now passes all 20 points and all three coefficients
  at f64, 256 and 384 bits. This covers regular one-mass quadratic limits,
  negative masses, three-mass cyclic relabelings and an exact zero Gram
  determinant. A polynomial resultant detects the relabeling condition;
  evaluating the discriminant directly as `lambda/(m1*m2*m3)` preserves exact
  coincident roots. Neither fix uses a numerical proximity tolerance.
  The added independent `C0([-4,-1,-1],[0,0,1])=-ln(2)` control additionally
  agrees with a native 256-/384-bit logarithm to below 1e-60 absolute error.

The first follow-up construction attempt hit an internal
Symbolica representation error before evaluating points: recursively expanded
triangle relabeling grew into multi-gigabyte atoms. Compact non-inlined native
stages fixed construction; the triangle sweep above subsequently completed.

The new `tests/box_limits.rs` initially reported **31 failed comparisons** across
28 original-accepted cases, f64/128-bit evaluation and signed-zero variants.
After preserving continued products in IR06/12/13 and exact first-/second-order
contour directions, these cases all pass. Rectangular complex square roots avoid
losing tiny-width information through polar-angle rounding. Exact zero-endpoint
limits additionally fix the newly exposed three-mass box NaN. The test now
passes **63 cases** at f64/128 bits, including 31 three-mass acceptance fixtures,
exact endpoint-chart controls and signed-zero variants. The complete expanded
suite has additionally passed at 256 bits as recorded above. All three master
tests now pass, including 879 coefficients through both hook/map routes. The
final expanded fixed-f64 rerun has only the 14 bubble failures recorded above.

Additional concrete findings include lower-lip mass square roots, a missing
exact-zero-channel D0 pre-permutation, and accepted C0 quadratic-degenerate
points previously returned as zero. These repairs pass the focused family tests
above.
Cross-review identified and repaired removable coincident-contour-root dilog
limits and checked the ordinary contour chain rules against the source. The
remaining repeated-root counterexample is recorded in the
[independent threshold checks](#independent-finite-box-threshold-and-removable-limit-checks):
at `D0([0,0,0,0,4,-3],[1,1,1,1])` the native expression has an incorrect
imaginary part at both 128 and 256 bits. Its explicitly ignored diagnostic still
fails when invoked; it is not included in passing coverage. Infinitesimal
directions in the implementation are symbolic coefficients, not finite numeric
regulators or proximity-based switches.

### Exact formulas versus legacy numerical artifacts

The original can emit finite values without `ERROR` even when an internal
regulator/tolerance controls the answer. For example, finite C0 with
`p=[-1,-2,-3]`, `m=[1,2,3]` gives about -0.2090612533590607 at `mu²=1`, but
-0.2616240718822739 at `mu²=1e20` after a scale-dependent vacuum approximation.
A finite scalar triangle with zero Laurent poles is analytically independent
of that renormalization scale. Likewise, the exactly pinched one-mass C0
`p=[-2,3,2]`, `m=[0,0,1]` gets a large finite imaginary component from the
original's finite regulator rather than a regulator-independent finite limit.

Per the agreed exact-expression/no-near-degeneracy-gate design, these artifacts
are documented separately, not silently approximated or counted as proof that
the exact formula is wrong. Literal equality with *every* finite original
printed result is therefore a different contract from analytic scalar parity.

A further unresolved example is the unequal-mass vacuum triangle
`p=[0,0,0]`, `m=[-1,2,3]`. Original Fortran returns approximately
`-0.3618610961277855 - 0.2617993877991497 i`, but exchanging the first two
masses changes its imaginary part to `+0.7853981633974484 i`. The native
formula currently gives `+0.2617993877991494 i` in both orders. This is both a
native analytic-continuation gap and an inconsistency in the legacy oracle;
agreement with either arbitrary mass ordering cannot establish correctness.
The symmetric lower-lip vacuum formula
`-sum_i m_i*log(m_i-i0)/product_{j!=i}(m_i-m_j)` gives
`(2/3)*log(2)-(3/4)*log(3)-i*pi/12`, i.e. the negative imaginary part in
both orders. Choosing the intended analytic-versus-literal-legacy contract
and repairing this case remain pending; it is not counted as verified parity.

## Baseline audit — 2026-09-07

Audit date: 2026-09-07. Verdict: **full scalar parity is not established, and the
expanded audit has found reproducible disagreements.** Passing the earlier
acceptance suite is not evidence that every supported kinematic region works.
This report separates incorrect results, fixed-precision numerical failures,
legacy interface behavior, and points rejected by the original implementation.

## Original scope and conventions

The authoritative local interface description is
`../OneLOopBridge/oneloop/DESCRIPTION`, at original revision
`402bf0d70bfa2a9d8ae66b8251cc4395b7b5ab99`. Relevant locations are lines 26–42
(overloads), 82–100 (scale), 103–116 (on-shell threshold), 119–157 (scalar
functions and derivative), and 176–221 (tensor ranks). The original wrapper has
not been changed by this audit.

The author's [OneLOop paper, sections 1 and 2](https://arxiv.org/pdf/1007.4716)
describes scalar one- through four-point integrals, UV/IR dimensional
regularization, and collider-physics kinematics with nonpositive imaginary
squared masses. It explicitly does not claim complete treatment of exceptional
Landau-singular points. The newer checked-in interface/source is authoritative
for the particular version being compared here.

| Property | Original contract | Native expression contract / audit consequence |
| --- | --- | --- |
| Scalar families | A0, B0, C0, D0; separate dB0 derivative | All constructors exist; correctness is not complete. |
| Laurent output | Finite, simple pole, double pole | Same ordering and normalization. Higher epsilon orders are not promised by either scalar interface. |
| Mass inputs | Squared masses, real or complex, with Im(m²) ≤ 0 | Same intended mathematical domain. |
| External inputs | Real squared invariants, including spacelike, timelike, zero and on-shell cases | Same intended domain. Complex-typed original overloads discard nonzero imaginary external parts and report an error; they do not supply general complex-external analytic continuation. |
| Scale | Real unsquared `rmu`, default 1, optional per-call override | Explicit positive `mu_squared`; compare with `rmu=sqrt(mu_squared)`, not `rmu=mu_squared`. |
| Invalid mass imaginary parts | Documented error and sign reversal | Not a promised native coercion API. Actual A0 source lacks the B/C/D coercion block, so exact invalid-input compatibility needs per-family tests. |
| IR classification | Exact by default; optional absolute `olo_onshell` tolerance | Exact symbolic conditions. No approximate on-shell coercion API. |
| Numerical/backend controls | Precision variants, diagnostic units, errorcode, switches | Symbolica evaluator configuration; not legacy-runtime API parity. |

No positive-real-part mass restriction is stated in `DESCRIPTION` or enforced
by the public B/C/D input wrappers. The paper's collider-physics restriction
also means that arbitrary negative-real-mass behavior should not be advertised
as globally validated. Negative-real and purely imaginary masses must be tested
and their corner cases reported, not silently excluded by a new assumption.

The original forms `rmu*rmu` internally without checking its sign. Thus ±rmu
give the same squared scale; zero scale is not a meaningful finite-value
requirement for divergent integrals. Similarly, an original error-path return
value must not be treated as an independently correct analytic result.

The separate `ONI` numerical-integration routines provide an alternative
evaluation/checking method, not additional scalar families that require opaque
callbacks in this expression crate. Tadpole and bubble tensors through rank
four exist in both APIs, but this scalar audit does not independently certify
all tensor coefficients.

## Dispatch audit

| Family | Original paths inspected | Native status |
| --- | --- | --- |
| A0 | `avh_olo_a0.h90`, `tadp` in `avh_olo_bub.f90` | Zero and massive expressions present; negative-real-axis numerical disagreement found. |
| B0 | `avh_olo_b0.h90`, `bub0` | Massive/massless, zero momentum and scaleless expressions present. Small nonzero momentum has substantial f64 cancellation. |
| dB0 | `avh_olo_db0.h90`, `dbub0` | On-shell IR and zero-momentum cases present; regular coincident-root pseudothreshold limit missing. |
| C0 | `avh_olo_c0.h90:154`, `tria0–4`, `trif0–3`, `trif3HV` | All mass masks and IR dispatch present. Mixed-sign finite massless continuation fails. Generic zero-quadratic guards are not established analytic limits. |
| D0 | `avh_olo_d0.h90:195`, `box00`, `box06–16`, `boxf1–5`, `boxc` | All mass masks and named IR sectors present. Finite complex-mass and IR15/16 continuation failures found. |

In the original, finite three-/four-mass boxes select `boxc` when widths are
nonzero or the specified kinematic conditions require it
(`avh_olo_d0.h90:195–228`). Its T/T13/planar reductions and kinematic permutations
are not fully transcribed in the native crate. Alternate formulas need not be
ported line by line if equivalent expressions can be proved, but the current
failing samples rule out claiming equivalent coverage.

The original `trif3` also rotates masses according to their widths
(`avh_olo_tri.f90:531`); the corresponding generic native representation does
not. The three deliberately asymmetric-width cyclic samples tested so far
passed. This source difference is therefore an unproven risk, not a confirmed
additional defect.

## Reproduced disagreements

All quantities below are squared except the implicit Fortran `rmu=1`.
Every row uses `mu_squared=1`. Values are rounded here; generated fixture rows
retain the full input precision. Except where noted, the table displays the
finite coefficient. The original oracle returned finite values without an
`ERROR` diagnostic for these samples.

| Case and inputs | Native f64 result | Original Fortran result |
| --- | --- | --- |
| dB0: p=1, m=[1,4] (also [4,1]) | NaN | 0.07944154167983575 |
| C0: p=[−4.062768000049295, 2.4508858786367513, 4.004606533896626], m=[0,0,0] | −4.196483354270513 − 0.8237311565786827i | 0.476436102924887 − 0.8237311565786822i |
| D0 IR15: p=[−1,2,3,−4,8,−6], m=[0,2,0,3] | −0.04160609346860533 | −0.041606093468605325 + 0.11549261732149718i |
| D0 IR16: p=[2,−2,9,4,−5,12], m=[0,2,3−0.3i,4], **simple pole** | 0.014277663646144823 + 0.19688495323609248i | 0.02898319954998221 − 0.19526267086623783i |
| Finite four-mass D0, inputs immediately below | 0.29521879685976965 + 0.23433537622641903i | 0.046731572847697384 + 0.0780104270205361i |
| A0: m=−2 | −0.6137056388801094 + 6.283185307179586i | −0.6137056388801094 − 6.283185307179586i |

Finite four-mass D0 inputs:

```text
p = [-3.3906256351404696, -4.086683736180616,
     -2.880037208949436, -0.9160468158580961,
      8.562635640376207, -5.784338529774539]
m = [0.5542034093325066 - 0.000525812454592841i,
     0.3303230946406446 - 0.0006709837183118287i,
     1.4393815762877777 - 0.0001816189030170589i,
     1.5647400023647011 - 0.0006709386631579254i]
```

### Analytic versus numerical causes

- **dB0 pseudothreshold:** `db0_with_two_masses` in `src/two_point.rs` divides
  an identical-log difference by a zero discriminant. This is a regular
  pseudothreshold, not the divergent normal threshold. Its exact value is
  `3*log(2)-2`. Original `dbub0` handles coincident roots explicitly at
  `avh_olo_bub.f90:211`. Higher precision cannot supply a missing 0/0 limit.
- **IR15 branch loss:** `box_two_opposite_15` takes a principal logarithm of
  `(r13/r23)^2/r24`. Squaring loses the continuation phase of negative `r13`.
  In the displayed sample the missing imaginary part equals exactly
  `2*pi` times the simple-pole coefficient. Continued logarithms or explicit
  native sheet expressions are necessary; increasing precision is insufficient.
- **Massless C0:** the displayed real-part error equals exactly
  `-4*pi²/sqrt(lambda)`, with `lambda=71.37446769320898`. At this real-root
  point the native quadratic helper reverses the named roots relative to
  original `solabc` followed by `x -> -x`, but the triangle attaches opposite
  upper/lower lips as though their order were unchanged. This produces a
  continuation-sheet error, not an ordinary accuracy loss. Original `trif0`
  also derives the root lips from `-p1`, instead of using fixed signs.
- **IR16 and finite D0:** wrong finite and, for IR16, pole coefficients
  establish actual disagreements away from the displayed pseudothreshold.
  Raw products inside logarithms/dilogarithms and root-sheet assignments need
  targeted review. The present evidence does not prove a single common cause
  for every failing case.
- **A0 on the negative real axis:** the f64 evaluator disagrees with the
  intended lower-lip formula. Audit literal construction versus symbolic
  substitution, signed-zero treatment, and higher-precision evaluation before
  attributing this exclusively to the formula or the Symbolica evaluator.
- **Small nonzero bubble momentum:** at p=−1e−10, m=[1,2], dB0 evaluates in f64
  to approximately 2.07350913536e12 rather than 0.1137056388785528. The explicit
  1/p and 1/p² cancellations explain a fixed-precision risk. At p=+1e-10,
  m=[1,2], f64 gives about 2.063509138432e12; 128 bits gives
  0.11370563904907949 (still about 1.7e-10 off), while 256 bits gives
  0.11370563888166610, agreeing with Fortran to machine precision. The diagnostic
  `audit_near_zero_bubble_precision` records this. That confirms recoverable
  cancellation for this sample, not every failed bubble case. A `map_coeff`
  f64 evaluator is not itself a demonstration of automatic adaptive precision.

The separate analytic auditor independently cross-checked the pseudothreshold
limit and the IR15 phase diagnosis against both the native and original source.
In particular, original `avh_olo_box.f90:199–213` retains the squared Q phase that
the current native `box_two_opposite_15` loses. The cross-review ran no concurrent
Symbolica evaluator and made no production changes.

## Follow-up: three parallel mismatch investigations

The follow-up assigned A0/B0/dB0, C0 and D0 to separate auditors. The bubble
and triangle auditors cross-reviewed each other's proposed branch corrections;
the lead independently inspected the original and native IR15/IR16 expressions.
Production formulas were not changed in this follow-up. Independent ordinary-
Python mathematical evaluations below are diagnostic checks of proposed
corrections, not passing Rust regression tests or replacements for the oracle.

| Group | Failing coefficients | Follow-up conclusion |
| --- | ---: | --- |
| A0 negative real mass | 1 | Proven signed-zero-dependent logarithmic lip. |
| Small nonzero B0/dB0 momentum | 14 | Severe cancellation; independent Feynman-parameter checks support the Fortran values. |
| dB0 regular pseudothreshold | 2 | Missing coincident-root analytic limit. |
| Finite massless C0 | 4 | One wrong real-root/lip pairing, not four unrelated sectors. |
| D0 IR15 | 11 | Lost continued product phases and wrong real-axis logarithmic lip. |
| D0 IR16 | 11 | Same classes of continuation loss, affecting finite and pole terms. |
| Finite four-mass D0 | 1 | Reproduced; a complete causal correction remains unproven. |

### A0/B0/dB0

The shared `physical_log` in `src/lib.rs:108` is `conj(log(conj(z)))`.
Numerica's `Complex::conj` negates the imaginary component and its argument uses
`atan2` (`lib/numerica/src/domains/float/complex.rs:114`, `:145`, `:790`).
Consequently runtime `z=-2+0i` becomes `-2-0i` inside the log, and the outer
conjugation produces the upper, not lower, lip. Original `tadp` explicitly uses
`qonv(...,-1)` (`avh_olo_bub.f90:52`). Independent IEEE complex arithmetic
reproduced the failure and the opposite result for `-2-0i`.

A proposed exact native primitive is
`if(z-conj(z), log(z), log(abs(z))-i*pi*if(z-abs(z),1,0))`.
For nonzero real z this forces the intended axis value independently of signed
zero; off-axis it leaves the principal log unchanged. At zero it remains
singular. Because the helper is shared, a correction needs whole-suite testing,
including literal versus substituted arguments, both zero signs, shrinking
negative widths, and lazy conditional evaluation. The triangle auditor
independently checked these semantics.

Six B0 and eight dB0 failures occur at `p²=+/-1e-10`, covering equal, unequal,
complex and one-massless masses. Independent Simpson integration with 10,000
panels reproduced all eight pairs of Fortran values to about 1e-15, apart from
an approximately 4e-13 difference in one equal-mass derivative. This supports
the oracle and diagnoses the cancellation; it does not demonstrate that every
native failing case has been repaired by a particular precision setting.
The two exact pseudothreshold orders independently reproduce `3*log(2)-2`.

### C0

Every failing C0 is finite and fully massless, with `p1<0`, `p2,p3>0` and
positive Källén discriminant. Fixture lines are 159, 321, 325 and 329.
For the current root order, the proposed local change is
`qx1=SheetAtom::with_sign(root_1,-sign(r23))` and
`qx2=SheetAtom::with_sign(root_2,sign(r23))`, where `r23=-p1`.
It belongs in `triangle_finite_massless` (`src/triangle.rs:114`), not in the
shared quadratic solver. Original `trif0` (`avh_olo_tri.f90:293`) negates the
oppositely ordered roots and derives their real-axis signs from `r23`.

Independent continued-log/dilog evaluation with the corrected pairing matches
all four Fortran values within 1.2e-15. Twelve additional Fortran controls cover
positive-Källén Euclidean, all-timelike, mixed-sign and permuted inputs, also
within 1.2e-15. The bubble auditor independently verified the root-order/sign
argument for either sign of the quadratic coefficient. Nonreal roots keep
their actual imaginary signs; exact zero-discriminant limits are not solved
by this correction. These proposed changes still require Rust regression tests.

### D0

All 22 infrared disagreements are in fixture lines 655–685: IR15 and IR16 each
contribute seven finite and four simple-pole failures. Unlike original
`box15`/`box16`, native `src/box_integral.rs:577` and `:642` form raw products
and quotients before applying principal continued logs/dilogs. The original
retains `qmplx_type` sheets through these operations
(`avh_olo_box.f90:101–130`, `:197–214`). Negative-real `r24` also exposes the
same lower-lip problem identified in A0.

An independent mathematical re-evaluation using explicit original-style sheet
products and lower-lip logs matches all 16 IR15/IR16 fixture points (including
their passing coefficients) within 5e-16, without changing the scalar roots.
For IR16's displayed pole, using the lower lip gives
`0.02898319954998221-0.1952626708662378i`; using the upper lip reproduces the
current native disagreement. This isolates analytic continuation defects;
higher precision alone cannot fix them.

Original IR15 additionally orders by `abs(m2-p2)` and IR16 by `abs(m2)`
(`avh_olo_box.f90:162`, `:77`). Those swaps are absent in the native functions.
The diagnostic sheet expressions pass these samples both with and without
the swaps, so the omission is a robustness/source-alignment issue, not a proven
cause of these 22 failures. The single finite four-mass failure at line 587
needs a separate `boxc`/root-continuation investigation; the infrared correction
does not establish a remedy for it.

For that finite point, independent evaluation of the existing four-mass formula
with continued sheet arithmetic reproduces the wrong native value
`0.29521879685976965+0.2343353762264191i`. The roots have imaginary parts of
order 1e-4, not tiny real-axis roundoff. Four cyclic rotations give inconsistent
values, so neither an isolated signed-zero fix nor an arbitrarily chosen
rotation is an established repair. Original `avh_olo_d0.h90:199` routes these
complex widths to `boxc`. This is evidence of an inadequate continuation
representation in the port, rather than a Symbolica-only evaluation problem.

## Checks performed and remaining coverage

The previously passing acceptance corpus contains 293 points / 879 coefficients
(4 A0, 13 B0, 12 dB0, 125 C0, 139 D0), with 31 passing tests reported in
`AUDIT.md`. It covers every C0/D0 mass mask but not every kinematic branch within
each mask.

The final expanded exploratory run generated **342 accepted original-Fortran
points / 1,026 coefficients and found 44 coefficient disagreements**. Eleven
additional original-oracle points were rejected or nonfinite and were not
counted as successes. Family point counts are A0=13, B0=49, dB0=50, C0=100,
D0=130. Failing coefficient counts are A0=1, B0=6, dB0=10, C0=4, D0=23.
The scan includes positive-Källén massless triangles, permutations, negative/
zero-real-part complex masses, and independent mu_squared from 1e-12 to 1e12.

The retained inputs and Fortran outputs are `tests/data/scalar_audit.txt`, with
failure results in `tests/data/scalar_audit_failures.txt`. Fixture SHA-256:
`913e5b81f8bff8349c1bb3c75232c9c8c12a7f01b70b14e2abb136cacde4ca2d`.
These failing exploratory cases are deliberately separate from the previously
passing acceptance inventory; ignoring the optional diagnostic is not a parity
claim. The run used kinematic mass-dimension normalization independent of mu.

```sh
ONELOOP_AUDIT_FIXTURES=tests/data/scalar_audit.txt cargo test --test parity \
  audit_extra_fortran_fixtures -- --ignored --test-threads=1 --nocapture
```

The expanded checks exercise:

- A0/B0/dB0/C0/D0 separately, both bubble mass orders, normal-threshold
  neighborhoods, exact regular pseudothresholds, and small/zero momentum;
- all C0/D0 internal mass masks with real masses and multiple width scales,
  including finite complex boxes and unequal width orderings;
- spacelike and mixed-sign/timelike invariants, cyclic permutations, and IR15/16
  sign changes with both finite and pole coefficients checked;
- independent renormalization-scale changes with fixed kinematics, rather than
  only rescaling masses, invariants and mu together;
- explicit rejection of nonfinite native outputs and separate accounting of
  original diagnostic/error cases.

Reproduction uses `tests/support/oracle.f90` linked against the unchanged
original library, `tests/generate_audit_fixtures.py`, and the explicit ignored
`audit_extra_fortran_fixtures` test in `tests/parity.rs`. The generated expectation
values are original Fortran results, never results from the numerical Rust port.
The main audit records commands and retained fixture locations.

Before full-parity or stability sign-off, the remaining work is:

1. Preserve each failing sample as a regression and repair its analytic
   branch/limit or evaluator behavior; do not loosen comparison tolerances.
2. Cover every original scalar dispatch leaf, each zero/equality subcase,
   allowed mass ordering, and reflected/cyclic external ordering. Checking all
   mass masks is necessary but insufficient.
3. Generate physically consistent momentum configurations as well as algebraic
   stress inputs. Sweep cuts, positive/zero/negative Källén regions, shrinking
   widths, and hierarchies of kinematic and mass scales.
4. Recheck failures and cancellation neighborhoods at 128/256 bits against an
   independently higher-precision oracle or parameter integration where the
   double-precision original is ambiguous. Precision agreement alone does not
   verify the analytic branch.
5. Test scale identities at fixed kinematics. For
   `L=log(mu_squared_new/mu_squared_old)`, the coefficients must obey
   `F2'=F2`, `F1'=F1+L*F2`, `F0'=F0+L*F1+L²*F2/2`.
6. Keep legacy input/error behavior distinct from the native expression API.
   Do not add opaque numerical OneLOop helpers or claim unimplemented runtime
   controls are automatically supplied by Symbolica.

The previously noted C0 p=[−4,−1,−1], m=[0,0,1] is a mathematical finite-limit
gap, but the original oracle itself emits an error and returns zero there.
Likewise the tested all-zero-momentum boxes trigger original `boxc` errors;
the native equal-mass vacuum formula is supported by its analytic limit, not
by those invalid original outputs. Original `boxc` also rejects some arbitrary
external invariant tuples when no acceptable positive Källén permutation is
available (`avh_olo_boxc.f90:176–191`). These must remain separate from clean
Fortran-parity counterexamples.

## Independent finite-box threshold and removable-limit checks

The following checks use equal squared masses `[1,1,1,1]`, squared scale `1`,
and external invariants `[0,0,0,0,s,t]`. Both Laurent poles vanish. The
quadrature below is an independent **test oracle only**: it is not used by any
native expression, master-symbol callback, or production implementation.

The four-dimensional Feynman-parameter representation is
`D0 = integral_simplex [1-s*x1*x3-t*x2*x4]^-2`. Apply the Cheng-Wu choice
`x1+x2=1`, write `x1=x`, `x2=1-x`, `x3=r*v`, `x4=r*(1-v)`, and integrate
`v` and then `x` analytically. With `z=r/(1+r)` and `w=z*(1-z)`, this gives

```text
D0(s,t) = -integral_0^1 log[(1-s*w)*(1-t*w)] / (s+t-s*t*w) dz.
```

For the two real points below the logarithms have real limiting values. At
`s=4,t=-3`, the zero at `z=1/2` is an integrable logarithmic endpoint after
splitting the interval. At `s=2,t=-4`, numerator and denominator vanish
together, so their quotient must be evaluated by its removable limit.
Symmetry and `u=abs(1-2*z)` give cancellation-free one-dimensional formulas:

```text
D0(4,-3) = -integral_0^1 [2*log(u)+log((7-3*u^2)/4)]/(4-3*u^2) du
D0(2,-4) =  integral_0^1 log1p(u^2*(1-u^2)/2)/(2*u^2) du
          = pi/4 - atanh(1/sqrt(2))/sqrt(2).
```

The second integrand is `1/4` at `u=0`; its closed form follows immediately
by integration by parts. The independently obtained targets are

```text
D0(4,-3) = 0.437508120479433839303015369526783505697336155799893173453710894022789135269725375157557290 + 0i
D0(2,-4) = 0.162172923257217796221640765569307718398597037497209202345021371980784097024538069189339633 + 0i.
```

Reproduction with Python and `mpmath==1.3.0`:

```python
import mpmath as mp
for digits in (80, 120):
    mp.mp.dps = digits
    threshold = mp.quad(
        lambda u: -(2*mp.log(u) + mp.log((7-3*u*u)/4))/(4-3*u*u),
        [0, mp.mpf(1)/4, mp.mpf(1)/2, 1])
    removable = mp.quad(
        lambda u: mp.log1p(u*u*(1-u*u)/2)/(2*u*u) if u else mp.mpf(1)/4,
        [0, mp.mpf(1)/2, 1])
    exact = mp.pi/4 - mp.atanh(1/mp.sqrt(2))/mp.sqrt(2)
    print(digits, mp.nstr(threshold, digits), mp.nstr(removable, digits))
    print("closed-form error:", mp.nstr(removable-exact, 10))
```

The two precision runs agree through the displayed 80-digit result; the
closed-form errors were `-5.27e-82` and `-1.21e-121`, respectively.

The unchanged original build also provides `avh_olo_qp` with
`olo_kind=16`, `digits=113`, and machine epsilon
`1.92592994438723585305597794258492732e-34`. A separate temporary driver was
linked against its existing `libavh_olo.a` using gfortran 15.2.0. No original
source, module, archive, or wrapper build was changed. The zero-width results
were:

| Point | Original DP finite term | Original QP finite term |
| --- | --- | --- |
| `s=4,t=-3` | `0.437508112203920163 + 8.27551422065440851e-9 i` | `0.437508120479433831595842035686652605 + 7.70717333384013151e-18 i` |
| `s=2,t=-4` | `-0.134889937937259674 - 1.54027890413999557 i` | `0.127142510038012190864975536896963604 - 0.0838448301598513322241412737412247225 i` |

At the threshold, `Re+Im` from QP is
`0.437508120479433839303015369526784134`, agreeing with the independent
target to about 32 decimal digits. Explicit common negative widths `1e-24`,
`1e-30`, and `1e-32` converge to the same sum. The individual real and
imaginary displacements are order `sqrt(width)`; the original finite-regulator
imaginary part is not the exact threshold value. In contrast, the zero-width
`s=2,t=-4` original result is unstable even in QP. A common width `1e-16`
gives approximately `0.1621729232572177966441 + 3.77745e-17 i`, consistent
with the independent closed form, but shrinking widths reintroduce severe
cancellation. Merely obtaining finite original output is therefore insufficient
to make either zero-width result a strict regression target.

To reproduce QP with the existing original build, copy
`tests/support/oracle.f90` to a temporary directory, import `avh_olo_qp`
instead of `avh_olo`, and declare its momenta, masses and results as
`complex(olo_kind)`. Keep the scale as `real(kind(1d0))`, matching the original
`olo_scale` API. Print `real(result(i),kind=olo_kind)` and `aimag(result(i))`
using `2es45.35`, and link the driver with `-I /path/to/original/OneLOop`
and `/path/to/original/OneLOop/libavh_olo.a`. This uses the native QP module;
it does not mix promoted driver types with a DP archive.

The native exact `s=4,t=-3` result is a **confirmed remaining gap**: its real
part matches the independent target, but it has the spurious imaginary part
`+3.5830460942481976` at both 128 and 256 bits. With common masses `1-i*delta`
and strictly positive width magnitude `delta`, the native values instead
converge correctly. Independent quadrature at 80 and 120 decimal digits gives,
for example,

```text
delta=1e-8:
  0.437452584442429213196616652862811396
  + 0.0000555328011958603090209985543129131 i
delta=1e-30:
  0.437508120479433283942648099731002629
  + 5.55360367269795457351641215380708e-16 i.
```

For this complex-mass check the same parameter integral applies with
`D=m*(s+t)-s*t*w` and integrand `-log1p(-D*w/m^2)/D`, where `m=1-i*delta`
denotes the squared mass. Splitting the transformed `u` interval at
`sqrt(delta)` resolves the threshold boundary layer. The exact repeated inner
polynomial is `4*(z-1/2)^2`; perturbation by `-i*delta` splits its roots at
order `sqrt(delta)`, not order `delta`. A Taylor-only continuation cannot
silently assign this term the same order as ordinary root directions.
Documenting this discrepancy is not a production repair or a passing exact
threshold regression.

The native exact `s=2,t=-4` value, on the other hand, agrees with its independent
closed form at 128 and 256 bits; a test expecting either unstable zero-width
original output would be testing the wrong reference. At width `1e-30`,
128-bit evaluation exhibits about `6.6e-9` cancellation drift while 256-bit
evaluation agrees with the independent value. This is distinct from the
precision-independent exact-threshold branch gap above.

A bounded native repair has been derived and independently cross-reviewed,
but is **not implemented** in this snapshot. In `s3`, for real `a>0,b,c` with
`b^2-4*a*c=0`, repeated root `z=-b/(2*a)` strictly inside `(0,1)`, and both
outer poles outside the closed real interval `[0,1]`, replace the coalescing-root
logarithms by `log(a)+2*log(abs(x-z))` and integrate the two intervals with
native dilogarithm divided differences. The kernel is bounded on this domain,
so the integrable logarithmic limit can be taken directly. This avoids a full
fractional-power root expansion for the confirmed counterexample. Simultaneous
inner/outer poles inside the interval are excluded from that argument and
require a limit of the combined contour expression, not the individual `s3`.
