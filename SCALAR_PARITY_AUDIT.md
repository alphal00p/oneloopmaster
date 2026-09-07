# Scalar-integral parity audit

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
