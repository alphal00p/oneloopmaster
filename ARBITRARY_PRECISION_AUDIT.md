# Focused arbitrary-precision audit

This independent audit, completed on 2026-09-08, traced precision through the
Python/Rust interfaces, scalar-master hooks, generated arithmetic, and the
patched Symbolica/Numerica operations actually used by the masters. It is a
focused numerical audit, **not** a full integration, publication, or
all-kinematics accuracy certificate. See [PRECISION.md](PRECISION.md) for the API
and its working-precision contract.

The subsequent [clean integration audit](NATIVE_RELEASE_AUDIT.md) passes 114
rebuilt release Rust tests and 27 standalone Python tests. It also records the
extreme-point reference-conditioning issue exposed by the corrected precision
metadata, without relaxing the comparison tolerance.

## Findings and repairs

Here, `value@N` means a stored Float with `N` binary precision bits. A component
can retain few relative bits after cancellation while its absolute uncertainty
is negligible compared with another component.

| Path | Reproduced problem | Repair and resulting behavior |
| --- | --- | --- |
| Complex square root | A tiny real component at 2 bits could reduce an otherwise accurate root to 2 bits through the half-angle divisor. Even after fixing that divisor, `sqrt(2@512 + i·1e-154@2)` returned real part `1.5@2`. | A Float-specific scaled Cartesian root avoids the low-precision phase/cosine path. Tests retain the accurate dominant component and the appropriate branch lip. |
| Phase/logarithm | `arg(-1@512 + i·1e-154@2)` returned `3.0@2`, rather than a phase close to π. | Dominant-axis reduction computes only the small correction at its available precision, then combines it with an accurately constructed π or π/2. The probe now returns a phase near `3.141592653589793…` at 512 bits; a genuinely uncertain tiny angle is not padded. |
| Integer powers | `(1e-154@2 + 1.25i@512)^2` returned real part `-1.0@2`, instead of `-1.5625`. | Nonzero powers start from their first factor, avoiding an unnecessary multiplication by a weak-precision identity and its signed-zero effects. The probe now returns `-1.5625@512`. Zero powers construct the exact identity from the strongest component. |
| Half-integer powers | General polar evaluation remained reachable for negative half powers; `(2@512 + i·1e-154@2)^(-3/2)` returned real part `0.375@2`, rather than approximately `0.353553390593273762…`. | Explicit ±1/2 and ±3/2 paths use square roots, products, and reciprocals. They avoid reintroducing the phase/cosine precision loss. |
| Reciprocal ordering | The first `-3/2` specialization formed an overflowing intermediate for binary64 `1e210`, returning zero instead of the representable result near `1e-315`. | Native evaluation takes reciprocals before multiplying; the generic expression path cubes the inverse root. Real/complex subnormal controls pass. |
| Precision metadata | Adding a 512-bit exact zero to stored `1.5@2` advertised `1.5@512`. Raw max-precision `hypot` could similarly inflate a weak dominant component's precision, including off-axis cases. | Zero addition/subtraction retains the nonzero operand's precision. Finite native Float norms use tracked scaled arithmetic; nonzero axes preserve the original component precision. These repairs do not create an error-bound certificate. |

The native operations are in [primitives.rs](src/native/primitives.rs). The
Numerica changes to `domains/float/{complex,multiprecision}.rs` and Symbolica
changes are distributed in [symbolica-dev.patch](patches/symbolica-dev.patch);
[patch documentation](patches/README.md) records the dependency baseline.
These findings concern the port/evaluation dependencies, not a newly claimed
Fortran defect.

## Reviewed precision boundaries

- [Python conversion](python/src/precision.rs): Decimal components and large
  integers are parsed directly at explicit working bits, without binary64
  conversion. Decimal output strings avoid ambient Decimal-context rounding.
  Ordinary Python floats retain their supplied binary64 values; requesting
  more precision cannot recover digits absent from those inputs.
- [PrecisionEvaluator](src/precision.rs) and
  [NativeEvaluator](src/native/mod.rs): decimal-to-bit sizing, input rounding,
  batch/domain validation, and output handling. The fixed guard is not an
  adaptive accuracy guarantee; outputs are not padded after evaluation.
- [Master hooks](src/masters.rs): Float cache selection includes every real and
  imaginary component, including the final `mu_squared` argument. DoubleFloat
  conversion to the Li2 kernel retains both compensated components.
- [Generated constants](src/native/generated_constants.rs) and
  [generator](examples/generate_native.rs): exact rational strings and symbolic
  π metadata are initialized at the requested precision. Expression mapping
  resolves symbolic constants at explicit working bits rather than accepting
  placeholder values.
- Reachable Li2 routines in the dependency patch: exact Float predicates,
  precision-scaled checked convergence limits, guarded continuation, and MPFR
  real Li2/integer-zeta calculations. No finite-value f64 narrowing or fixed
  900-bit stopping cap was found in the reviewed arbitrary-precision path.
  Unconverged partial sums are not silently returned as successful results.

## Generated-operation inventory

Static inspection and constant propagation through the generated helper calls
resolved every power instruction in the audited snapshot:

- **917** native integer-power calls, all exponent **-1**. Positive integer
  scalar powers are explicit multiplication chains.
- **30** noninteger-power calls: **25** exponent **-1/2**, **5** exponent
  **-3/2**. There were no unresolved/dynamic exponents in this inventory.

Thus the scalar-master calls do not require the remaining general polar-power
fallback. This inventory should be repeated if the generated formulas change.

## Focused test provenance

[complex_precision_provenance.rs](tests/complex_precision_provenance.rs)
passed **5/5** against the fresh normal-release artifact
`libnumerica-86f02b6d5e2fd079.rlib`, built at approximately 16:36 UTC on
2026-09-08. The direct `rustc --test` run used that artifact as
`--extern symbolica=…`: these tests use only the Numerica-reexported `domains`
API and do not initialize Symbolica State. Runtime was 0.01 seconds. Coverage
includes 512/3456-bit component imbalances, large phases versus uncertain tiny
angles, integer-power identities, signed axes, and ±1/2/±3/2 powers.

Three additional focused tests in
[native_primitives.rs](tests/native_primitives.rs) passed against the then-current
native source, directly compiled with the existing debug Symbolica artifact
`libsymbolica-5e04a1d83362a6d5.rlib`: unbalanced-precision half powers,
binary64 reciprocal-order/subnormal controls, and integer-power signed lips
for f64/Float. These were not a fresh full dependency/integration build. The
subsequent tracked Float-norm repair was cross-reviewed separately and has its
own value/precision regression in that test file. The owning release audit
records the complete rebuilt Rust/Python suite results separately.

## Remaining bounds on the conclusions

No unaddressed finite-value narrowing/provenance defect was identified in the
reviewed ArbPrec scalar paths. This does not prove 1,000 correct digits for every
input: cancellation in composite expressions, very small output components,
singular conditioning, finite exponent ranges, and already-rounded inputs remain
limitations of fixed working precision. Float precision metadata is approximate,
not a rigorous error enclosure. SymJIT remains an explicitly binary64 backend.

The generic Symbolica binary64 complex polar-square-root extreme-range behavior
is distinct from the native f64 scaled implementation and the repaired Float
implementation. DoubleFloat also has independent signed-zero normalization in
compensated arithmetic: the new private non-inverse integer-power lip tests do
not claim universality for that type. Those integer-power instructions are not
present in the audited native master DAG. Neither caveat is a claim that the
ArbPrec masters silently fall back to binary64.
