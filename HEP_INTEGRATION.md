# Integration with Symbolica community HEP

Inspected on 2026-09-21: community main
`b7a2b327dd92181e46941a7ce14aa0a94d24dc2c` and alphal00p/gammaloop's Feynkit branch
`58aaf5380023d118a91aa4083c01395586246f14`.
The local community checkout has additional work, so the data-ownership
assessment uses the upstream Git trees. The local host integration described
below additionally uses that checkout's Feynkit paths.

## Data ownership

| Interface | Recommendation | Reason |
| --- | --- | --- |
| `Atom`, `Symbol`, `Complex<T>`, `Float`, `DoubleFloat` | Share Symbolica/Numerica types directly | OneLOop already uses these; Feynkit tensor operations consume/return `Atom` and its kinematics uses Numerica traits. |
| `FourMomentum<T>`, `ThreeMomentum<T>`, `MomentumSignature` | Reuse Feynkit when adding momentum or diagram frontends | `FourMomentum::dot` and `mass_squared` use the mostly-minus metric and generic scalar arithmetic. Convert momenta to OneLOop's existing squared-invariant arguments at the boundary. |
| `LaurentSeries` and three numeric coefficients | Keep the small OneLOop representation for now | Feynkit has no equivalent exact Laurent-series container. OneLOop's order is explicitly `(0, -1, -2)`. |
| Vakint `NumericalEvaluationResult` | An explicit future conversion, not a core dependency | It stores arbitrary `(epsilon power, Complex<Float>)` pairs, not exact `Atom` coefficients; importing Vakint solely for this type would couple unrelated engines and precision policies. |
| `TensorReducer`, `TensorReduction` | Reuse for a future vacuum tensor frontend | Their scope is vacuum tensor reduction. General external-momentum B/C/D tensor reduction is not supplied by this API. |
| `FeynmanDiagram`, model `ModelExpression` | Keep outside the integral engine | Diagrams need routing and propagator selection; `ModelExpression` is a named string transport object, not an exact expression replacement. |
| Evaluators and `DecimalComplex` | Keep OneLOop-specific adapters | Feynkit does not provide the matching evaluator-cache or lossless Python decimal contract. Its current Python momenta store `f64`; Rust generic momenta can preserve higher precision. |

If multiple engines need a common Laurent container, extract a small generic
`LaurentSeries<C>` with explicit epsilon powers into a shared HEP crate. Do not
make OneLOop depend on Vakint to reuse one result type. Preserve a checked
conversion to OneLOop's fixed three coefficients, and keep function-map ownership
separate from coefficient storage. This is a proposed future extraction, not an
API change in this update.

`NativeFloat` also needs more than Feynkit's `KinematicScalar`: it supplies
precision construction and the dilogarithm callback. The common Numerica traits
are already shared, so replacing this local evaluator trait would not remove
those requirements.

Do not infer a shared epsilon normalization from common types: measure factors,
scale conventions, mass widths and the coefficient power ordering must be agreed
at each integral-engine boundary. No graph-to-master mapping is implemented here.

## One kernel and one public namespace

The adapter now provides `register_hep_module(&hep)` for a host to attach
`hep.oneloop` to the existing HEP module. Classes identify themselves as
`symbolica.community.hep.oneloop`. The minimal development host exercises this
path and preserves `symbolica.community.oneloop` as an alias with identical
class objects. It does not include Feynkit itself.

In community `src/hep.rs`, immediately after Feynkit registration, add:

```rust
oneloop_native::register_hep_module(module)?;
```

Initialize `oneloop_native::CommunityModule` through
`SymbolicaCommunityModule::initialize(py)` in `HepModule::initialize`, alongside
Feynkit initialization. OneLOop's hook performs no numerical work: formulas and
each requested numerical backend remain lazy. The Symbolica state inventory
registers only the lightweight master and JIT callbacks.
Add the `oneloop-python` dependency with `default-features = false` and
`features = ["community"]`; the Rust library name is `oneloop_native`.
The public Python HEP package already re-exports the native module's attributes.
Add corresponding `hep.oneloop` stubs when merging into the full distribution.

OneLOop's manifests now declare a versioned Symbolica dependency and select its
Git revision through root-level `[patch.crates-io]` entries. This is deliberate:
community main also uses versioned dependencies plus a root patch. Direct `rev`
and `branch` Git dependencies are different Cargo source identities even when
they resolve to the same commit. Mixing them can produce incompatible `Atom`
and `PythonExpression` types and independent Symbolica state.

A consuming host must select the same tested Symbolica and Numerica revision as
OneLOop's cache assets. Dependency-manifest patches do not propagate. Check
`cargo tree -d` / `cargo metadata` for a single Symbolica and Numerica package;
matching Python version strings alone does not establish shared state. A future
community-wide dependency update must also regenerate OneLOop's version-bound
cache assets. Native JIT is currently required by OneLOop; a WASM host must gate
this module until that dependency is made optional.

## Composition with Feynkit

After registration in the full community host, existing APIs suffice:

```python
from symbolica import S
from symbolica.community import hep

p = hep.FourMomentum(3.0, 1.0, 0.0, 0.0)
finite, pole, double_pole = hep.oneloop.B0(p.mass_squared, 4, 4)

s = S("s", is_real=True)
coefficients = hep.oneloop.master_coefficients(S("oneloopmaster::B0")(s, 4, 4, 1))
amplitude = hep.oneloop.compile_native([coefficients[0]], [s])
```

The first route intentionally has binary64 momentum inputs. For arbitrary
precision use Symbolica expressions/Decimal inputs or Rust
`feynkit_kinematics::FourMomentum<symbolica::domains::float::Float>` and pass its
squared invariants without an intermediate `f64` conversion. Tensor-reducer
outputs are already Symbolica expressions and can be multiplied by individual
Laurent coefficients before compiling; keep the shared function map via
`OneLoopExpressions`/`compile_native`.

## Local full community host

The checkout at `/common/dev/symbolica-community/main` now links the native
OneLOop adapter into `HepModule` alongside Feynkit. Its root manifest pins
Symbolica, Numerica, and Graphica to `a19c760`; its lockfile uses SymJIT 2.26.0.
OneLOop's registration and initialization are gated out on WASM. The full host
provides `symbolica.community.hep.oneloop`, with Python type hints; the legacy
`symbolica.community.oneloop` alias belongs only to the minimal development host.
The full host uses the system allocator after its original mimalloc build
crashed during the Python suite's thread transitions; see the
[allocator report](patches/community-host-allocator-2026-09-22.md).

Use the existing community environment (including NumPy):

```sh
cd /common/dev/symbolica-community/main
.venv-feynkit/bin/python examples/oneloop_smoke.py
```

The example passes Feynkit's `FourMomentum.mass_squared` into `B0`, compares
native and SymJIT output, composes OneLOop coefficients with Feynkit's tensor
reducer using shared Symbolica expressions, and evaluates an arbitrary-precision
integral. Restart an existing notebook kernel after rebuilding the extension.

To run OneLOop's Python tests against this full host:

```sh
cd /common/dev/oneloopmaster
ONELOOP_PYTHON_MODULE=symbolica.community.hep.oneloop \
  /common/dev/symbolica-community/main/.venv-feynkit/bin/python \
  -m unittest discover -s python/tests -v
```

The full host passes all 28 applicable OneLOop Python tests (two
standalone-only checks skipped), including all benchmark fixtures and the
worker transitions that crashed its original mimalloc configuration.
The 13 host HEP/community checks also pass after refreshing four missing
integral-family declarations from Feynkit's existing stubs. The runnable smoke
example verifies momentum inputs, backend agreement, shared expressions, and
50-digit output in the installed `.venv-feynkit` environment.
