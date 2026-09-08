# Inspecting the complete expressions

`get_expression` returns the actual three Laurent coefficient bodies, in
`[epsilon^0, epsilon^-1, epsilon^-2]` order. It recursively substitutes the
transparent FunctionMap definitions, including every sheet, triangle and box
helper. A bounded one-mass triangle sector also uses its analytically
integrated compact identity, described below. It does not return compact
evaluatable master calls or evaluator IR.
The last argument is always the squared renormalization scale `mu_squared`.

In Rust:

```rust
use oneloop::{get_expression, ScalarIntegral};
use symbolica::{atom::AtomCore, printer::PrintOptions};

let m = symbolica::symbol!("m_squared"; Positive).to_atom();
let mu = symbolica::symbol!("mu_squared"; Positive).to_atom();
let coefficients = get_expression(ScalarIntegral::A0, &[m, mu])?;
println!("{}", coefficients.coefficients()[0].printer(PrintOptions::full()));
```

`ScalarIntegral::get_expression` is an equivalent convenience method. Run
`cargo run --example get_expression` to print A0 and massless C0 bodies. As with
the other examples, Symbolica work stays on one explicitly enlarged-stack thread;
normal startup still prepares the eager numerical backends.

For the standalone Python extension:

```python
import oneloop_native as olo

finite, pole, double_pole = olo.get_expression(
    "A0", ["m_squared", "mu_squared"],
    positive=["m_squared", "mu_squared"],
)
print(finite)  # complete, parseable Symbolica text
pole_only = olo.get_expression("A0", ["m_squared", "mu_squared"], coefficient=-1)
```

Standalone arguments and returned coefficients are **strings**, not Expression
objects from a separately loaded Symbolica extension. Unqualified input names
belong to `oneloop_input`. Output includes namespaces and native attributes for
faithful parsing. Alternatively, input strings can carry native attribute syntax,
for example `"kinematics::{real}::s"` or `"kinematics::{positive}::mu_squared"`.
`real=[...]` and `positive=[...]` create attributed variables before parsing;
existing variables must already satisfy the requested attribute. The API never
retags an existing symbol. A successful declaration persists in this kernel.

In a shared-kernel community build, the same function instead accepts and returns
actual Symbolica Expressions:

```python
from symbolica import S
from symbolica.community import oneloop

m, mu = S("m_squared", "mu_squared", is_positive=True)
finite, pole, double_pole = oneloop.get_expression("A0", [m, mu])
```

Use the host's `S(..., is_real=True)` / `S(..., is_positive=True)` when creating
community variables; there are no separate `real` / `positive` keyword lists in
that build. The existing `master_coefficients` remains the compact Symbol API.

## Assumptions, limits and numerical use

The complete one-scale triangle is compact for a proven-real invariant `s` and a
possibly complex mass squared `m`, in the configuration
`C0(0, 0, s, 0, m, 0, mu_squared)` or either cyclic relabeling. For `s != 0` and
`Im(m) < 0`, its finite coefficient is

```text
(pi^2/6 - Li2(1 + s/m)) / s
```

For an exactly real nonzero mass, native conditionals select the Li2 upper lip
when `s > 0` and its lower lip when `s < 0`; the latter is essential for negative
real masses. Off the mass axis, ordinary complex Li2 is used. This continuation
follows from `d[s*C0]/ds = log(-s/m)/(m+s)` and `s*C0 -> 0` as `s -> 0`.
The full returned expression also retains native `if` branches for `s = 0`
(the infrared vacuum Laurent coefficients) and `m = 0` (the massless Laurent
coefficients); it is not just the finite formula above. At `s = -m`, the Li2
argument is exactly zero and no artificial endpoint division is introduced.
No real or positive mass attribute is required; an untyped mass is allowed to be
genuinely complex. Physical inputs must still satisfy `Im(m) <= 0`. The general
definitions remain available subject to
the explicit expansion budget; fully symbolic generic C0/D0 expressions are not
promised to expand compactly within the default budget.

The representative complete three-coefficient output has 220 expression nodes
with a positive-real mass and 316 with an untyped, possibly complex mass. Both
fit the default expansion budget without opaque helpers.

The [real-sector regression](tests/inspection.rs) compares all three coefficients
on 864 rows, keeping the compact expression at 384 working bits and the original
FunctionMap reference at 512 bits, with its unchanged `1e-60` comparison target.
The reference needs extra precision near zero: at `s/m = -1e-30`, a 384-bit
evaluation retains only 188 output bits and differs from a 768-bit reference by
about `3.3e-55`. The 512-bit reference differs by `1.9e-93`, while the compact
384-bit result differs by `4e-84`. A separate 384/512/768-bit convergence
regression preserves this evidence. This is fixed-working-precision cancellation
in the original expression, not an analytic-continuation mismatch or an automatic
precision guarantee.

The [complex continuation grid](tests/inspection_complex.rs) passes all 4,410
rows at 512 bits, including positive/negative real masses, purely imaginary
masses, widths down to `1e-30`, three mass positions, three scales, both invariant
signs, exact/near zero, and the retained zero-mass branches. Another 12-point
complex-mass check compares the complete expression at 3,456 bits with the
original expression at 3,840 bits using a `1e-1000` target, for all three
coefficients. These finite test sets are not proofs of global analytic or
numerical parity.

For example, the standalone interface can construct that full sector directly:

```python
finite, pole, double_pole = olo.get_expression(
    "C0", ["0", "0", "s", "0", "m", "0", "mu"],
    real=["s"], positive=["mu"],
)
```

Known real subexpressions simplify conjugation; native normalization then removes
conditions that become exact numbers. Undecidable conditions remain native `if`
expressions. String tags alone do not imply realness or positivity. In particular,
Symbolica's positivity predicate includes zero in some structural cases, so it is
**not** used to discard zero-test branches. Real masses are optional: complex
squared masses with nonpositive imaginary parts remain supported by the formulas.
Physical inputs still require real external invariants and positive real scale;
constructing an expression does not prove its eventual numeric inputs obey these
domain requirements.

Full inlining duplicates shared subexpressions and generic boxes can become very
large. Rust `get_expression_with_options(..., ExpressionOptions { max_nodes,
max_depth })` and Python's corresponding keywords control the expansion budget.
Defaults are 1,000,000 cumulative visited/materialized nodes and depth 512.
Memoized copies count toward work; this is not a final-output-only node limit.
If a limit or definition cycle is encountered, the call fails explicitly. It
never silently returns an incomplete expression with opaque OneLOop helpers.
Supplying actual kinematic relations or attributed inputs can prune branches
before their bodies are expanded. Selecting one Python coefficient currently
selects from the complete three-coefficient expansion, so the same budget applies.

For numerical work, the compact `OneLoopExpressions` FunctionMap route preserves
shared evaluation and argument binding. The fully inlined formula is exact, but
does not inherit that route's floating-point conditioning or performance. Neither
inspection nor simplification implies a new numerical-stability/parity guarantee.
