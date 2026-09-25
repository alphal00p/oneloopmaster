# Inspecting the complete expressions

`get_expression(master)` returns the actual three Laurent coefficient bodies,
in `[epsilon^0, epsilon^-1, epsilon^-2]` order. It recursively substitutes the
transparent FunctionMap definitions, including every sheet, triangle and box
helper. A bounded one-mass triangle sector also uses its analytically integrated
compact identity, described below. The input is a Symbolica master call with the
squared scale last and **no leading Laurent tag**; the result contains all three
coefficients, not compact evaluatable master calls or evaluator IR.

In Rust:

```rust
use oneloop::{get_expression, A0};
use symbolica::{atom::AtomCore, printer::PrintOptions};

let m = symbolica::symbol!("m_squared"; Positive).to_atom();
let mu = symbolica::symbol!("mu_squared"; Positive).to_atom();
let coefficients = get_expression(A0().call((&m, &mu)))?;
println!("{}", coefficients.coefficients()[0].printer(PrintOptions::full()));
```

`ScalarIntegral::get_expression` and `get_expression_for_family` are equivalent
enum-based conveniences. Run `cargo run --example get_expression` for A0 and
massless C0. Examples use an explicitly enlarged-stack thread; normal startup
still prepares the eager numerical backends.

For Python, build the [single-kernel host](python/README.md) or register the
adapter in a compatible Symbolica community host:

```python
from symbolica import S, E, N, Replacement
from symbolica.community import oneloop as olo

m, mu = S("m_squared", "mu_squared", is_positive=True)
master = S("oneloopmaster::A0")(m, mu)
finite, pole, double_pole = olo.get_expression(master)
print(finite)  # A genuine Symbolica Expression.
pole_only = olo.get_expression(master, coefficient=-1)
```

All symbolic arguments and results are actual Symbolica expressions. Arbitrary
composite arguments are supported, not just variables. Use Symbolica's `S` or
`E` to construct them; OneLOop never parses string arguments, retags variables,
or maintains a separate assumption system. For example:

```python
other_master = E("oneloop::B0(psq,mz_masses::{real,positive}::m2+12/3,m3,mu_r)")
other_coefficients = olo.get_expression(other_master)
```

The `oneloop::` namespace is an inspection shorthand for `oneloopmaster::`.
`master_coefficients(master)` returns the canonical compact coefficient calls
with native evaluation hooks. Only a context that requires exactly a Symbol,
such as the evaluator family selector, converts a variable Expression to Symbol;
integral arguments are never restricted to symbols.

The legacy standalone numeric extension has no symbolic API. A separately
loaded Symbolica extension has separate internal state even at the same revision;
its expressions must not be passed across that boundary.

## Selecting branches while retaining parameters

```python
psq = S("branch_psq", is_real=True)
m2 = S("branch_m2", is_positive=True)
master = S("oneloopmaster::B0")(psq, m2, m2, 1)
all_branches = olo.get_expression(master)
rules = [Replacement(psq, N("3.23")), Replacement(m2, N(1))]
selected = olo.select_branch(all_branches, rules)
finite_only = olo.select_branch(all_branches[0], rules)
assert selected[0] == finite_only
# selected has no ifs in this example; its finite part still contains psq/m2.
```

Rust uses Symbolica's native `Replacement` objects too:

```rust
use oneloop::{B0, get_expression, select_branch};
use symbolica::prelude::*;

let s = symbol!("branch_psq"; Real).to_atom();
let m = symbol!("branch_m2"; Positive).to_atom();
let all = get_expression(B0().call((&s, &m, &m, 1)))?;
let rules = [Replacement::new(s.clone(), Atom::num((323, 100))),
             Replacement::new(m.clone(), 1)];
let selected = all.coefficients().each_ref().map(|c| select_branch(c, &rules));
```

The helper applies native replacement rules only to a temporary copy of each
`if` condition, and uses Symbolica numerical evaluation to resolve zero/nonzero.
It never applies probe values to the returned branches. Every nested `if`,
including those inside conditions and arbitrary surrounding expressions, is
processed using an explicit work stack, without a helper-imposed depth cap.
Inactive branches are skipped. Unknown or undefined predicates remain `if`s;
their original parametric conditions are retained, with nested selections only.
Rules use native matching, namespaces, ordering and restrictions. Misspelled or
wrongly namespaced variables do not match; partial rules may leave conditions.
Nested conditions are selected first; each resulting condition is then sampled
with one simultaneous native replacement pass.

Integer, fraction, floating-point and complex numeric atoms are supported, with
no conversion to binary64 or zero tolerance for decisions. Finite floats are
represented by their exact binary rational value during substitution; temporary
native guards defer arithmetic until numerical evaluation, avoiding premature
rounding and expensive symbolic integer factorization. These guards never enter
the returned expression. See the [actual B0 before/after output](examples/expressions/b0_equal_mass.txt)
and run `cargo run --release --example select_branch` to reproduce it.
Probe points must respect the existing assumptions and physical domain. Numerical branch
selection is not a formal proof of transcendental equalities or a guarantee at
ill-conditioned boundaries. The selected formula is valid only in the selected
analytic region (or on a boundary if the probe lies there), not globally. There
is no implicit algebraic rewrite to a particular textbook form.

## Selecting before expansion (large C0/D0 expressions)

`select_branch(get_expression(master), rules)` must finish expanding all branches
before selection starts. For large masters this can exhaust `max_nodes` even
when the selected expression is small. Pass the rules into the expansion instead:

```python
mass2, mass2B, mu2 = S("mass2", "mass2B", "mu2")
master = S("oneloopmaster::C0")(0, -mass2, mass2, mass2, mass2, mass2B, mu2)
rules = [Replacement(mass2, N(2)), Replacement(mass2B, N(1)),
         Replacement(mu2, N(1))]
selected = olo.get_expression(master, branch_rules=rules, max_nodes=100_000_000)
assert olo.select_branch(selected, rules) == selected
value = selected[0].evaluate({mass2: 2, mass2B: 1, mu2: 1},
                             decimal_digit_precision=60)
```

The Rust entry points are `get_expression_on_branch(master, &rules)` and
`get_expression_on_branch_with_options(master, &rules, options)`. They return
the same `LaurentSeries` type as the all-branches entry points.

This uses the existing transparent FunctionMap definitions, not a separate
numerical backend or a hard-coded formula for the probe. Each expanded condition
is sampled using the same protected native replacements as `select_branch`;
only the chosen arm is expanded if decidable. Nested selections happen inside
conditions as well as function arguments and returned bodies. Undecidable
conditions keep their parametric form and both arms are expanded. All caches
are local to one expansion and its fixed rules; no result is reused across
different probes. Expansion work and recursion limits still apply.

This massive C0 still exceeds the default one-million-node **cumulative work**
budget. With branch probes, the original 100-million-node budget is sufficient;
without probes it is not. The returned generic expression is not automatically
rewritten to the compact four-dilogarithm identity.

At `(mass2, mass2B, mu2) = (2, 1, 1)`, the expression obtained through these
repo APIs gives

```text
finite      = -0.313413985805891556410851729326115326931925200907803748216227
simple pole = 0
double pole = 0
```

Run `python examples/c0_selected_branch.py --full` in the shared Python host to
print the actual if-free parametric expression and verify its value. The
independently derived four-dilogarithm identity is only a numerical check in
that script, never the source of the returned expression. The Python regression
also compares the same selected formula at `(2,1,7)`, `(4,2,3)` and `(3,1,1)`
with that identity and the native backend, and checks positive-tagged inputs.
These checks establish this case and nearby selected-region behavior, not
global correctness of every possible branch probe.

Focused regression commands (from the repository root, with the shared Python
host installed) are:

```sh
RUST_MIN_STACK=134217728 cargo test --release --lib inspection::tests -- --test-threads=1
cargo test --release --test branch_expansion --test branch_selection --test inspection --test inspection_complex --test inspection_master -- --test-threads=1
PYTHONPATH=python/tests ONELOOP_PYTHON_MODULE=symbolica.community.oneloop python -m unittest test_inspection test_precision -v
```

The explicit Rust harness stack makes isolated initializer tests independent of
which test first initializes Symbolica. This is not a request to change the
library's runtime threading. The wider, untargeted library suite was not completed
in this fix's validation. Rust 1.98 Clippy also reports pre-existing
`chunks_exact_to_as_chunks` suggestions in the unchanged backend, native and
precision batching code; the scoped lint check allows that lint only.

`branch_rules=None` (the default) retains the original all-branches behavior;
an explicit list enables branch-aware expansion. `coefficient=` can be combined
with `branch_rules`. Parameters are never replaced in the returned arithmetic.
Rules must be native `Replacement` objects, not strings or pairs. As with
`select_branch`, a chosen formula is only valid in the sampled region, and
partial or non-finite probes may leave undecidable conditions.

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

For example, construct that full sector directly from a master Expression:

```python
s = S("triangle_s", is_real=True)
m = S("triangle_m")  # May be complex.
mu = S("triangle_mu", is_positive=True)
finite, pole, double_pole = olo.get_expression(
    S("oneloopmaster::C0")(0, 0, s, 0, m, 0, mu)
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
