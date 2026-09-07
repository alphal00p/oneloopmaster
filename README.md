# Native OneLOop expressions

An independent Rust library whose only direct dependency is Symbolica. It does
not link the Fortran wrapper or the standalone numerical Rust implementation.
All functions construct native symbolic expressions; continuation definitions
are transparent `FunctionMap` entries, never numerical OneLOop callbacks.

Development status: the expanded scalar audit found 44 coefficient mismatches
against accepted original-Fortran results. See `SCALAR_PARITY_AUDIT.md`; this
snapshot does not have full scalar parity and is not publication-ready.

This work is based on **Andreas van Hameren's OneLOop package**. The original
algorithms, analytic continuation conventions and Fortran reference results are
his work; see [Credits and references](#credits-and-references).

Inputs are squared invariants, squared masses with `Im(m²) <= 0`, and a positive
real squared renormalization scale. External invariants are real. Coefficients
are ordered `[finite, 1/epsilon, 1/epsilon²]`, with OneLOop normalization. This is
a Laurent expansion through the finite term, not an unexpanded function of epsilon.

## Development dependency

The manifest currently uses the sibling `../symbolica-dev` checkout for three
local fixes: non-inlined evaluator external-function indices, exact native
absolute values on the real axis, and native polylogarithm convergence near the
unit circle. Symbolica is based on dev revision
`0b57776b`. See [patches/README.md](patches/README.md) for a fresh-checkout
bootstrap. Keep Symbolica, numerica and graphica compatible. This temporary path
setup is **not a publishable dependency configuration**. A consuming project's
root manifest must supply these overrides too; library-level patches do not
propagate to consumers.

## Quick start

Clone this repository as `oneloopmaster`, then follow the sibling Symbolica
[bootstrap instructions](patches/README.md). A recent Rust toolchain supporting
edition 2024 and a working Symbolica license/setup are required (the audited
toolchain is Rust 1.91.1). No Fortran compiler is needed for ordinary use or tests.

```sh
git clone https://github.com/alphal00p/oneloopmaster.git
# Set up the patched sibling symbolica-dev as described above, then:
cd oneloopmaster
cargo run --example basic
```

[examples/basic.rs](examples/basic.rs) is a complete runnable example. Its core is:

```rust,ignore
use symbolica::prelude::*;

let p = parse!("p");
let mass_squared = Atom::num(1);
let mu_squared = Atom::num(1);
let bubble = oneloop::b0(&p, &mass_squared, &mass_squared, &mu_squared);
let views = bubble.coefficients().each_ref().map(Atom::as_view);
let mut evaluator = Atom::evaluator_multiple(&views, &[p])
    .build()
    .expect("compile B0")
    .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
let mut output = [Complex::new(0.0, 0.0); 3];
evaluator.evaluate(&[Complex::new(-1.0, 0.0)], &mut output);
println!("[finite, simple pole, double pole] = {output:?}");
```

At `p²=-1`, `m0²=m1²=mu²=1`, this gives approximately
`[-0.1520447048200202, 1, 0]`, with zero imaginary parts. The example also prints
the exact A0 expressions. To evaluate complex masses, make the mass arguments
symbolic parameters and supply `Complex::new(real_part, nonpositive_imag_part)`.
Use complex-valued evaluation even for real inputs: timelike integrals can have
an imaginary part. See the precision and stack requirements below.

The scalar entry points are `a0(m², mu²)`, `b0(p², m0², m1², mu²)`,
`db0(p², m0², m1², mu²)`, `c0([p1²,p2²,p3²], [m0²,m1²,m2²], mu²)`,
and `d0([p1²,p2²,p3²,p4²,s12,s23], [m0²,m1²,m2²,m3²], mu²)`.
Here `s12=(p1+p2)²`, `s23=(p2+p3)²`; ordering follows the original OneLOop API.
`db0` differentiates B0 with respect to `p²`, not the renormalization scale.
Arguments are references to Symbolica `Atom`s. The routines construct the three
Laurent coefficients; they do not numerically integrate Feynman parameters.

## Combining integrals

`a0`, `an`, `b0`, `db0` and bubble tensor coefficients return plain native Laurent
series. `c0` and `d0` return `MappedLaurentSeries`, including the definitions needed
to evaluate them. Their `_mapped` names are aliases. For amplitudes, construct a
single `OneLoopExpressions` context and share its definitions:

```rust,ignore
use oneloop::OneLoopExpressions;
use symbolica::prelude::*;

let context = OneLoopExpressions::new();
let p = [parse!("p1"), parse!("p2"), parse!("p3")];
let m = [parse!("m1"), parse!("m2"), parse!("m3")];
let mu2 = parse!("mu2");
let c = context.c0([&p[0], &p[1], &p[2]], [&m[0], &m[1], &m[2]], &mu2);
let combined = c.coefficients()[0].clone() + oneloop::a0(&m[0], &mu2).coefficients()[0].clone();
let parameters = [p[0].clone(), p[1].clone(), p[2].clone(),
                  m[0].clone(), m[1].clone(), m[2].clone(), mu2];
let exact_evaluator = context.evaluator(&[combined], &parameters)?;
// Map exact coefficients with map_coeff for f64, or map_coeff_with_prec for Float.
```

Keep the function map when exporting expressions or building an evaluator through
Symbolica directly. Compact function calls alone are not self-contained formulas.
The context caches the native definitions once per process; cyclic mass sectors
reuse canonical formulas instead of constructing copies.

The current constructors are **not** public A0/B0/C0/D0 Symbols with attached
`EvaluationInfo` hooks. A0/B0/dB0 expand their formulas; C0/D0 use internal sector
Symbols. The requested master-Symbol API remains an open audit item, documented
with the relevant Symbolica APIs in [MASTER_SYMBOL_AUDIT.md](MASTER_SYMBOL_AUDIT.md).

## Stack and precision

Generic C0/D0 compilation still traverses a large conditional expression graph.
The regression tests run the entire Symbolica workload on a thread with a 128 MiB
stack. For a standalone program, create that thread **before constructing any
Symbolica values**, and run construction/evaluation inside it. A restricted
Symbolica license permits only one active Symbolica thread per user. Do not
construct atoms in one thread and then transfer work to a second active thread.
Applications already using Symbolica must arrange an adequate calling-thread
stack rather than silently spawning another Symbolica worker.

Use exact rational inputs when testing precision convergence. Mapping coefficients
and supplying values as `Complex<f64>` is fixed-precision evaluation; it cannot
recover digits already rounded out of the input. `Complex<Float>` evaluators can
be built at a chosen precision using `map_coeff_with_prec`. No custom adaptive
precision or near-degeneracy gate is implemented here.

## Verification

Run local checks from this directory:

```sh
cargo test -- --test-threads=1
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

The acceptance suite currently passes 31 tests. The broader retained audit has
342 accepted Fortran points and 44 failing coefficients; run it explicitly with:

```sh
ONELOOP_AUDIT_FIXTURES=tests/data/scalar_audit.txt cargo test --test parity \
  audit_extra_fortran_fixtures -- --ignored --test-threads=1 --nocapture
```

That command currently **fails by design to expose the known mismatches**; the
failures have not been fixed or converted into passing assertions. Passing only
the default test suite is not a full-parity certificate.

The permanent Fortran fixture generator is standard-library-only Python. Compile
the small test driver against an existing, unchanged original OneLOop build:

```sh
gfortran -I /path/to/OneLOop tests/support/oracle.f90 \
  /path/to/OneLOop/libavh_olo.a -o /tmp/oneloop-oracle
python3 tests/generate_fixtures.py /tmp/oneloop-oracle
```

It prints the fixture table to standard output and rejects rejected/incomplete or
nonfinite Fortran results. The original example program does not print the full
A0/B0/dB0 outputs needed by this table, so use this test driver. The original
library is not a production dependency. Set its compiler runtime library path if required
by your local toolchain. The checked-in fixture table makes normal tests independent
of Fortran, Python, and the standalone numerical port.

See `AUDIT.md` and `STATUS.md` for verified coverage and known incorrect or
unsupported cases. This is a development implementation, not a claim of global
Fortran equivalence. Cargo publishing remains disabled. See
`COMMUNITY_INTEGRATION.md` for the inspected Spenso-style integration plan.

## Credits and references

The mathematical and algorithmic basis of this transcription is
[OneLOop, by Andreas van Hameren](https://helac-phegas.web.cern.ch/OneLOop.html).
The original Fortran implementation is the independent reference used by the
scalar parity tests. This repository is a separate, still-incomplete Rust /
Symbolica transcription, not an official replacement or a claim of endorsement
by the original author. The original
[OneLOopBridge wrapper](https://github.com/SecretGmG/OneLOopBridge) remains unchanged.

Please credit the original work when using or building on these expressions:

- A. van Hameren, **OneLOop: for the evaluation of one-loop scalar functions**,
  *Computer Physics Communications* **182** (2011) 2427–2438,
  [arXiv:1007.4716](https://arxiv.org/abs/1007.4716),
  [doi:10.1016/j.cpc.2011.06.011](https://doi.org/10.1016/j.cpc.2011.06.011).
- A. van Hameren, C. G. Papadopoulos and R. Pittau,
  **Automated one-loop calculations: a proof of concept**,
  *JHEP* **09** (2009) 106,
  [arXiv:0903.4665](https://arxiv.org/abs/0903.4665),
  [doi:10.1088/1126-6708/2009/09/106](https://doi.org/10.1088/1126-6708/2009/09/106).

[Symbolica](https://symbolica.io/) supplies the symbolic expressions, function
maps and numerical evaluation infrastructure. It has its own licensing terms.
Permission for the intended distribution/relicensing of this transcription is
being handled separately; it has not been recorded as granted here. The manifest's
license identifier is provisional, not a grant of rights to the original code.
Cargo publication is disabled pending permission and technical readiness.
