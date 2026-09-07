# Native OneLOop expressions

An independent Rust library whose only direct dependency is Symbolica. It does
not link the Fortran wrapper or the standalone numerical Rust implementation.
All functions construct native symbolic expressions; continuation definitions
are transparent `FunctionMap` entries, never numerical OneLOop callbacks.

Development status: the expanded scalar audit found 44 coefficient mismatches
against accepted original-Fortran results. See `SCALAR_PARITY_AUDIT.md`; this
snapshot does not have full scalar parity and is not publication-ready.

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
