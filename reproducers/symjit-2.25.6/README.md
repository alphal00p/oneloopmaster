# SymJIT 2.25.6: standalone bug reproducers

Six small Rust programs reproduce bugs encountered while porting OneLOop.
They depend only on the unmodified crates.io `symjit = "=2.25.6"`, using its
public `Translator`, `Composer`, and `Applet` APIs. No Symbolica, OneLOop,
license, custom callbacks, or dependency patches are required.

Each source file is self-contained and can accompany a separate bug report.
`Cargo.lock` pins the complete dependency resolution. The published SymJIT
crate records upstream commit `3fc04010f69db954463f9666fffb652b244ccc52`.

## Run everything

From this directory, with Rust/Cargo and Python 3 on PATH:

```sh
python3 run.py
```

The runner builds the six binaries, runs twelve cases in separate processes,
and writes [observed.txt](observed.txt), including the environment and complete
output. It continues after assertion failures or a child process's segmentation
fault and disables core dumps. A reproducer's nonzero exit is expected on the
affected version; the four control runs should exit zero. The runner itself is
a transcript collector, not a pass/fail test harness.

Individual commands are below. Run them from this directory without a
`symjit.toml` configuration file or `SYMJIT_TOML` override. Each program uses
the crate's default configuration plus the explicit options visible in its
source. All JITs use O2 and `fastmath=false`; SIMD and threads are disabled
except for the explicit SIMD batch case. Cargo's debug profile is sufficient:
these options control the generated JIT code independently of the Rust profile.

## Observed environment and results

Verified on 2026-09-19 with Rust 1.98.1, Linux x86-64, and an AMD EPYC 9754
(AVX, AVX2, FMA, AVX-512 available). The SIMD reproducer explicitly disables
AVX-512. These results describe this release and machine; other versions and
architectures have not been tested.

| Report | Observed failure | Passing control |
| --- | --- | --- |
| Complex square root | Wrong branch for negative imaginary input in generic mode; wrong signed-zero lip in both modes | Positive imaginary input and upper cut lip |
| Complex condition | `if(i, 7, 11)` returns 11 | Real nonzero and exact zero conditions |
| Strict complex multiplication | `z*conj(z)` acquires an imaginary residue | `fast_complex=true` |
| Direct fractional power | Compilation panics: `label _func_power_ not found` | Indirect translation |
| Direct negative powers | `2^-2 = 4`, `2^-3 = 8` | Indirect translation |
| Zero-parameter SIMD batch | `evaluate_matrix` raises SIGSEGV | Scalar evaluation and SIMD-disabled batches |

Seven reproducer invocations fail with Rust exit code 101; the eighth raises
SIGSEGV. All four separate control invocations exit zero. Most exit-101 cases
are the MRE's own assertions detecting incorrect results. The fractional-power
case instead panics inside SymJIT's assembler during compilation.

## 1. Complex square root chooses the wrong branch

Source: [complex_sqrt.rs](src/bin/complex_sqrt.rs)

```sh
cargo run --locked --bin complex_sqrt
cargo run --locked --bin complex_sqrt -- generic
```

The only generated operation is `symbolica_sqrt` on one complex parameter.

With `fast_complex=false`, `sqrt(0 - 4i)` returns approximately
`-1.414213562373095 + 1.4142135623730951i`. The principal root is
`+sqrt(2) - sqrt(2)i`.

With either value of `fast_complex`, `sqrt(-4 - 0i)` returns `-0 + 2i`.
Preserving the lower lip of the branch cut requires `+0 - 2i`. Here `-0i`
is IEEE negative zero, supplied directly as a runtime argument.

The same program checks `sqrt(0 + 4i)` and `sqrt(-4 + 0i)` as passing controls.
Expected values are explicit mathematical constants, not another JIT result.

## 2. Complex conditions ignore a nonzero imaginary component

Source: [complex_if.rs](src/bin/complex_if.rs)

```sh
cargo run --locked --bin complex_if
cargo run --locked --bin complex_if -- join-only
```

The first variant emits `IfElse`, branch assignments, and `Join` for
`if(z, 7, 11)`. The second emits only `Join` with constant alternatives,
isolating conditional selection from branch execution.

Both return 11 for `z=+i` and `z=-i`, where a complex nonzero predicate should
select 7. They correctly return 7 for `z=1` and 11 for `z=0`.
All arguments are declared complex; none are marked real.

## 3. `fastmath=false` still produces an imaginary residue in `z*conj(z)`

Source: [complex_rounding.rs](src/bin/complex_rounding.rs)

```sh
cargo run --locked --bin complex_rounding
cargo run --locked --bin complex_rounding -- packed
```

For `z=0.1 + 0.2i`, with `fast_complex=false`:

```text
actual   = 0.05000000000000001 + 1.6653345369377347e-18 i
expected = 0.05000000000000001 + 0 i
```

The imaginary part is an exact cancellation of two identical binary64
products under separate multiplication/addition rounding. The MRE also prints
the ordinary Rust complex product as a reference. This is a strict-rounding
report: the residue is small, but changes exact-zero predicates and branch-cut
selection in consumers. `fastmath` is explicitly disabled.

Changing only `fast_complex` to true yields the expected zero imaginary part.

## 4. Direct translation cannot compile a fractional power

Source: [direct_power.rs](src/bin/direct_power.rs)

```sh
cargo run --locked --bin direct_power
cargo run --locked --bin direct_power -- indirect
```

The program appends one `Powf` instruction for `z^(1/4)`, with the exponent
stored as a constant. Direct translation panics inside `assembler.rs`:

```text
label _func_power_ not found
```

Expected: compilation succeeds and evaluation at `z=16+0i` returns `2+0i`.
The same instruction stream compiles and evaluates correctly with indirect
translation.

## 5. Direct translation computes positive powers for exponents -2 and -3

Source: [direct_negative_powers.rs](src/bin/direct_negative_powers.rs)

```sh
cargo run --locked --bin direct_negative_powers
cargo run --locked --bin direct_negative_powers -- indirect
```

Two `Pow` instructions request `z^-2` and `z^-3`. At `z=2+0i`, direct
translation returns `[4+0i, 8+0i]`; the expected results are
`[0.25+0i, 0.125+0i]`. Changing only to indirect translation produces the
expected results. This reproduces the reciprocal-power specialization problem
independently of the missing fractional-power target.

## 6. SIMD batch evaluation segfaults for a function with no parameters

Source: [zero_input_batch.rs](src/bin/zero_input_batch.rs)

```sh
ulimit -c 0
cargo run --locked --bin zero_input_batch
cargo run --locked --bin zero_input_batch -- scalar
```

The function has **zero parameters**, two constant complex outputs `[2, pi]`,
and no callbacks or branches. It first passes a scalar evaluation check, then
calls the safe API `app.evaluate_matrix(&[], &mut output, 5)` with ten allocated
complex output elements: five rows of two outputs.

With SIMD enabled, the program confirms a SIMD kernel exists and then raises
SIGSEGV. Disabling SIMD produces `[2, pi]` in all five rows. This is a function
with an empty parameter list, not a batch of zero-valued parameters or a batch
with zero rows. The reproducer itself contains no unsafe code.

The crash requires a machine on which SymJIT actually compiles a SIMD kernel;
the program prints whether that occurred. If it prints `compiled=false`, that
run has not exercised the failing SIMD path.
