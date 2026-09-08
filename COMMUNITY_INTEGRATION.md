# Community integration: inspected implementation

This is a plan, not an integration or publishing change. The original Fortran
wrapper and the standalone numerical crate remain separate and unchanged.

## Current adapter progress (2026-09-08)

The thin binding crate now exists under [python/](python/README.md), with an
optional `community` feature implementing `SymbolicaCommunityModule`. It returns
genuine host expressions for the compact master Symbols and offers `compile_native`
to retain the transparent function map for combined expressions. Its standalone
PyO3 extension exchanges numbers only; separately loaded Symbolica binaries must
not exchange their internal expression objects.

Both adapter build modes compile. The standalone debug extension passes all
seven API tests on CPython 3.13.12, including original fresh-worker mode
(107.040 s), and the unchanged two-worker regression passes. Registration in a
real community host, stub generation, and shared-host runtime verification are
still future steps; these standalone results do not establish shared expression
state in a host. Numerical performance measurements remain pending.
No community-root files have been changed in this pass. The current dev snapshot
is `fb845d34`; a consuming host must also carry the documented numerica/graphica
and SymJIT overrides. See [patches/README.md](patches/README.md).

Core and adapter manifests disable Symbolica default features and explicitly
retain `tracing_max_level_info`, `integer-gmp`, `float-mpfr`,
`native_code_generation`, and `bincode`; community mode additionally enables
`python_export`. This excludes the optional `faster_alloc` global mimalloc
allocator. A sequential-thread crash was traced to that allocator during PyO3
argument extraction, before evaluator construction; the system-allocator rebuild
passes the unchanged checks. Cargo features are unified across the host graph,
so any dependency requesting Symbolica defaults or `faster_alloc` can undo this
selection. The consuming root must inspect its resolved features and verify the
actual host allocator configuration, not assume dependency-level
`default-features = false` disables a feature requested elsewhere.

## What Spenso actually does

The community checkout depends on the Python binding crate `spynso3`, not just
the numerical/symbolic `spenso` core. Its manifest currently overrides that
dependency with gammaloop revision `fff1d66ee7ca039b9e165fe8d29da91f4c27113c`.
`spynso3::SpensoModule` implements Symbolica's `SymbolicaCommunityModule` trait:

- `get_name()` returns `spenso`.
- `register_module()` adds Python classes/functions without registering Symbolica
  symbols. The trait explicitly reserves symbolic initialization for `initialize()`.
- `initialize()` initializes the package when the Python package is imported.

The community root calls `register_module!(m, spynso3::SpensoModule)` inside
`src/lib.rs::core`. This creates `symbolica.community.spenso_native`, registers it
in `sys.modules`, and exports `initialize_module`. The Python facade in
`python/symbolica/community/spenso/__init__.py` imports that native module and
calls `initialize_module()`.

The root `python_stubgen` feature forwards `spynso3/python_stubgen`; the shared
`src/bin/stub_gen.rs` gathers the bindings' annotated type information. Maturin
builds one `symbolica.core` extension containing the common Symbolica engine.
The README's `register_extension` example is outdated: follow the actual
`register_module!` macro and `example_extension/src/lib.rs` instead.

## Proposed OneLOop steps

1. Keep `oneloop` as a separate Rust library with **only Symbolica** as a direct
   dependency. Do not add Python dependencies or the Fortran wrapper here.
2. Add a thin sibling binding crate (for example `oneloop-python`) depending on
   `oneloop`, compatible Symbolica/PyO3, and optional stub generation support.
   Follow `example_extension` and Spenso's core/bindings split.
3. Expose native expressions, coefficient order, and a reusable expression
   context/function map. Returning a bare C0/D0 expression without its definitions
   is insufficient. Preserve sharing when combining integrals into an amplitude;
   do not replace definitions with opaque numerical callbacks.
4. Implement `CommunityModule` with name `oneloop`. Register Python types first;
   create any Symbolica symbols/definitions only during `initialize()` or later.
   Respect Symbolica's thread/license policy when deciding where construction runs.
5. Add the binding crate to the community root manifest, forward its stub feature,
   and register it in `core`. Add `python/symbolica/community/oneloop/__init__.py`
   using the existing native-import/initialization pattern.
6. Align the entire dependency graph to **one compatible Symbolica revision**.
   The community root currently patches Symbolica to `main`, whereas OneLOop uses
   `dev` plus local evaluator fixes. Dependency patches in a nested library manifest
   do not propagate to a consuming root. Move the agreed overrides to the build
   root and check for duplicate Symbolica/numerica/graphica versions before binding.
7. Generate Python stubs with the existing `python_stubgen` binary, build a local
   wheel with Maturin, and test imports, shared expression state, mixed C0/D0
   expressions, complex masses, Laurent ordering, and the Fortran regressions.
   These are local validation steps; no CI is requested yet.
8. Before release, replace local development paths with an upstream revision
   containing the required evaluator fixes (or an explicitly agreed maintained
   patch), settle version/package metadata, and obtain the licensing permission
   the project owner is handling. Nothing here assumes permission has already
   been granted or changes the original Fortran licensing.

## Sources inspected

- Community `Cargo.toml`, `src/lib.rs`, `src/bin/stub_gen.rs`, `pyproject.toml`.
- `example_extension/Cargo.toml` and `example_extension/src/lib.rs`.
- `python/symbolica/community/spenso/__init__.py` and its generated stub.
- Gammaloop's pinned `crates/spynso3/Cargo.toml` and `src/lib.rs`.
- Symbolica dev `src/api/python.rs::SymbolicaCommunityModule`.

This describes the checked-out sources, not a claim about a later upstream HEAD.
