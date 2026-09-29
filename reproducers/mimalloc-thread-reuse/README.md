# mimalloc: retired initial thread and reused thread identity

This reproducer uses only C, pthreads, dlopen, and mimalloc. It has no Rust,
Python, Symbolica, SymJIT, licensing, or numerical-expression dependency.

1. A worker loads a DSO containing mimalloc and allocates/frees 64 bytes.
2. The worker exits, running mimalloc's pthread cleanup.
3. A subsequent worker reuses the original pthread identity and allocates.

Bundled mimalloc 3.3.2 crashes at step 3. The `main` control first loads and
uses the DSO on the process main thread, which stays alive; this passes.

## Reproduce

Use the C source in `libmimalloc-sys-0.1.49/c_src/mimalloc/v3` as `SOURCE`.
The build deliberately retains `-ftls-model=local-dynamic`, matching the
community host's `local_dynamic_tls` feature.

```sh
bash build.sh "$SOURCE" /tmp/mimalloc-before
python3 check.py /tmp/mimalloc-before/driver /tmp/mimalloc-before/libprobe.so --expect-crash
python3 check.py /tmp/mimalloc-before/driver /tmp/mimalloc-before/libprobe.so --main
```

Copy the source into a temporary directory, apply the
[backport](../../patches/mimalloc-3.3.2-thread-reuse.patch) there, and repeat:

```sh
cp -a "$SOURCE" /tmp/mimalloc-patched-source
git -C /tmp/mimalloc-patched-source apply /absolute/path/to/mimalloc-3.3.2-thread-reuse.patch
bash build.sh /tmp/mimalloc-patched-source /tmp/mimalloc-after
python3 check.py /tmp/mimalloc-after/driver /tmp/mimalloc-after/libprobe.so
```

The backport contains the relevant `src/init.c` and `src/alloc.c` changes from
upstream commit [`b92c116b67d0db73f09653695061e5838cf82b8e`](https://github.com/microsoft/mimalloc/commit/b92c116b67d0db73f09653695061e5838cf82b8e),
with added-line trailing whitespace removed. Other commit hunks concern tests,
comments, or unrelated cleanup. The source is MIT-licensed by Microsoft/Daan
Leijen. Apply this patch only to a separate copy, not Cargo's registry cache.

## Cause and evidence

`_mi_theap_free` explicitly clears `theap_main.tld` when the initial worker
exits. However, `tld_main.thread_id` keeps that worker's identity. A later
pthread can reuse the identity; `_mi_is_main_thread()` then incorrectly
selects the retired static main heap. Its initialization is skipped because
the process-global heap is already initialized, and `_mi_subproc` dereferences
the null thread-data pointer.

The [GDB watchpoint transcript](watchpoint.txt) catches the actual store to
null in `_mi_theap_free`, followed by the failing read in `_mi_subproc`.
This is an ordinary software write, not an inferred RAM bit flip.

Capture it with `gdb -batch -x watch.gdb --args /tmp/mimalloc-before/driver
/tmp/mimalloc-before/libprobe.so` (one shell command).

The upstream fix clears the retired thread identity and distinguishes it from
the initial/uninitialized sentinel. Future workers consequently get fresh
thread metadata. It is part of C mimalloc 3.4.0 and later. See
[upstream issue #1287](https://github.com/microsoft/mimalloc/issues/1287).

Twenty fresh processes were tested for each case:

| Allocator | First use | Result |
| --- | --- | --- |
| Bundled 3.3.2 | Worker | 20 SIGSEGVs |
| Bundled 3.3.2 | Main | 20 passes |
| Upstream fix commit | Worker | 20 passes |
| Exact bundled 3.3.2 + backport | Worker | 20 passes |
| Bundled v2 alternative | Worker | 20 passes |

The same before/after result held with both CPU and memory bound to NUMA node
0 (CPU socket 0) and node 4 (CPU socket 1): five crashes before and five passes
after on each node. This is not a claim that the hardware is healthy; the
debugger identifies the ordinary software store that causes this particular
failure.

The [Rust cdylib probe](rust-probe/) independently confirms the published-crate
behavior: 20/20 crashes with default v3, and 20/20 passes with `--features v2`.
Build it with Cargo and pass the resulting `libmimalloc_thread_reuse_probe.so`
to the same C driver/check script. The probe uses `MiMalloc` as Rust's global
allocator and an ordinary `Box` allocation.

As checked on 2026-09-22, crates.io's latest `mimalloc` is 0.1.52 and
`libmimalloc-sys` is 0.1.49. Both the release and the Rust wrapper's current
Git revision still use C v3.3.2. The published `v2` feature is an alternative:

```toml
mimalloc = { version = "0.1.52", features = ["local_dynamic_tls", "v2"] }
```

That selects an older allocator implementation; it does not backport the v3
fix. A fixed v3 currently requires patching `libmimalloc-sys` or updating its
vendored C source until the Rust wrapper publishes an update. A host must also
enable Symbolica's `faster_alloc` if it wants Symbolica to use mimalloc globally.

## Full community-host confirmation

A separate test build of the complete community host used the exact
`libmimalloc-sys` 0.1.49 sources with this backport, Symbolica `faster_alloc`,
and mimalloc's `local_dynamic_tls` (without `v2`). With the supplied license
accepted, the complete OneLOop Python suite passed: 28 tests passed and two
standalone-only checks were skipped. In particular, the worker transition that
crashed the unpatched build passed. This build was kept under
`/tmp/oneloop-mimalloc-patched-host`; the installed community host continues to
use its tested system-allocator configuration.
