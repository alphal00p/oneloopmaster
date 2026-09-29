# Community host allocator crash

**Diagnosis established:** this is mimalloc's initial-thread retirement/thread-ID
reuse bug, already fixed upstream in C mimalloc 3.4.0. An isolated C reproducer
and a hardware watchpoint show the exact software store to null and subsequent
invalid reuse; Symbolica, SymJIT, Python, and licensing are not required.
See the [reproducer, controls, and backport](../reproducers/mimalloc-thread-reuse/README.md)
and [upstream issue #1287](https://github.com/microsoft/mimalloc/issues/1287).
The earlier hardware uncertainty below records the investigation context;
the ECC problem is real but does not explain this demonstrated software defect.
A separate full-host build with mimalloc v3 restored, `local_dynamic_tls`
enabled, and the upstream fix backported also passes all 28 applicable Python
tests (two skipped). The installed host was not replaced by this test build.

The full local community host reproduced a native crash with Symbolica
`a19c760dd567c239f30d87e4e924ca2f8b8457ab`, SymJIT 2.26.0, mimalloc 0.1.52 /
libmimalloc-sys 0.1.49 (`local_dynamic_tls`), and CPython 3.11.16 on Linux x86-64.
The host was built with Maturin's dev profile and Symbolica at opt-level 2.

From the OneLOop checkout:

```sh
ONELOOP_PYTHON_MODULE=symbolica.community.hep.oneloop \
  /common/dev/symbolica-community/main/.venv-feynkit/bin/python \
  -m unittest discover -s python/tests -v
```

The API and backend tests passed. The process then exited with SIGSEGV (139)
on `test_all_nested_ifs_in_functions_powers_and_deep_branches`, after moving to
a new test worker. The workers run sequentially, with 128 MiB stacks. The
systemd core dump from PID 1016612 captured this stack:

```text
_mi_subproc
mi_heap_main
mi_thread_init
_mi_malloc_generic
...
mi_malloc_aligned
alloc::raw_vec::RawVecInner::try_allocate_in
...
directories::ProjectDirs::from
symbolica::license::try_acquire_lock
symbolica::license::LicenseManager::acquire_restricted_thread_permit
symbolica::state::State::get_global_state
symbolica::atom::DefaultNamespace::attach_namespace
symbolica::api::python::expression::PythonExpression::symbol
symbolica::api::python::symbol_shorthand
```

This is an allocator-path failure, not an incorrect integral result or a
license-rejection exception. It does not yet establish the underlying cause in
Symbolica versus mimalloc, nor exclude earlier corruption. The lightweight
[thread control](../reproducers/community_thread_handoff.py) passed with
mimalloc, so that control alone is not a reproducer.

The local community host selects Symbolica's native numeric/code-generation
features explicitly, without `faster_alloc`, and removes its direct mimalloc
dependency. Its example extension also avoids `symbolica/default`, which would
otherwise re-enable the allocator through Cargo feature unification. No
license checks or upstream Symbolica source are changed.

## TLS and licensed controls

Cargo's compiled fingerprints confirm `local_dynamic_tls` on both mimalloc
crates. The sys crate selects `-ftls-model=local-dynamic` for that feature.
This failure is not the loader's static-TLS-space exhaustion error.

GDB resolves the failing instruction to mimalloc v3 `src/init.c:435`:
`return theap->tld->subproc`. Here `theap == &theap_main` and
`theap_main.tld == NULL`, so the faulting read is at address `0x18`.
The follow-up watchpoint identified `_mi_theap_free` as the writer during
thread cleanup; the subsequent thread incorrectly reuses the retired main
heap because `tld_main.thread_id` was not cleared.

The first run was restricted. Two subsequent runs of the saved mimalloc
binary also crashed at the same test transition with licenses configured.
For the user's supplied key, the runner asserted `symbolica.is_licensed()`
inside the test worker before each test action; it returned true. That run
still exited 139 (PID 1107697). No license key is included in this report.

With the same supplied key and the system allocator, the full Python suite
finished successfully: **28 passed, 2 standalone-only tests skipped**. It used
the same worker initialization and retirement pattern. This validates the
local host workaround; it does not establish a general mimalloc fix.

## Hardware qualification

After the user reported possible physical memory faults, Linux EDAC counters
on this machine showed **83,306 corrected errors**, all under
`mc1/csrow0/ch8_ce_count` (label `mc#1csrow#0channel#8`). Both controllers
reported zero uncorrected errors; controller `mc0` also reported zero corrected
errors. The counters covered about 17.9 days and did not change across two
reads approximately 30 seconds apart. System kernel-journal access was denied
to this account, so event timestamps could not be correlated with the crashes.

These are corrected-error reports, not proof of corrupted application data.
They establish a separate hardware reliability concern. Initially this left
the software attribution uncertain. The subsequent isolated reproducer,
watchpoint, upstream diagnosis, and exact backport comparison resolve the
specific crash: it is mimalloc bug #1287. Twenty fresh-process runs crashed
before the fix and twenty passed after it. The same result held with CPU and
memory bound to NUMA nodes on each CPU socket. No other hardware faults have
been ruled out by these checks.
