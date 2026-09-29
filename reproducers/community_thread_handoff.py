"""Control probe: load Symbolica on a worker, then use successive workers.

Run with the full community host's Python. This lightweight control passed
with mimalloc; the full OneLOop unittest suite reproduced a mimalloc
crash in _mi_subproc between test modules. Keep both checks when changing the
host allocator: passing this control alone does not establish that it is safe.
"""
import threading


def work():
    from symbolica import S

    x = S("thread_handoff_probe::x")
    assert (x + x).expand() == 2 * x


threading.stack_size(128 * 1024 * 1024)
for index in range(8):
    errors = []

    def run():
        try:
            work()
        except BaseException as error:
            errors.append(error)

    worker = threading.Thread(target=run)
    worker.start()
    worker.join()
    if errors:
        raise errors[0]
    print(f"Thread {index + 1}: OK", flush=True)
