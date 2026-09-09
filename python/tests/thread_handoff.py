"""Standalone diagnostic: use eager caches after their initializing thread exits.

Run separately from the test suite, with the native extension on PYTHONPATH.
Both workers have explicit 128 MiB stacks and never overlap.
"""
import importlib

from test_api import family_selector
import threading


def run(stage, action):
    errors = []

    def work():
        print(f"{stage}: entered worker", flush=True)
        try:
            module = importlib.import_module("oneloop_native")
            print(f"{stage}: imported", flush=True)
            action(module)
            print(f"{stage}: completed", flush=True)
        except BaseException as error:
            errors.append(error)

    worker = threading.Thread(target=work)
    worker.start()
    worker.join()
    if errors:
        raise errors[0]
    print(f"{stage}: joined", flush=True)


def clone_all(module):
    for family in ("A0", "B0", "dB0", "C0", "D0"):
        print(f"constructing {family}", flush=True)
        evaluator = module.Evaluator(family_selector(module, family))
        print(repr(evaluator), flush=True)
        del evaluator


if __name__ == "__main__":
    threading.stack_size(128 * 1024 * 1024)
    run("first", lambda module: print(module.A0(2), flush=True))
    run("second", clone_all)
