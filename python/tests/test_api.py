"""Optional adapter tests: python -m unittest discover -s python/tests -v.

Build/install the local adapter first. Set ONELOOP_PYTHON_MODULE to the actual
shared-kernel module for community integration tests. ONELOOP_PYTHON_REBUILD=1
passes explicit rebuild=True, backend="symjit" to numeric calls without disabling
eager startup. Native workspace rebuilding is tested independently as well.
Cold-start subprocess tests run first, sequentially on their main threads with
at most an 8 MiB stack on POSIX. Remaining Symbolica work runs on one explicitly
created, persistent 128 MiB test thread. Never run this suite concurrently with
another Symbolica process using a restricted license. The adapter spawns no thread.
"""

import importlib
import json
import math
import os
from pathlib import Path
import queue
import subprocess
import sys
import tempfile
import threading
import unittest


MODULE = os.environ.get("ONELOOP_PYTHON_MODULE", "oneloop_native")
REBUILD = os.environ.get("ONELOOP_PYTHON_REBUILD", "0") == "1"
FIXTURES = Path(__file__).resolve().parents[2] / "tests/data/benchmark.txt"
FAMILIES = {1: "A0", 2: "B0", -2: "dB0", 3: "C0", 4: "D0"}


def read_fixture_groups():
    """Decode numeric data only; importing this module never loads Symbolica."""
    groups = {kind: [] for kind in FAMILIES}
    for line_number, line in enumerate(FIXTURES.read_text().splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        values = list(map(float, line.split()))
        kind = int(values[0])
        np = 0 if kind == 1 else 1 if abs(kind) == 2 else 3 if kind == 3 else 6
        nm = abs(kind)
        inputs = values[2:2 + np]
        inputs += [complex(*values[i:i + 2]) for i in range(2 + np, 2 + np + 2 * nm, 2)]
        inputs += [values[1]]
        expected = [complex(*values[i:i + 2]) for i in range(len(values) - 6, len(values), 2)]
        scale = max(map(abs, values[2:-6]), default=1) or 1
        normalization = scale ** {1: -1, 2: 0, -2: 1, 3: 1, 4: 2}[kind]
        groups[kind].append((line_number, inputs, expected, normalization))
    return groups


def new_evaluator(module, family, **options):
    options.setdefault("rebuild", REBUILD)
    if REBUILD:
        options.setdefault("backend", "symjit")
    return module.Evaluator(family, **options)


def scalar(module, name, *arguments, **options):
    options.setdefault("rebuild", REBUILD)
    if REBUILD:
        options.setdefault("backend", "symjit")
    return getattr(module, name)(*arguments, **options)


def on_fresh_symbolica_thread(action):
    """Opt-in diagnostic for work after the initializing thread has exited."""
    result = []
    errors = []

    def work():
        try:
            result.append(action(importlib.import_module(MODULE)))
        except BaseException as error:
            errors.append(error)

    previous = threading.stack_size()
    threading.stack_size(128 * 1024 * 1024)
    try:
        worker = threading.Thread(target=work)
        worker.start()
        worker.join()
    finally:
        threading.stack_size(previous)
    if errors:
        raise errors[0]
    return result[0]


class SymbolicaTestWorker:
    """Keep eager initialization and all subsequent test work on the same thread."""

    def __init__(self):
        self.requests = queue.Queue()
        previous = threading.stack_size()
        threading.stack_size(128 * 1024 * 1024)
        try:
            self.worker = threading.Thread(target=self.work, name="oneloop-python-tests")
            self.worker.start()
        finally:
            threading.stack_size(previous)

    def work(self):
        while True:
            request = self.requests.get()
            if request is None:
                return
            action, response = request
            try:
                response.put((True, action(importlib.import_module(MODULE))))
            except BaseException as error:
                response.put((False, error))
            finally:
                del request, action, response

    def call(self, action):
        response = queue.Queue(maxsize=1)
        self.requests.put((action, response))
        succeeded, value = response.get()
        if not succeeded:
            raise value
        return value

    def close(self):
        self.requests.put(None)
        self.worker.join()


WORKER = None


def on_symbolica_thread(action):
    global WORKER
    if os.environ.get("ONELOOP_PYTHON_THREAD_HANDOFF") == "1":
        # Separate thread-handoff diagnostic. This reproduced a native allocator
        # segfault before evaluator construction in the earlier debug build.
        return on_fresh_symbolica_thread(action)
    if WORKER is None:
        WORKER = SymbolicaTestWorker()
    return WORKER.call(action)


def tearDownModule():
    global WORKER
    if WORKER is not None:
        WORKER.close()
        WORKER = None


class AdapterTests(unittest.TestCase):
    def run_standalone_child(self, family, entry, *, guard=None, backend="auto"):
        environment = os.environ.copy()
        environment.pop("SYMJIT_TOML", None)
        if "PYTHONPATH" in environment:
            environment["PYTHONPATH"] = os.pathsep.join(
                str(Path(path or ".").resolve())
                for path in environment["PYTHONPATH"].split(os.pathsep)
            )
        command = [
            sys.executable, "-X", "faulthandler",
            str(Path(__file__).with_name("cold_start.py").resolve()), family, entry,
            "--backend", backend,
        ]
        if REBUILD:
            command.append("--rebuild-only")
        with tempfile.TemporaryDirectory(prefix="oneloop-python-cold-") as directory:
            if guard == "environment":
                environment["SYMJIT_TOML"] = ""
            elif guard == "file":
                Path(directory, "symjit.toml").write_text("# Disallowed portable-API override.\n")
            completed = subprocess.run(
                command, cwd=directory, env=environment, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300,
            )
        self.assertEqual(
            completed.returncode, 0,
            f"{family}/{entry} child failed ({completed.returncode})\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
        )
        records = [line.removeprefix("ONELOOP_COLD_JSON ")
                   for line in completed.stdout.splitlines() if line.startswith("ONELOOP_COLD_JSON ")]
        self.assertEqual(len(records), 1, completed.stdout)
        return json.loads(records[0])

    @unittest.skipIf(REBUILD, "cold prebuilt loading requires embedded evaluator assets")
    @unittest.skipUnless(MODULE == "oneloop_native", "standalone extension cold-load test")
    def test_00_main_thread_cold_start(self):
        # Run before this process imports the extension in any worker test.
        # Eager import already loads all five families. Two fresh processes
        # distinguish first-call entry routes without ten redundant startups.
        for entry in ("direct", "evaluator"):
            with self.subTest(entry=entry):
                backend = "auto" if entry == "direct" else "symjit"
                record = self.run_standalone_child("D0", entry, backend=backend)
                self.assertEqual(record["backend_requested"], backend)
                if backend == "symjit":
                    self.assertEqual(record["backend_resolved"], "symjit")
                self.assertTrue(record["main_thread"])
                self.assertEqual(record["family"], "D0")
                self.assertEqual(record["entry"], entry)
                self.assertTrue(record["import_completed_before_master_calls"])
                self.assertTrue(record["initialized_before_master_calls"])
                self.assertEqual(record["ready_families"], list(FAMILIES.values()))
                self.assertEqual(record["checked_rows"], 31)
                self.assertEqual(record["batch_rows"], 45)
                self.assertEqual(
                    record["checks_by_family"],
                    {name: {"single_rows": 7 if name == "B0" else 6, "batch_rows": 9}
                     for name in FAMILIES.values()},
                )
                if record["stack_soft_limit_bytes"] is not None:
                    self.assertLessEqual(record["stack_soft_limit_bytes"], 8 * 1024 * 1024)

    @unittest.skipUnless(MODULE == "oneloop_native", "standalone extension environment guards")
    def test_01_portable_environment_guards(self):
        for entry, guard in (("guard-env", "environment"), ("guard-file", "file")):
            with self.subTest(guard=guard):
                record = self.run_standalone_child("D0", entry, guard=guard)
                self.assertEqual(record["guarded_imports"], 1)

    def test_elementary_calls_and_order(self):
        def check(module):
            self.assertTrue(module.is_initialized())
            self.assertEqual(module.COEFFICIENT_ORDER, (0, -1, -2))
            self.assertEqual(scalar(module, "A0", 0, mu_squared=1), (0j, 0j, 0j))
            result = scalar(module, "A0", 2, mu_squared=1)
            self.assertTrue(all(isinstance(value, complex) for value in result))
            self.assertAlmostEqual(result[0].real, 2 * (1 - math.log(2)), places=12)
            self.assertEqual(result[1:], (2 + 0j, 0j))
            self.assertEqual(scalar(module, "a0", 2, mu_squared=1), result)
            self.assertEqual(scalar(module, "A0", 2 - 0.1j)[1], 2 - 0.1j)
            scaled = scalar(module, "A0", 2, mu_squared=4)
            self.assertAlmostEqual((scaled[0] - result[0]).real, 2 * math.log(4), places=12)
        on_symbolica_thread(check)

    def test_reusable_single_batch_and_rebuild(self):
        def check(module):
            evaluator = new_evaluator(module, "B0")
            self.assertEqual(evaluator.family, "B0")
            self.assertEqual(evaluator.arity, 4)
            rows = [
                [-1, 1, 1, 1],
                [7, 1, 2, 1],
                [3, 0.7 - 0.03j, 1.4 - 0.08j, 4],
                [2, 0, 2, 1],
                [0, 0, 0, 1],
                [7, 1, 2, 0.01],
                [7, 1, 2, 10000],
                [-3, 1, 1, 4],
            ]
            singles = [evaluator.evaluate(row) for row in rows]
            self.assertEqual(evaluator.evaluate_batch([]), [])
            # Distinct lane predicates, a complete four-wide group, a tail, and
            # two full groups. These sizes do not assert that SIMD is enabled.
            layouts = (
                [0, 1, 2],
                [0, 1, 2, 3],
                [4, 2, 0, 6, 1],
                [0, 2, 5, 4, 7, 1, 3, 6],
            )
            for layout in layouts:
                actual_rows = evaluator.evaluate_batch([rows[index] for index in layout])
                self.assertEqual(len(actual_rows), len(layout))
                for actual, index in zip(actual_rows, layout):
                    self.assertEqual(len(actual), 3)
                    self.assertEqual(len(singles[index]), 3)
                    for coefficient, (a, b) in enumerate(zip(actual, singles[index])):
                        with self.subTest(rows=len(layout), input=index, coefficient=coefficient):
                            self.assertTrue(math.isfinite(abs(a - b)))
                            self.assertLess(abs(a - b), 1e-12)
            self.assertAlmostEqual(singles[0][0].real, -0.1520447048200202, places=12)
            # Explicitly distinguish SymJIT source rebuilding from a fresh
            # Native constant/workspace setup, independent of DEFAULT_BACKEND.
            for backend in ("symjit", "native"):
                rebuilt = new_evaluator(module, "B0", rebuild=True, backend=backend)
                for actual, expected in zip(rebuilt.evaluate(rows[0]), singles[0]):
                    self.assertLess(abs(actual - expected), 1e-12)
                rebuilt.rebuild()
                for actual, expected in zip(rebuilt.evaluate(rows[0]), singles[0]):
                    self.assertLess(abs(actual - expected), 1e-12)
                for actual, expected in zip(module.B0(-1, 1, 1, rebuild=True, backend=backend), singles[0]):
                    self.assertLess(abs(actual - expected), 1e-12)
        on_symbolica_thread(check)

    def test_every_shared_benchmark_fixture(self):
        def check(module):
            families = FAMILIES
            groups = read_fixture_groups()

            def check_output(actual_row, case):
                _, _, expected, normalization = case
                self.assertEqual(len(actual_row), 3)
                for actual, reference in zip(actual_row, expected):
                    error = abs(actual - reference) * normalization
                    self.assertTrue(math.isfinite(error))
                    self.assertLessEqual(error, 2e-10 + 2e-8 * abs(reference) * normalization)

            self.assertEqual(sum(map(len, groups.values())), 31)
            self.assertEqual(
                {kind: len(cases) for kind, cases in groups.items()},
                {1: 6, 2: 7, -2: 6, 3: 6, 4: 6},
            )
            for kind, cases in groups.items():
                evaluator = new_evaluator(module, families[kind])
                for case in cases:
                    with self.subTest(family=families[kind], line=case[0], mode="single"):
                        check_output(evaluator.evaluate(case[1]), case)
                # Cycle each six-/seven-fixture family to nine rows. Every
                # fixture is present; batch sizes 4/5/8 all have a partial tail,
                # and size 8 includes two complete four-wide groups.
                ordered = [cases[index % len(cases)] for index in range(9)]
                for batch_size in (4, 5, 8):
                    for offset in range(0, len(ordered), batch_size):
                        chunk = ordered[offset:offset + batch_size]
                        actual_rows = evaluator.evaluate_batch([case[1] for case in chunk])
                        self.assertEqual(len(actual_rows), len(chunk))
                        for actual, case in zip(actual_rows, chunk):
                            with self.subTest(
                                family=families[kind], line=case[0],
                                batch=batch_size, offset=offset,
                            ):
                                check_output(actual, case)
        on_symbolica_thread(check)

    def test_invalid_inputs(self):
        def check(module):
            with self.assertRaises(ValueError):
                module.Evaluator("not_an_integral")
            with self.assertRaises(ValueError):
                module.A0(1, mu_squared=0)
            with self.assertRaises(ValueError):
                module.A0(1 + 0.1j)
            with self.assertRaises(ValueError):
                module.B0(1 + 1j, 1, 2)
            with self.assertRaises(ValueError):
                module.A0(float("nan"))
            with self.assertRaises(TypeError):
                module.A0("not a number")
            evaluator = new_evaluator(module, "A0")
            with self.assertRaises(ValueError):
                evaluator.evaluate([1])
            with self.assertRaises(ValueError):
                evaluator.evaluate([1, 1 + 1j])
            with self.assertRaises(ValueError):
                evaluator.evaluate_batch([[1, 1], [2]])
        on_symbolica_thread(check)

    def test_expression_interop_only_in_shared_kernel(self):
        def check(module):
            self.assertTrue(module.is_initialized())
            if not module.EXPRESSION_INTEROP:
                self.assertFalse(hasattr(module, "master_coefficients"))
                self.assertFalse(hasattr(module, "compile_native"))
                return
            symbolica = importlib.import_module("symbolica")
            x = symbolica.S("oneloop_python_test_x")
            coefficients = module.master_coefficients("A0", [x, 1])
            self.assertEqual(len(coefficients), 3)
            evaluator = module.compile_native([coefficients[0] + coefficients[1]], [x])
            result = evaluator.evaluate_complex([2 + 0j])
            value = complex(result.reshape(-1)[0])
            self.assertLess(abs(value - (2 * (1 - math.log(2)) + 2)), 1e-12)
        on_symbolica_thread(check)


if __name__ == "__main__":
    unittest.main()
