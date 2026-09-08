"""Pure-stdlib checks of Python survey bookkeeping using fake numeric evaluators."""

import copy
import csv
import importlib.util
import io
import json
import math
from pathlib import Path
import tempfile
import types
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location("python_performance_survey", Path(__file__).resolve().parents[1] / "performance_survey.py")
survey = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(survey)


class FakeEvaluator:
    def __init__(self, cases):
        self.values = {tuple(survey.arguments(case)): tuple(complex(*value) for value in case["expected"])
                       for case in cases}
        self.singles = 0
        self.batch_sizes = []

    def evaluate(self, arguments):
        self.singles += 1
        return self.values[tuple(arguments)]

    def evaluate_batch(self, rows):
        self.batch_sizes.append(len(rows))
        return [self.values[tuple(row)] for row in rows]


class SurveyTests(unittest.TestCase):
    def setUp(self):
        self.cases = survey.support.read_fixtures(survey.support.DEFAULT_FIXTURES, ["B0"])

    def test_warmup_exact_calls_tails_and_all_checksum_outputs(self):
        evaluator = FakeEvaluator(self.cases)
        samples, summary = survey.workload(evaluator, self.cases, calls=11, repetitions=2,
                                           batch=4, mode="HETERO", name="all")
        self.assertEqual(evaluator.singles, 256)
        self.assertEqual(evaluator.batch_sizes, [4, 4, 3] * 3)  # Validation plus two measurements.
        self.assertEqual(summary["calls"], 11)
        self.assertEqual(summary["tail_rows"], 3)
        self.assertEqual(summary["batch_calls"], 3)
        self.assertEqual(summary["distinct_rows"], 7)
        self.assertEqual(len(samples), 2)
        expected = [[math.fsum(self.cases[index % 7]["expected"][coefficient][part]
                              for index in range(11)) for part in range(2)] for coefficient in range(3)]
        for sample in samples:
            self.assertEqual(sample["laurent_checksums"], expected)
            self.assertEqual(sample["finite_checksum"], expected[0])
            self.assertGreater(sample["ns_per_call"], 0)
        rows = list(csv.DictReader(io.StringIO(survey.as_tsv(samples)), delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertEqual(tuple(rows[0]), survey.comparison.HEADER)
        self.assertEqual(rows[0]["batch"], "4")

    def test_bad_pole_and_missing_tail_are_rejected(self):
        evaluator = FakeEvaluator(self.cases)
        output = list(evaluator.evaluate(survey.arguments(self.cases[0])))
        output[2] += 1j
        with self.assertRaisesRegex(ValueError, "coefficient 2 disagrees"):
            survey.checked_output(self.cases[0], output)
        stream = [evaluator.evaluate(survey.arguments(self.cases[index % 7])) for index in range(11)]
        with self.assertRaisesRegex(ValueError, "incorrect output row count"):
            survey.check_batches(self.cases, [stream[:4], stream[4:8], stream[8:10]], 11, 4)

    def test_fortran_comparison_matches_inputs_not_only_names(self):
        reference = {"fixtures": copy.deepcopy(self.cases)}
        survey.matching_reference(reference, self.cases)
        reference["fixtures"][0]["mu_squared"] *= 4
        with self.assertRaisesRegex(ValueError, "input/output fixture does not match"):
            survey.matching_reference(reference, self.cases)

    def test_full_cli_control_flow_with_fake_module_only(self):
        fake_module = types.SimpleNamespace(
            is_initialized=lambda: True,
            SYMBOLICA_REVISION="synthetic-test-only",
            Evaluator=lambda family: FakeEvaluator([case for case in self.cases if case["family"] == family]),
        )
        with tempfile.TemporaryDirectory(prefix="oneloop-python-survey-unit-") as directory:
            output, tsv = Path(directory, "report.json"), Path(directory, "report.tsv")
            argv = [str(SPEC.origin), "--module", "fake_oneloop", "--family", "B0",
                    "--iterations", "7", "--repetitions", "2", "--batch", "4",
                    "--output", str(output), "--tsv", str(tsv)]
            with mock.patch.object(survey.sys, "argv", argv), \
                    mock.patch.object(survey.importlib, "import_module", return_value=fake_module) as importer, \
                    mock.patch.object(survey.sys, "stderr", io.StringIO()):
                self.assertEqual(survey.main(), 0)
            importer.assert_called_once_with("fake_oneloop")
            report = json.loads(output.read_text())
            self.assertEqual(report["schema"], "oneloop_python_benchmark_v1")
            self.assertEqual(report["warmup"], 256)
            self.assertEqual(report["batch"], 4)
            self.assertEqual(len(report["fixtures"]), 7)
            self.assertEqual(len(report["summaries"]), 8)
            self.assertEqual(len(report["samples"]), 16)
            self.assertEqual(len(list(csv.DictReader(io.StringIO(tsv.read_text()), delimiter="\t"))), 16)
            self.assertTrue(report["timing_contract"]["python_boundary_conversion_included"])
            self.assertTrue(report["timing_contract"]["output_validation_and_checksums_excluded"])


if __name__ == "__main__":
    unittest.main()
