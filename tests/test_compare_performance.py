"""Stdlib-only report validation tests; no integral evaluation or Symbolica import.

Run from the crate root with:
    python3 -m unittest discover -s tests -p test_compare_performance.py -v

The synthetic numbers test bookkeeping, not numerical parity or performance.
There are 31 fixtures, 36 SAME/HETERO workloads and 180 workload/batch groups.
"""

import copy
import csv
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import compare_performance as comparison


CALLS = 1031  # Partial tails for all nonunit survey batches and both family sizes.
REPETITIONS = 3


def synthetic_reports(ratio=1.25):
    """Build independent fixtures/checksums without invoking comparison helpers."""
    fixtures = []
    labels = ("euclidean", "timelike", "complex", "ir", "mu_0.01", "mu_10000")
    for family_index, family in enumerate(comparison.FAMILIES):
        for case_index, label in enumerate(labels + (("extra",) if family == "B0" else ())):
            finite = [(family_index + 1) * (case_index + 1) / 8, (case_index - 2) / 16]
            if family == "A0" and label == "ir":
                finite = [0.0, 0.0]
            expected = [finite, [float(family_index % 2), 0.0], [0.0, 0.0]]
            if family in ("C0", "D0") and label == "ir":
                expected[1:] = [[0.0, 0.5], [0.25, 0.0]]
            fixtures.append({
                "family": family, "name": f"{family}_{label}",
                "expected": expected, "actual": copy.deepcopy(expected),
                "normalization": (case_index + 1) / 4,
            })

    workloads = [("SAME", case["family"], case["name"], [case]) for case in fixtures]
    workloads += [("HETERO", family, "all", [case for case in fixtures if case["family"] == family])
                  for family in comparison.FAMILIES]
    reference = {
        "schema": "oneloop_fortran_benchmark_v1", "correctness": comparison.CORRECTNESS,
        "warmup": 256, "iterations": CALLS, "repetitions": REPETITIONS,
        "fixtures": fixtures, "summaries": [], "samples": [],
    }
    native = []
    for index, (mode, family, name, cases) in enumerate(workloads):
        # Sum the explicit input stream, independently of the production
        # cycles/tail checksum formula, including a nonzero complex component.
        checksum = [math.fsum(cases[row % len(cases)]["expected"][0][part]
                              for row in range(CALLS)) for part in (0, 1)]
        median = 100.0 + index * 4
        dimensions = {
            "mode": mode, "family": family, "name": name,
            "calls": CALLS, "distinct_rows": len(cases),
        }
        reference["summaries"].append({**dimensions, "wall_ns_per_call_median": median})
        # Deliberately not sorted: the comparator must calculate the median.
        for repetition, offset in enumerate((16.0, -8.0, 0.0), 1):
            timing = median + offset
            reference["samples"].append({
                **dimensions, "repetition": repetition,
                "wall_ns_per_call": timing, "finite_checksum": checksum[:],
            })
            for batch in comparison.BATCHES:
                native.append({
                    "mode": mode, "family": family, "name": name,
                    "batch": batch, "calls": CALLS, "repetition": repetition,
                    "ns_per_call": timing * ratio,
                    "checksum_re": checksum[0], "checksum_im": checksum[1],
                })
    return reference, native


def native_text(rows):
    output = io.StringIO()
    output.write("Synthetic license/log preamble, not a TSV record.\n")
    writer = csv.DictWriter(output, fieldnames=comparison.HEADER, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    output.write("Synthetic log epilogue.\n")
    return output.getvalue()


class ComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.valid_reference, cls.valid_native = synthetic_reports()

    def setUp(self):
        self.reference = copy.deepcopy(self.valid_reference)
        self.native = copy.deepcopy(self.valid_native)

    def compare(self, **options):
        return comparison.compare(self.reference, native_text(self.native), **options)

    def test_complete_inventory_and_independent_medians(self):
        result = self.compare()
        self.assertEqual(len(self.reference["fixtures"]), 31)
        self.assertEqual(len(self.reference["summaries"]), 36)
        self.assertEqual(len(result), 180)
        self.assertEqual(len({(row[0], row[1]) for row in result}), 180)
        self.assertEqual({row[1] for row in result}, {1, 4, 32, 256, 1024})
        self.assertTrue(all(row[-1] for row in result))
        for key, _, original, native, ratio, _ in result:
            summary = next(row for row in self.reference["summaries"]
                           if (row["mode"], row["family"], row["name"]) == key)
            self.assertEqual(original, summary["wall_ns_per_call_median"])
            self.assertEqual(native, original * 1.25)
            self.assertEqual(ratio, 1.25)

    def test_explicit_family_and_batch_selection(self):
        result = self.compare(family="B0", batch=4)
        self.assertEqual(len(result), 8)  # Seven fixtures and one HETERO workload.
        self.assertTrue(all(row[0][1] == "B0" and row[1] == 4 for row in result))
        self.assertEqual(len(self.compare(batch=32)), 36)

    def test_ratio_threshold_is_inclusive_and_workload_specific(self):
        self.assertTrue(all(row[-1] for row in self.compare(maximum_ratio=1.25)))
        self.assertTrue(all(not row[-1] for row in self.compare(maximum_ratio=1.2)))
        for row in self.native:
            if row["mode"] == "HETERO" and row["family"] == "D0" and row["batch"] == 1024:
                row["ns_per_call"] *= 1.6 / 1.25
        failures = [row for row in self.compare() if not row[-1]]
        self.assertEqual(len(failures), 1)
        self.assertEqual(failures[0][:2], (("HETERO", "D0", "all"), 1024))
        self.assertAlmostEqual(failures[0][4], 1.6)

    def test_missing_entire_native_batch_is_not_silently_ignored(self):
        self.native = [row for row in self.native if row["batch"] != 1024]
        with self.assertRaisesRegex(ValueError, "inventory"):
            self.compare()
        self.assertEqual(len(self.compare(batch=4)), 36)

    def test_native_missing_or_duplicate_repetitions(self):
        for mutation in ("missing", "duplicate"):
            with self.subTest(mutation=mutation):
                self.native = copy.deepcopy(self.valid_native)
                if mutation == "missing":
                    self.native.pop(0)
                else:
                    self.native.append(self.native[0].copy())
                with self.assertRaisesRegex(ValueError, "missing/duplicate native repetitions"):
                    self.compare()

    def test_original_missing_entire_repetition_or_duplicate_sample(self):
        for mutation in ("missing", "duplicate"):
            with self.subTest(mutation=mutation):
                self.reference = copy.deepcopy(self.valid_reference)
                if mutation == "missing":
                    self.reference["samples"] = [row for row in self.reference["samples"]
                                                 if row["repetition"] != 2]
                else:
                    self.reference["samples"].append(self.reference["samples"][0].copy())
                with self.assertRaisesRegex(ValueError, "missing/duplicate original repetitions"):
                    self.compare()

    def test_original_missing_workload_summary(self):
        self.reference["summaries"].pop()
        with self.assertRaisesRegex(ValueError, "incomplete original workload inventory"):
            self.compare()

    def test_original_summary_must_match_sample_median(self):
        self.reference["summaries"][0]["wall_ns_per_call_median"] *= 2
        with self.assertRaisesRegex(ValueError, "median disagrees"):
            self.compare()

    def test_both_finite_checksum_components_are_checked(self):
        for source in ("original", "native"):
            for part in (0, 1):
                with self.subTest(source=source, part=part):
                    self.reference = copy.deepcopy(self.valid_reference)
                    self.native = copy.deepcopy(self.valid_native)
                    if source == "original":
                        self.reference["samples"][0]["finite_checksum"][part] += 1
                    else:
                        self.native[0][("checksum_re", "checksum_im")[part]] += 1
                    with self.assertRaisesRegex(ValueError, f"{source}.*checksum mismatch"):
                        self.compare()

    def test_all_three_original_laurent_coefficients_are_checked(self):
        for coefficient in range(3):
            with self.subTest(coefficient=coefficient):
                self.reference = copy.deepcopy(self.valid_reference)
                self.reference["fixtures"][0]["actual"][coefficient][1] += 1
                with self.assertRaisesRegex(ValueError, "fixture correctness mismatch"):
                    self.compare()
        self.reference = copy.deepcopy(self.valid_reference)
        self.reference["fixtures"][0]["actual"].pop()
        with self.assertRaisesRegex(ValueError, "missing original Laurent coefficients"):
            self.compare()

    def test_original_warmup_must_be_exactly_256(self):
        for warmup in (None, 0, -1, 255, 257, "256", True):
            with self.subTest(warmup=warmup):
                self.reference = copy.deepcopy(self.valid_reference)
                if warmup is None:
                    del self.reference["warmup"]
                else:
                    self.reference["warmup"] = warmup
                with self.assertRaisesRegex(ValueError, "warmup.*256"):
                    self.compare()

    def test_original_schema_and_success_record_are_required(self):
        for field, message in (("schema", "schema"), ("correctness", "correctness record")):
            for value in (None, "failed"):
                with self.subTest(field=field, value=value):
                    self.reference = copy.deepcopy(self.valid_reference)
                    if value is None:
                        del self.reference[field]
                    else:
                        self.reference[field] = value
                    with self.assertRaisesRegex(ValueError, message):
                        self.compare()

    def test_nonfinite_and_nonpositive_timings_are_rejected(self):
        for value in (0.0, -1.0, float("nan"), float("inf")):
            with self.subTest(value=value):
                self.native[0]["ns_per_call"] = value
                with self.assertRaisesRegex(ValueError, "invalid native timing/result"):
                    self.compare()
        self.native = copy.deepcopy(self.valid_native)
        for source in ("fixtures", "samples"):
            with self.subTest(source=source):
                self.reference = copy.deepcopy(self.valid_reference)
                pair = (self.reference[source][0]["actual"][0] if source == "fixtures"
                        else self.reference[source][0]["finite_checksum"])
                pair[0] = float("nan")
                with self.assertRaisesRegex(ValueError, "finite complex pair"):
                    self.compare()

    def test_native_headers_and_row_dimensions_are_strict(self):
        header = "\t".join(comparison.HEADER)
        valid = native_text(self.native)
        for text in (
            valid.replace(header + "\n", "", 1),
            header + "\n" + valid,
            valid.replace(header, header + "\textra", 1),
            valid.replace("\nSAME\t", "\nSAME\textra\t", 1),
        ):
            with self.subTest(text=text[:80]):
                with self.assertRaisesRegex(ValueError, "header|malformed native TSV row"):
                    comparison.compare(self.reference, text)

    def test_call_counts_and_workload_names_must_match(self):
        for field, value, message in (
            ("calls", CALLS - 1, "different call counts"),
            ("name", "unmatched_fixture", "unmatched native workload"),
        ):
            with self.subTest(field=field):
                self.native = copy.deepcopy(self.valid_native)
                self.native[0][field] = value
                with self.assertRaisesRegex(ValueError, message):
                    self.compare()
        self.native = copy.deepcopy(self.valid_native)
        next(row for row in self.native if row["mode"] == "HETERO")["name"] = "not_all"
        with self.assertRaisesRegex(ValueError, "unexpected native heterogeneous workload name"):
            self.compare()

    def test_invalid_selection_and_ratio_arguments(self):
        for options in (
            {"batch": 0}, {"batch": -1}, {"batch": True}, {"family": "E0"},
            {"maximum_ratio": 0}, {"maximum_ratio": -1},
            {"maximum_ratio": float("nan")}, {"maximum_ratio": float("inf")},
        ):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.compare(**options)

    def test_cli_exit_status_and_no_partial_pass_output_for_invalid_data(self):
        script = Path(comparison.__file__).resolve()
        with tempfile.TemporaryDirectory(prefix="oneloop-compare-unit-") as directory:
            reference_path = Path(directory, "fortran.json")
            native_path = Path(directory, "native.tsv")
            reference_path.write_text(json.dumps(self.reference))
            native_path.write_text(native_text(self.native))
            command = [sys.executable, str(script), str(reference_path), str(native_path)]
            passed = subprocess.run(command, capture_output=True, text=True, check=False, timeout=30)
            self.assertEqual(passed.returncode, 0, passed.stderr)
            self.assertEqual(len(passed.stdout.splitlines()), 181)
            self.assertEqual(passed.stdout.count("\tPASS\n"), 180)
            slow = subprocess.run(command + ["--maximum-ratio", "1.2"], capture_output=True,
                                  text=True, check=False, timeout=30)
            self.assertEqual(slow.returncode, 2, slow.stderr)
            self.assertEqual(slow.stdout.count("\tFAIL\n"), 180)
            # Reject a malformed final group before printing any earlier PASS.
            native_path.write_text(native_text(self.native[:-1]))
            invalid = subprocess.run(command, capture_output=True, text=True,
                                     check=False, timeout=30)
            self.assertEqual(invalid.returncode, 1)
            self.assertEqual(invalid.stdout, "")
            self.assertIn("compare_performance:", invalid.stderr)


if __name__ == "__main__":
    unittest.main()
