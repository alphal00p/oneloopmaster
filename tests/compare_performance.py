#!/usr/bin/env python3
"""Compare paired Fortran JSON and native performance_survey TSV; stdlib only.

No result is a universal-domain guarantee. Report each matched workload/batch;
never hide a slow family in an aggregate average. With no --batch selection,
require the current survey's complete 1/4/32/256/1024 batch inventory. Exit 2 if
any selected median exceeds the requested ratio (default 1.5); exit 1 on
incomplete/mismatched data. Finite checksums are validated, but are not a
replacement for the survey's untimed checks of all three Laurent coefficients.
"""
import argparse
import csv
import io
import json
import math
from pathlib import Path
import statistics
import sys


FAMILIES = ("A0", "B0", "dB0", "C0", "D0")
BATCHES = (1, 4, 32, 256, 1024)
HEADER = (
    "mode", "family", "name", "batch", "calls", "repetition", "ns_per_call",
    "checksum_re", "checksum_im",
)
CORRECTNESS = "all fixture outputs, batch-final outputs and finite checksums passed"


def positive_integer(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def workload_key(row):
    if row["mode"] not in ("SAME", "HETERO") or row["family"] not in FAMILIES:
        raise ValueError("invalid workload mode/family")
    return (row["mode"], row["family"], row["name"] if row["mode"] == "SAME" else "all")


def finite_complex(value):
    if len(value) != 2 or not all(math.isfinite(v) for v in value):
        raise ValueError("expected a finite complex pair")
    return complex(*value)


def checksum_contract(cases, calls):
    cycles, tail = divmod(calls, len(cases))
    weighted = [(case, cycles + int(index < tail)) for index, case in enumerate(cases)]
    expected = complex(*(
        math.fsum(count * case["expected"][0][part] for case, count in weighted)
        for part in (0, 1)
    ))
    absolute_sum = math.fsum(
        count * abs(finite_complex(case["expected"][0])) for case, count in weighted
    )
    # Match the per-coefficient fixture tolerance, plus summation roundoff.
    tolerance = math.fsum(
        count * (2e-10 / case["normalization"] + 2e-8 * abs(finite_complex(case["expected"][0])))
        for case, count in weighted
    ) + max(1e-12, 1e-9 * absolute_sum)
    if not all(map(math.isfinite, (expected.real, expected.imag, tolerance))):
        raise ValueError("nonfinite expected workload checksum")
    return expected, tolerance


def baseline_data(reference, family):
    if reference.get("schema") != "oneloop_fortran_benchmark_v1":
        raise ValueError("unsupported original-backend report schema")
    if reference.get("correctness") != CORRECTNESS:
        raise ValueError("missing successful original-backend correctness record")
    if reference.get("warmup") != 256:
        raise ValueError("original warmup must match the native survey's 256 calls")
    calls = positive_integer(reference["iterations"], "original calls")
    repetitions = positive_integer(reference["repetitions"], "original repetitions")
    cases = {}
    for case in reference["fixtures"]:
        if family and case["family"] != family:
            continue
        if case["family"] not in FAMILIES:
            raise ValueError("invalid original fixture family")
        normalization = case["normalization"]
        if not math.isfinite(normalization) or normalization <= 0:
            raise ValueError("invalid fixture normalization")
        if len(case["expected"]) != 3 or len(case["actual"]) != 3:
            raise ValueError("missing original Laurent coefficients")
        for expected, actual in zip(case["expected"], case["actual"]):
            expected, actual = finite_complex(expected), finite_complex(actual)
            if abs(actual - expected) * normalization > 2e-10 + 2e-8 * abs(expected) * normalization:
                raise ValueError("original fixture correctness mismatch")
        key = ("SAME", case["family"], case["name"])
        if key in cases:
            raise ValueError("duplicate original fixture name")
        cases[key] = [case]
        cases.setdefault(("HETERO", case["family"], "all"), []).append(case)
    if not cases:
        raise ValueError("empty original fixture selection")

    baseline = {}
    for row in reference["summaries"]:
        if family and row["family"] != family:
            continue
        key = workload_key(row)
        if key in baseline or key not in cases:
            raise ValueError("duplicate/unmatched reference workload")
        if row["calls"] != calls or row["distinct_rows"] != len(cases[key]):
            raise ValueError("incorrect original workload dimensions")
        original = row["wall_ns_per_call_median"]
        if not math.isfinite(original) or original <= 0:
            raise ValueError("invalid original median")
        baseline[key] = row
    if set(baseline) != set(cases):
        raise ValueError("incomplete original workload inventory")

    checksums = {key: checksum_contract(rows, calls) for key, rows in cases.items()}
    samples = {key: [] for key in cases}
    for row in reference["samples"]:
        if family and row["family"] != family:
            continue
        key = workload_key(row)
        if key not in samples or row["calls"] != calls or row["distinct_rows"] != len(cases[key]):
            raise ValueError("unmatched/malformed original sample")
        timing = row["wall_ns_per_call"]
        if not math.isfinite(timing) or timing <= 0:
            raise ValueError("invalid original sample timing")
        expected, tolerance = checksums[key]
        if abs(finite_complex(row["finite_checksum"]) - expected) > tolerance:
            raise ValueError("original sample checksum mismatch")
        samples[key].append((positive_integer(row["repetition"], "original repetition"), timing))
    expected_repetitions = list(range(1, repetitions + 1))
    for key, values in samples.items():
        if sorted(index for index, _ in values) != expected_repetitions:
            raise ValueError(f"missing/duplicate original repetitions: {key}")
        computed = statistics.median(value for _, value in values)
        if not math.isclose(computed, baseline[key]["wall_ns_per_call_median"], rel_tol=1e-12):
            raise ValueError(f"original summary median disagrees with samples: {key}")
    return baseline, checksums, calls, expected_repetitions


def compare(reference, native_text, *, batch=None, family=None, maximum_ratio=1.5):
    if batch is not None:
        positive_integer(batch, "batch")
    if family is not None and family not in FAMILIES:
        raise ValueError("unknown scalar family")
    if not math.isfinite(maximum_ratio) or maximum_ratio <= 0:
        raise ValueError("maximum ratio must be finite and positive")
    baseline, checksums, calls, repetitions = baseline_data(reference, family)
    # Permit library license/log messages around the explicit TSV records.
    lines = [line for line in native_text.splitlines()
             if line.startswith(("mode\tfamily\t", "SAME\t", "HETERO\t"))]
    header = "\t".join(HEADER)
    if not lines or lines[0] != header or sum(line.startswith("mode\tfamily\t") for line in lines) != 1:
        raise ValueError("missing, malformed or duplicate native TSV header")
    groups = {}
    for row in csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"):
        if None in row or any(value is None for value in row.values()):
            raise ValueError("malformed native TSV row")
        row_batch = positive_integer(int(row["batch"]), "native batch")
        if batch is not None and row_batch != batch:
            continue
        if family and row["family"] != family:
            continue
        key = workload_key(row)
        if row["mode"] == "HETERO" and row["name"] != "all":
            raise ValueError("unexpected native heterogeneous workload name")
        if key not in baseline:
            raise ValueError(f"unmatched native workload: {key}")
        if int(row["calls"]) != calls:
            raise ValueError(f"different call counts: {key}")
        numeric = [float(row[field]) for field in ("ns_per_call", "checksum_re", "checksum_im")]
        if not all(map(math.isfinite, numeric)) or numeric[0] <= 0:
            raise ValueError(f"invalid native timing/result: {key}")
        expected, tolerance = checksums[key]
        if abs(complex(*numeric[1:]) - expected) > tolerance:
            raise ValueError(f"native checksum mismatch: {key}, batch={row_batch}")
        repetition = positive_integer(int(row["repetition"]), "native repetition")
        groups.setdefault((key, row_batch), []).append((repetition, numeric[0]))
    batches = (batch,) if batch is not None else BATCHES
    if set(groups) != {(key, size) for key in baseline for size in batches}:
        raise ValueError("incomplete/unexpected native family/workload/batch inventory")
    result = []
    for (key, size), values in groups.items():
        if sorted(index for index, _ in values) != repetitions:
            raise ValueError(f"missing/duplicate native repetitions: {key}, batch={size}")
        original = baseline[key]["wall_ns_per_call_median"]
        native = statistics.median(value for _, value in values)
        ratio = native / original
        result.append((key, size, original, native, ratio, ratio <= maximum_ratio))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fortran", type=Path)
    parser.add_argument("native", type=Path)
    parser.add_argument("--batch", type=int, help="select one positive native batch size")
    parser.add_argument("--family", choices=FAMILIES)
    parser.add_argument("--maximum-ratio", type=float, default=1.5)
    args = parser.parse_args()
    result = compare(
        json.loads(args.fortran.read_text()), args.native.read_text(),
        batch=args.batch, family=args.family, maximum_ratio=args.maximum_ratio,
    )
    print("mode\tfamily\tname\tbatch\tfortran_ns\tnative_ns\tratio\tstatus")
    for key, batch, original, native, ratio, passed in result:
        status = "PASS" if passed else "FAIL"
        print("\t".join(key) + f"\t{batch}\t{original:.3f}\t{native:.3f}\t{ratio:.4f}\t{status}")
    return 0 if all(row[-1] for row in result) else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, ValueError, TypeError, OverflowError, OSError) as error:
        print(f"compare_performance: {error}", file=sys.stderr)
        raise SystemExit(1) from error
