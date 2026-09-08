#!/usr/bin/env python3
"""Warmed Python batch throughput; stdlib only, default batch size 1024.

Build/install the standalone extension first, preferably in release mode:
  python python/performance_survey.py --output /tmp/python.json --tsv /tmp/python.tsv

All three Laurent outputs are computed. The timed loop includes Python boundary
conversion and result allocation, but excludes input-list preparation, startup,
validation, checksums, and final output destruction. This differs from the Rust
and Fortran surveys, which include a native finite-coefficient checksum inside
their timed loops. Optional comparisons report that distinction, not a universal
performance guarantee. No worker thread or precision fallback is introduced.
"""

import argparse
import csv
import gc
import hashlib
import importlib
import importlib.util
import io
import json
import math
from pathlib import Path
import platform
import statistics
import sys
import threading
import time


ROOT = Path(__file__).resolve().parent.parent
WARMUP = 256


def load_support(name):
    """Load only the repository's pure-stdlib report/fixture helpers."""
    spec = importlib.util.spec_from_file_location(name, ROOT / "tests" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


support = load_support("benchmark_fortran")
comparison = load_support("compare_performance")


def arguments(case):
    return [*case["momenta"], *(complex(*mass) for mass in case["masses"]), case["mu_squared"]]


def checked_output(case, output):
    if len(output) != 3 or not all(isinstance(value, complex) for value in output):
        raise ValueError("expected exactly three Python complex Laurent coefficients")
    pairs = [[value.real, value.imag] for value in output]
    support.validate_output(case, pairs, case["name"])
    return pairs


def check_batches(cases, outputs, calls, batch):
    if len(outputs) != (calls + batch - 1) // batch:
        raise ValueError("incorrect output batch count")
    offset = 0
    for values in outputs:
        if len(values) != min(batch, calls - offset):
            raise ValueError("incorrect output row count, including partial tail")
        for output in values:
            checked_output(cases[offset % len(cases)], output)
            offset += 1
    if offset != calls:
        raise ValueError("incorrect total output row count")


def workload(evaluator, cases, *, calls, repetitions, batch, mode, name):
    prepared_rows = [arguments(case) for case in cases]
    stream = [prepared_rows[index % len(cases)] for index in range(calls)]
    chunks = [stream[offset:offset + batch] for offset in range(0, calls, batch)]
    # Check the exact mixed-lane layout and partial last batch outside timing.
    validation = [evaluator.evaluate_batch(chunk) for chunk in chunks]
    check_batches(cases, validation, calls, batch)
    del validation
    for index in range(WARMUP):
        evaluator.evaluate(prepared_rows[index % len(cases)])

    samples = []
    for repetition in range(1, repetitions + 1):
        outputs = [None] * len(chunks)
        started = time.perf_counter_ns()
        for index, chunk in enumerate(chunks):
            outputs[index] = evaluator.evaluate_batch(chunk)
        elapsed = time.perf_counter_ns() - started
        if elapsed <= 0:
            raise ValueError("nonpositive elapsed time; increase --iterations")
        check_batches(cases, outputs, calls, batch)
        checksums = [[math.fsum(getattr(row[coefficient], part)
                                for values in outputs for row in values)
                      for part in ("real", "imag")] for coefficient in range(3)]
        if not all(math.isfinite(value) for pair in checksums for value in pair):
            raise ValueError("nonfinite Laurent checksum")
        samples.append({
            "mode": mode, "family": cases[0]["family"], "name": name,
            "batch": batch, "calls": calls, "distinct_rows": len(cases),
            "batch_calls": len(chunks), "tail_rows": calls % batch,
            "repetition": repetition, "wall_ns": elapsed,
            "ns_per_call": elapsed / calls, "finite_checksum": checksums[0],
            "laurent_checksums": checksums,
        })
        del outputs  # Do not charge previous-repetition output destruction to the next.
    timing = [row["ns_per_call"] for row in samples]
    summary = {
        key: samples[0][key]
        for key in ("mode", "family", "name", "batch", "calls", "distinct_rows", "batch_calls", "tail_rows")
    }
    summary.update({
        "wall_ns_per_call_median": statistics.median(timing),
        "wall_ns_per_call_min": min(timing), "wall_ns_per_call_max": max(timing),
    })
    return samples, summary


def as_tsv(samples):
    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(comparison.HEADER)
    for row in samples:
        writer.writerow([*(row[key] for key in comparison.HEADER[:6]),
                         row["ns_per_call"], *row["finite_checksum"]])
    return output.getvalue()


def matching_reference(reference, fixtures):
    originals = {(case["family"], case["name"]): case for case in reference["fixtures"]}
    for case in fixtures:
        original = originals.get((case["family"], case["name"]))
        if original is None or any(original[key] != case[key] for key in (
            "momenta", "masses", "mu_squared", "expected", "normalization",
        )):
            raise ValueError(f"Fortran input/output fixture does not match: {case['name']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--iterations", type=int, default=20000)
    parser.add_argument("--repetitions", type=int, default=7)
    parser.add_argument("--batch", type=int, default=1024)
    parser.add_argument("--family", choices=list(support.FAMILIES.values()))
    parser.add_argument("--fixtures", type=Path, default=support.DEFAULT_FIXTURES)
    parser.add_argument("--module", default="oneloop_native")
    parser.add_argument("--backend", choices=("auto", "native", "symjit", "expression"), default="auto",
                        help="explicit binary64 implementation; auto keeps the module's default")
    parser.add_argument("--build-label", default="unverified", help="caller-supplied label, e.g. release; not inferred or verified")
    parser.add_argument("--output", type=Path, required=True, help="JSON file; engine diagnostics may use stdout")
    parser.add_argument("--tsv", type=Path, help="optional comparator-compatible timing records")
    parser.add_argument("--fortran", type=Path, help="optional matched original benchmark JSON")
    parser.add_argument("--maximum-ratio", type=float, default=1.5)
    args = parser.parse_args()
    if min(args.iterations, args.repetitions, args.batch) < 1:
        parser.error("iterations, repetitions and batch must be positive")
    if not math.isfinite(args.maximum_ratio) or args.maximum_ratio <= 0:
        parser.error("maximum-ratio must be finite and positive")
    if args.output.resolve() in (args.fixtures.resolve(), args.fortran.resolve() if args.fortran else None):
        parser.error("output must not overwrite input files")
    if args.tsv and args.tsv.resolve() in (
        args.output.resolve(), args.fixtures.resolve(), args.fortran.resolve() if args.fortran else None,
    ):
        parser.error("TSV output must be distinct from JSON and input files")
    fixtures = support.read_fixtures(args.fixtures, [args.family] if args.family else None)
    reference = json.loads(args.fortran.read_text()) if args.fortran else None
    if reference is not None:
        matching_reference(reference, fixtures)
        comparison.baseline_data(reference, args.family)
        if reference["iterations"] != args.iterations or reference["repetitions"] != args.repetitions:
            parser.error("Fortran iterations and repetitions must match this survey")

    started = time.perf_counter_ns()
    module = importlib.import_module(args.module)
    startup_ns = time.perf_counter_ns() - started
    if not module.is_initialized():
        raise ValueError("module import did not eagerly initialize every scalar family")
    resolved_backend = (getattr(module, "DEFAULT_BACKEND", "unreported")
                        if args.backend == "auto" else args.backend)
    metadata = {
        "schema": "oneloop_python_benchmark_v1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "machine": support.machine_metadata(), "python_implementation": platform.python_implementation(),
        "module": args.module, "symbolica_revision": module.SYMBOLICA_REVISION,
        "build_label_unverified": args.build_label, "numeric_domain": "Python complex / binary64 components",
        "route": "Evaluator.evaluate_batch (reusable selected-backend workspace)",
        "backend_requested": args.backend, "backend_resolved": resolved_backend,
        "batch": args.batch, "iterations": args.iterations, "repetitions": args.repetitions,
        "warmup": WARMUP, "eager_import_wall_ns_excluded": startup_ns,
        "main_thread": threading.current_thread() is threading.main_thread(),
        "gc_enabled": gc.isenabled(), "fixtures_path": str(args.fixtures.resolve()),
        "fixtures_sha256": support.sha256(args.fixtures),
        "selected_fixture_sha256": hashlib.sha256(support.fixture_text(fixtures).encode()).hexdigest(),
        "survey_source_sha256": support.sha256(__file__),
        "timing_contract": {
            "all_three_laurent_outputs": True, "python_boundary_conversion_included": True,
            "python_result_allocation_included": True, "loop_and_result_retention_included": True,
            "input_list_preparation_excluded": True, "startup_and_backend_clone_excluded": True,
            "output_validation_and_checksums_excluded": True, "final_output_destruction_excluded": True,
            "explicit_mu_squared_last": True,
            "warmup": "256 scalar calls per workload, matching the Rust/Fortran surveys",
            "hetero_order": "cycle family fixtures in file order; restart at first row each repetition",
            "comparison_caveat": "Rust/Fortran include native finite checksums in timing; Python checksums are outside timing",
        },
    }
    module_path = Path(module.__file__).resolve() if getattr(module, "__file__", None) else None
    metadata["extension_path"] = str(module_path) if module_path else None
    metadata["extension_sha256"] = support.sha256(module_path) if module_path and module_path.is_file() else None
    all_samples, summaries, actual = [], [], {}
    for family in support.FAMILIES.values():
        cases = [case for case in fixtures if case["family"] == family]
        if not cases:
            continue
        # Preserve compatibility with older adapters for auto, while explicit
        # requests must be honored (never silently substituted on TypeError).
        evaluator = (module.Evaluator(family) if args.backend == "auto"
                     else module.Evaluator(family, backend=args.backend))
        for case in cases:
            actual[case["id"]] = checked_output(case, evaluator.evaluate(arguments(case)))
        for mode, name, rows in [*(('SAME', case["name"], [case]) for case in cases), ('HETERO', 'all', cases)]:
            samples, summary = workload(evaluator, rows, calls=args.iterations,
                                        repetitions=args.repetitions, batch=args.batch, mode=mode, name=name)
            all_samples.extend(samples)
            summaries.append(summary)
            print(f"{mode} {family} {name}: batch={args.batch}, "
                  f"median={summary['wall_ns_per_call_median']:.3f} ns/call", file=sys.stderr, flush=True)
    result = {
        **metadata, "correctness": "all fixture, untimed layout and timed-output Laurent coefficients passed",
        "fixtures": [{**case, "actual": actual[case["id"]]} for case in fixtures],
        "samples": all_samples, "summaries": summaries,
    }
    tsv = as_tsv(all_samples)
    passed = True
    if reference is not None:
        ratios = comparison.compare(reference, tsv, batch=args.batch, family=args.family, maximum_ratio=args.maximum_ratio)
        result["fortran_comparison"] = {
            "reference_sha256": support.sha256(args.fortran), "maximum_ratio": args.maximum_ratio,
            "metric": "Python-boundary median ns/call / original Fortran median ns/call",
            "workloads": [{"mode": key[0], "family": key[1], "name": key[2], "batch": batch,
                           "fortran_ns": original, "python_ns": native, "ratio": ratio, "passed": ok}
                          for key, batch, original, native, ratio, ok in ratios],
        }
        passed = all(row[-1] for row in ratios)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if args.tsv:
        args.tsv.write_text(tsv)
    return 0 if passed else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, KeyError, TypeError, OSError, ImportError) as error:
        print(f"python performance survey: {error}", file=sys.stderr)
        raise SystemExit(1) from error
