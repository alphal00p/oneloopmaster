#!/usr/bin/env python3
"""Optional warmed original-Fortran scalar throughput survey (stdlib only).

Compile the driver against an EXISTING, unchanged original DP build:
  python3 tests/benchmark_fortran.py /path/to/OneLOop --output /tmp/fortran.json
Or use a precompiled tests/support/bench_oracle.f90 executable:
  python3 tests/benchmark_fortran.py --executable /tmp/oneloop-bench

The shared tests/data/benchmark.txt format is the existing scalar fixture format.
--emit-fixtures PATH exports exactly the selected rows for another backend.
SAME repeats each row; HETERO cycles a family's rows in file order. A call must
compute finite/simple-pole/double-pole outputs. The timed checksum adds only the
finite complex coefficient, while every full output is validated outside timing.

A <=1.5x comparison means native median ns/call / original median ns/call <=1.5
for EACH agreed family/workload, on the same machine, with matched fixtures,
total calls, warmup, output requirements and numeric domain. Report native SIMD
batch size (including tails), release/backend/optimization settings separately.
Do not include JIT/graph construction in this warmed-throughput ratio or infer
cold-start latency, universal-domain performance or numerical correctness from it.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent.parent
FAMILIES = {1: "A0", 2: "B0", -2: "dB0", 3: "C0", 4: "D0"}
DEFAULT_FIXTURES = ROOT / "tests/data/benchmark.txt"
DRIVER = ROOT / "tests/support/bench_oracle.f90"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_fixtures(path, families):
    fixtures = []
    name = source = None
    for line_number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if line.startswith("# name "):
            words = line.split()
            name = words[2]
            source = words[4] if len(words) == 5 and words[3] == "source" else None
        if not line.strip() or line.startswith("#"):
            continue
        values = [float(word) for word in line.split()]
        if not values or not all(map(math.isfinite, values)):
            raise ValueError(f"Nonfinite/empty fixture at {path}:{line_number}")
        kind = int(values[0])
        if kind not in FAMILIES or kind != values[0]:
            raise ValueError(f"Unsupported kind at {path}:{line_number}")
        np = 0 if kind == 1 else 1 if abs(kind) == 2 else 3 if kind == 3 else 6
        nm = abs(kind)
        if len(values) != 2 + np + 2 * nm + 6 or values[1] <= 0:
            raise ValueError(f"Malformed fixture at {path}:{line_number}")
        current_name = name or f"{FAMILIES[kind]}_line_{line_number}"
        current_source = source or f"{path}:{line_number}"
        name = source = None
        if families and FAMILIES[kind] not in families:
            continue
        p = values[2:2 + np]
        mass_values = values[2 + np:2 + np + 2 * nm]
        scale = max(map(abs, p + mass_values), default=1.0) or 1.0
        fixtures.append({
            "id": len(fixtures) + 1,
            "name": current_name,
            "source": current_source,
            "kind": kind,
            "family": FAMILIES[kind],
            "mu_squared": values[1],
            "momenta": p,
            "masses": [mass_values[i:i + 2] for i in range(0, len(mass_values), 2)],
            "expected": [values[i:i + 2] for i in range(len(values) - 6, len(values), 2)],
            "normalization": scale ** {1: -1, 2: 0, -2: 1, 3: 1, 4: 2}[kind],
            "row": line.strip(),
        })
    if not fixtures:
        raise ValueError("No benchmark fixtures selected")
    return fixtures


def fixture_text(fixtures):
    lines = ["# Shared scalar benchmark subset; numeric rows retain original fixture values."]
    for case in fixtures:
        lines.extend([f"# name {case['name']} source {case['source']}", case["row"]])
    return "\n".join(lines) + "\n"


def validate_output(case, actual, label):
    if len(actual) != 3:
        raise ValueError(f"{label}: expected all three Laurent outputs")
    for coefficient, (value, expected) in enumerate(zip(actual, case["expected"])):
        if not all(map(math.isfinite, value)):
            raise ValueError(f"{label}: nonfinite coefficient {coefficient}: {value}")
        error = math.hypot(value[0] - expected[0], value[1] - expected[1]) * case["normalization"]
        tolerance = 2e-10 + 2e-8 * math.hypot(*expected) * case["normalization"]
        if not math.isfinite(error) or error > tolerance:
            raise ValueError(f"{label}: coefficient {coefficient} disagrees: {value} vs {expected}")


def parse_results(stdout, fixtures, repetitions, calls):
    if "ERROR" in stdout:
        raise ValueError("Original Fortran reported an ERROR:\n" + stdout)
    values = {}
    samples = []
    clock = None
    for line in stdout.splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == "VALUE":
            if len(words) != 9:
                raise ValueError("Malformed VALUE record: " + line)
            index, kind = map(int, words[1:3])
            if index in values or not 1 <= index <= len(fixtures):
                raise ValueError("Duplicate/invalid fixture result: " + line)
            case = fixtures[index - 1]
            if kind != case["kind"]:
                raise ValueError("Wrong fixture family: " + line)
            numeric = [float(v.replace("D", "E")) for v in words[3:]]
            output = [numeric[i:i + 2] for i in range(0, 6, 2)]
            validate_output(case, output, case["name"])
            values[index] = output
        elif words[0] == "CLOCK":
            if len(words) != 3 or clock is not None:
                raise ValueError("Malformed/duplicate CLOCK record")
            rate, maximum = map(int, words[1:])
            if rate <= 0:
                raise ValueError("Invalid clock rate")
            clock = {"ticks_per_second": rate, "maximum_tick": maximum, "resolution_ns": 1e9 / rate}
        elif words[0] == "SAMPLE":
            if len(words) != 17:
                raise ValueError("Malformed SAMPLE record: " + line)
            mode = words[1]
            kind, fixture_id, repetition, measured_calls, distinct_rows = map(int, words[2:7])
            numbers = [float(v.replace("D", "E")) for v in words[7:]]
            if mode not in ("SAME", "HETERO") or kind not in FAMILIES:
                raise ValueError("Invalid workload: " + line)
            if measured_calls != calls or not 1 <= repetition <= repetitions:
                raise ValueError("Unexpected benchmark call/repetition count")
            if not all(map(math.isfinite, numbers)) or numbers[0] <= 0 or numbers[1] < 0:
                raise ValueError("Invalid/nonfinite timing; increase --iterations if clock resolution is insufficient")
            indices = ([fixture_id] if mode == "SAME" else
                       [case["id"] for case in fixtures if case["kind"] == kind])
            if mode == "HETERO" and fixture_id != 0:
                raise ValueError("HETERO must use fixture id zero")
            if (not indices or len(indices) != distinct_rows or
                    any(not 1 <= i <= len(fixtures) or fixtures[i - 1]["kind"] != kind for i in indices)):
                raise ValueError("Invalid workload fixture selection")
            last = indices[(calls - 1) % len(indices)]
            last_output = [numbers[i:i + 2] for i in range(4, 10, 2)]
            validate_output(fixtures[last - 1], last_output, f"{mode} last output")
            cycles, tail = divmod(calls, len(indices))
            counts = [cycles + int(offset < tail) for offset in range(len(indices))]
            if any(i not in values for i in indices):
                raise ValueError("Timing precedes correctness output")
            expected_sum = [math.fsum(count * values[i][0][component] for i, count in zip(indices, counts))
                            for component in (0, 1)]
            absolute_sum = math.fsum(count * math.hypot(*values[i][0]) for i, count in zip(indices, counts))
            checksum_error = math.hypot(numbers[2] - expected_sum[0], numbers[3] - expected_sum[1])
            if checksum_error > max(1e-12, 1e-9 * absolute_sum):
                raise ValueError(f"{mode} checksum mismatch: {numbers[2:4]} vs {expected_sum}")
            samples.append({
                "mode": mode, "kind": kind, "family": FAMILIES[kind], "fixture_id": fixture_id,
                "name": fixtures[fixture_id - 1]["name"] if mode == "SAME" else FAMILIES[kind] + "_heterogeneous",
                "repetition": repetition, "calls": calls, "distinct_rows": distinct_rows,
                "heterogeneous_tail_rows": tail if mode == "HETERO" else 0,
                "wall_seconds": numbers[0], "cpu_seconds": numbers[1],
                "wall_ns_per_call": numbers[0] * 1e9 / calls,
                "cpu_ns_per_call": numbers[1] * 1e9 / calls,
                "finite_checksum": numbers[2:4], "last_output": last_output,
            })
    if len(values) != len(fixtures) or clock is None:
        raise ValueError("Missing correctness results or clock metadata")
    expected_workloads = len(fixtures) + len({case["kind"] for case in fixtures})
    if len(samples) != expected_workloads * repetitions:
        raise ValueError("Incomplete benchmark sample inventory")
    groups = {}
    for sample in samples:
        key = (sample["mode"], sample["kind"], sample["fixture_id"])
        groups.setdefault(key, []).append(sample)
    if len(groups) != expected_workloads:
        raise ValueError("Missing benchmark workload")
    summaries = []
    for group in groups.values():
        if sorted(sample["repetition"] for sample in group) != list(range(1, repetitions + 1)):
            raise ValueError("Duplicate/missing timing repetition")
        first = group[0]
        wall = [sample["wall_ns_per_call"] for sample in group]
        cpu = [sample["cpu_ns_per_call"] for sample in group]
        summaries.append({
            **{key: first[key] for key in ("mode", "kind", "family", "fixture_id", "name", "calls", "distinct_rows")},
            "wall_ns_per_call_median": statistics.median(wall),
            "wall_ns_per_call_min": min(wall), "wall_ns_per_call_max": max(wall),
            "cpu_ns_per_call_median": statistics.median(cpu),
        })
    return clock, values, samples, summaries


def machine_metadata():
    cpu = None
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.partition(":")[2].strip()
                break
    return {
        "platform": platform.platform(), "machine": platform.machine(), "processor": cpu or platform.processor(),
        "logical_cpus": os.cpu_count(), "python": platform.python_version(),
        "affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("build_dir", nargs="?", type=Path, help="existing original OneLOop module/archive directory")
    parser.add_argument("--executable", type=Path, help="precompiled bench_oracle.f90; skip compilation")
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--fflags", default="-O2", help="driver-only compiler flags; original archive is never rebuilt")
    parser.add_argument("--fixtures", type=Path, default=DEFAULT_FIXTURES)
    parser.add_argument("--family", choices=list(FAMILIES.values()), action="append", help="repeat to select families")
    parser.add_argument("--iterations", type=int, default=20000, help="scalar calls per timed repetition, not repetitions of whole table")
    parser.add_argument("--warmup", type=int, default=256, help="untimed scalar calls per workload")
    parser.add_argument("--repetitions", type=int, default=7)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--emit-fixtures", type=Path)
    parser.add_argument("--output", type=Path, help="JSON output path; stdout if omitted")
    args = parser.parse_args()
    if args.iterations < 1 or args.warmup < 0 or args.repetitions < 1 or args.timeout <= 0:
        parser.error("iterations/repetitions/timeout must be positive; warmup must be nonnegative")
    fixtures = read_fixtures(args.fixtures, args.family)
    selected_text = fixture_text(fixtures)
    if args.emit_fixtures:
        if args.emit_fixtures.resolve() == args.fixtures.resolve():
            parser.error("--emit-fixtures must not overwrite the input fixture table")
        args.emit_fixtures.write_text(selected_text)
    if not args.build_dir and not args.executable:
        if args.emit_fixtures:
            return
        parser.error("provide an existing original build directory or --executable")
    build_dir = args.build_dir.resolve() if args.build_dir else None
    library = build_dir / "libavh_olo.a" if build_dir else None
    metadata = {
        "schema": "oneloop_fortran_benchmark_v1", "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "machine": machine_metadata(), "numeric_domain": "original OneLOop DP (binary64 complex)",
        "fixtures_path": str(args.fixtures.resolve()), "fixtures_sha256": sha256(args.fixtures),
        "selected_fixture_sha256": hashlib.sha256(selected_text.encode()).hexdigest(),
        "driver_source_sha256": sha256(DRIVER),
        "library_path": str(library) if library else None,
        "library_sha256": sha256(library) if library else None,
        "library_build_flags": "unverified: prebuilt archive reused unchanged",
        "iterations": args.iterations, "warmup": args.warmup, "repetitions": args.repetitions,
        "timing_contract": {
            "all_three_laurent_outputs": True, "input_decode_and_sqrt_mu_outside_timing": True,
            "explicit_per_call_unsquared_scale": True, "timed_checksum": "sum of finite complex coefficient",
            "compile_and_initialization_excluded": True, "simd_batch_size": 1,
            "same_mode": "repeat one row", "hetero_mode": "cycle same-family rows in fixture order; restart at first row each repetition",
            "loop_and_accumulator_overhead_included": True,
        },
    }
    with tempfile.TemporaryDirectory(prefix="oneloop-fortran-benchmark-") as temporary:
        if args.executable:
            executable = args.executable.resolve()
            metadata["compiler"] = None
            metadata["driver_compile_command"] = None
        else:
            compiler = shutil.which(args.compiler)
            if compiler is None:
                raise ValueError(f"Compiler not found: {args.compiler}")
            executable = Path(temporary) / "oneloop-bench"
            command = [compiler, *shlex.split(args.fflags), "-I", str(build_dir), str(DRIVER), str(library), "-o", str(executable)]
            metadata["compiler"] = subprocess.run([compiler, "--version"], text=True, capture_output=True, check=True).stdout.splitlines()[0]
            metadata["driver_compile_command"] = command
            compiled = subprocess.run(command, text=True, capture_output=True, timeout=args.timeout)
            if compiled.returncode:
                raise ValueError("Driver compilation failed:\n" + compiled.stdout + compiled.stderr)
        metadata["executable_sha256"] = sha256(executable)
        request = ["OLO_BENCH_V1", f"{len(fixtures)} {args.warmup} {args.iterations} {args.repetitions}"]
        request.extend(case["row"] for case in fixtures)
        started = time.perf_counter()
        measured = subprocess.run([str(executable)], input="\n".join(request) + "\n", text=True,
                                  capture_output=True, timeout=args.timeout)
        metadata["subprocess_wall_seconds_not_benchmark_time"] = time.perf_counter() - started
        if measured.returncode or "ERROR" in measured.stderr:
            raise ValueError("Original benchmark failed:\n" + measured.stdout + measured.stderr)
        clock, values, samples, summaries = parse_results(measured.stdout, fixtures, args.repetitions, args.iterations)
    result = {**metadata, "correctness": "all fixture outputs, batch-final outputs and finite checksums passed",
              "clock": clock, "fixtures": [{**case, "actual": values[case["id"]]} for case in fixtures],
              "summaries": summaries, "samples": samples}
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        sys.stdout.write(encoded)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        print(f"benchmark_fortran: {error}", file=sys.stderr)
        raise SystemExit(1) from error
