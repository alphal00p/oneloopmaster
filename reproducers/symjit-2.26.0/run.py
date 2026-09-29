#!/usr/bin/env python3
"""Run each reproducer in a separate process and record failures and controls."""
import json
import os
from pathlib import Path
import platform
import signal
import subprocess

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
# Keep ambient configuration and large core dumps out of these reproducers.
os.environ.pop("SYMJIT_TOML", None)
if (ROOT / "symjit.toml").exists():
    raise SystemExit("Remove symjit.toml from this directory before reproducing.")
try:
    import resource
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
except ImportError:
    pass

subprocess.run(["cargo", "build", "--locked", "--bins"], check=True)
metadata = json.loads(subprocess.check_output(
    ["cargo", "metadata", "--locked", "--no-deps", "--format-version=1"], text=True
))
bin_dir = Path(metadata["target_directory"]) / "debug"
lines = ["SymJIT 2.26.0 reproductions", platform.platform()]
lines.append(subprocess.check_output(["rustc", "--version"], text=True).strip())
cpuinfo = Path("/proc/cpuinfo")
if cpuinfo.exists():
    fields = dict(line.split(":", 1) for line in cpuinfo.read_text().split("\n\n", 1)[0].splitlines() if ":" in line)
    fields = {key.strip(): value.strip() for key, value in fields.items()}
    lines.append("CPU: " + fields.get("model name", "unknown"))
    lines.append("Relevant CPU flags: " + " ".join(
        flag for flag in fields.get("flags", "").split()
        if flag in {"sse2", "avx", "avx2", "fma", "avx512f"}
    ))

cases = [
    ("complex_sqrt", [], "reproducer"),
    ("complex_sqrt", ["generic"], "reproducer"),
    ("complex_if", [], "reproducer"),
    ("complex_if", ["join-only"], "reproducer"),
    ("complex_rounding", [], "reproducer"),
    ("complex_rounding", ["packed"], "control"),
    ("direct_power", [], "reproducer"),
    ("direct_power", ["indirect"], "control"),
    ("direct_negative_powers", [], "reproducer"),
    ("direct_negative_powers", ["indirect"], "control"),
    ("zero_input_batch", [], "reproducer"),
    ("zero_input_batch", ["scalar"], "control"),
]
failures = []
for name, args, kind in cases:
    executable = bin_dir / (name + (".exe" if os.name == "nt" else ""))
    heading = f"\n### {kind}: {name} {' '.join(args)}".rstrip()
    result = subprocess.run(
        [str(executable), *args], stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, timeout=30,
    )
    status = (signal.Signals(-result.returncode).name if result.returncode < 0
              else str(result.returncode))
    entry = f"{heading}\nexit: {status}\n{result.stdout}"
    print(entry, flush=True)
    lines.append(entry)
    if result.returncode != 0:
        failures.append((name, args, status))

(ROOT / "observed.txt").write_text("\n".join(lines) + "\n")
print("Wrote observed.txt.")
if failures:
    raise SystemExit(f"Failed regression cases: {failures}")
print("All twelve regression/control cases passed.")
