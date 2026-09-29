"""Run each allocator/thread-lifetime check in a fresh process (Linux/POSIX)."""
import argparse
from collections import Counter
import resource
import signal
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("driver")
parser.add_argument("library")
parser.add_argument("--main", action="store_true")
parser.add_argument("--expect-crash", action="store_true")
parser.add_argument("--repeat", type=int, default=20)
args = parser.parse_args()
if args.repeat < 1:
    parser.error("--repeat must be positive")
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
command = [args.driver, args.library] + (["main"] if args.main else [])
expected = -signal.SIGSEGV if args.expect_crash else 0
counts = Counter()
for _ in range(args.repeat):
    result = subprocess.run(command, capture_output=True, text=True, timeout=10)
    counts[result.returncode] += 1
    if result.returncode != expected:
        raise SystemExit(
            f"Expected {expected}, got {result.returncode}\n"
            f"{result.stdout}{result.stderr}"
        )
print(f"expected={expected}, return_codes={dict(counts)}")
