"""Fresh-process standalone smoke helper; never start a Symbolica worker here.

Invoked by test_api.py, or directly as: python cold_start.py D0 evaluator.
Import must eagerly initialize all scalar families before any master call.
Only ready-workspace paths (Native or prebuilt SymJIT) are exercised. Source rebuilding
still belongs on the adequately stacked calling thread described in the README.
"""

import argparse
import importlib
import json
import math
import threading

from test_api import FAMILIES, read_fixture_groups


def bounded_main_stack():
    if threading.current_thread() is not threading.main_thread():
        raise AssertionError("cold loading must run on the process main thread")
    try:
        import resource
    except ImportError:
        return None  # Non-POSIX platforms retain their ordinary main-thread stack.
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    limit = 8 * 1024 * 1024
    if soft != resource.RLIM_INFINITY:
        limit = min(limit, soft)
    if hard != resource.RLIM_INFINITY:
        limit = min(limit, hard)
    # A test-only upper bound: never increase the process's inherited stack.
    resource.setrlimit(resource.RLIMIT_STACK, (limit, hard))
    return limit


def check_output(output, case):
    line, _, expected, normalization = case
    if len(output) != 3 or not all(isinstance(value, complex) for value in output):
        raise AssertionError("expected three Python complex coefficients")
    for coefficient, (actual, reference) in enumerate(zip(output, expected)):
        error = abs(actual - reference) * normalization
        tolerance = 2e-10 + 2e-8 * abs(reference) * normalization
        if not math.isfinite(error) or error > tolerance:
            raise AssertionError(f"fixture line {line}, coefficient {coefficient}: {actual} != {reference}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("family", choices=tuple(FAMILIES.values()))
    parser.add_argument("entry", choices=("direct", "evaluator", "guard-env", "guard-file"))
    parser.add_argument("--rebuild-only", action="store_true")
    parser.add_argument("--backend", choices=("auto", "native", "symjit"), default="auto")
    args = parser.parse_args()
    stack_limit = bounded_main_stack()
    kind = next(kind for kind, name in FAMILIES.items() if name == args.family)
    groups = read_fixture_groups()
    record = {
        "family": args.family, "entry": args.entry,
        "main_thread": threading.current_thread() is threading.main_thread(),
        "stack_soft_limit_bytes": stack_limit, "checked_rows": 0, "batch_rows": 0,
        "backend_requested": args.backend,
    }
    if args.entry.startswith("guard-"):
        expected = (
            "unset SYMJIT_TOML when building or loading portable OneLOop evaluators"
            if args.entry == "guard-env" else
            "remove the working-directory symjit.toml override when building or loading portable OneLOop evaluators"
        )
        try:
            importlib.import_module("oneloop_native")
        except RuntimeError as error:
            if str(error) != expected:
                raise AssertionError(f"wrong eager-startup environment error: {error}") from error
        else:
            raise AssertionError("eager module import accepted a portable environment override")
        record["guarded_imports"] = 1
        print("ONELOOP_COLD_JSON " + json.dumps(record))
        return
    elif args.rebuild_only:
        raise AssertionError("cold loading must not use a source-rebuild fallback")

    # No master function or Evaluator has been called before import completes.
    module = importlib.import_module("oneloop_native")
    if module.EXPRESSION_INTEROP:
        raise AssertionError("this smoke test requires the standalone extension")
    record["backend_resolved"] = module.DEFAULT_BACKEND if args.backend == "auto" else args.backend
    # This query is read-only in the core: it cannot make this assertion pass
    # by lazily creating a master symbol, function map, or scalar cache.
    if not module.is_initialized():
        raise AssertionError("import returned before all symbols and five caches were ready")
    record["import_completed_before_master_calls"] = True
    record["initialized_before_master_calls"] = True
    # Exercise the requested first entry route, then every family without
    # another expensive import. Both routes also check heterogeneous batches.
    ordered_kinds = [kind] + [other for other in FAMILIES if other != kind]
    checked = {}
    for current_kind in ordered_kinds:
        family = FAMILIES[current_kind]
        cases = groups[current_kind]
        if args.entry == "direct":
            function = getattr(module, family)
            for case in cases:
                inputs = case[1]
                check_output(function(*inputs[:-1], mu_squared=inputs[-1], backend=args.backend), case)
            evaluator = module.Evaluator(family, backend=args.backend)
        else:
            evaluator = module.Evaluator(family, backend=args.backend)
            for case in cases:
                check_output(evaluator.evaluate(case[1]), case)
        ordered = [cases[index % len(cases)] for index in range(9)]
        outputs = evaluator.evaluate_batch([case[1] for case in ordered])
        if len(outputs) != len(ordered):
            raise AssertionError("incorrect batch row count")
        for output, case in zip(outputs, ordered):
            check_output(output, case)
        record["checked_rows"] += len(cases)
        record["batch_rows"] += len(ordered)
        checked[family] = {"single_rows": len(cases), "batch_rows": len(ordered)}
    record["ready_families"] = list(FAMILIES.values())
    record["checks_by_family"] = checked
    print("ONELOOP_COLD_JSON " + json.dumps(record))


if __name__ == "__main__":
    main()
