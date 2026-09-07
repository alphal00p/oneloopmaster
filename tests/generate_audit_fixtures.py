#!/usr/bin/env python3
"""Exploratory scalar-domain scan, separate from the passing acceptance table.

Usage: python3 tests/generate_audit_fixtures.py /path/to/compiled/oracle
Stdout is a table for audit_extra_fortran_fixtures; stderr lists any original
Fortran rejections/nonfinite outputs. Such points are not counted as parity
successes. No Symbolica or numerical Rust code generates expected values.
"""
import math
import random
import re
import subprocess
import sys


def cases():
    yield "degenerate_triangle", 3, 1, [-4, -1, -1], [0, 0, 1]
    yield "vacuum_box_equal", 4, 1, [0] * 6, [1] * 4
    yield "vacuum_box_unequal", 4, 1, [0] * 6, [1, 2, 3, 4]
    for p in [-3, -1e-10, 0, 1e-10, 3, 4 - 2**-20, 4 + 2**-20, 7, 12]:
        for m in ([1, 1], [1, 2], [.7-.03j, 1.4-.08j], [0, 2]):
            for kind in [2, -2]:
                yield "bubble_threshold", kind, 1, [p], m
    yield "bubble_exact_threshold", 2, 1, [4], [1, 1]
    for masses in ([1, 4], [4, 1]):
        yield "derivative_pseudothreshold", -2, 1, [1], masses
    for p in ([12, 1, 2], [1, 12, 2], [1, 2, 12], [-12, 1, 2], [1, -12, 2], [1, 2, -12]):
        yield "massless_triangle_positive_kallen", 3, 1, p, [0]*3
    for mass in [-2, -2-.1j, -.1j, 2, 2-.1j]:
        yield "tadpole_mass_domain", 1, 1, [], [mass]
    for mu2 in [1e-12, 1e-8, 1e-2, 1, 3.7, 1e4, 1e8, 1e12]:
        for kind, p, m in [
            (1, [], [2-.1j]),
            (2, [7], [1, 2]),
            (-2, [7], [1, 2]),
            (3, [0, 1, -1], [0]*3),
            (3, [0, 1, 2], [0, 0, 1]),
            (3, [1, 2, 3], [.5-.01j, .7-.02j, 1.4-.05j]),
            (4, [0, 0, 0, 0, 12, -5], [0]*4),
            (4, [0, 2, -3, 3, 5, -6], [0, 0, 2, 3]),
            (4, [-1, -2, -3, -4, 12, -5], [.3-.01j, .5-.02j, .7-.03j, 1.4-.08j]),
        ]:
            yield "independent_mu", kind, mu2, p, m
    rng = random.Random(0xA11CA11)
    for legs in [3, 4]:
        for mask in range(1 << legs):
            for width in [0, .001, .5]:
                for physical in [False, True]:
                    m = [complex(.2 + 2*rng.random(), -width*(.1+rng.random()))
                         if mask & (1 << i) else 0 for i in range(legs)]
                    p = [-.3 - 4*rng.random() for _ in range(3 if legs == 3 else 6)]
                    if physical:
                        if legs == 3:
                            p[1:] = [2 + 8*rng.random(), 2 + 8*rng.random()]
                        else:
                            p[4:] = [8 + 8*rng.random(), -5 - 5*rng.random()]
                    yield "mass_masks_widths", legs, 1, p, m
    for p, m in [
        ([0, 0, -3], [0]*3), ([0, 1, -1], [0]*3),
        ([0, 1, 2], [0, 0, 1]), ([2, 12, 3], [0, 2, 3]),
        ([-1, 3, 7], [.3-.01j, .7-.04j, 1.4-.06j]),
    ]:
        for shift in range(3):
            yield "triangle_cycles", 3, 1, p[shift:]+p[:shift], m[shift:]+m[:shift]
    for first in [-1-.1j, -.5j]:
        for p in [-3, 3]:
            for kind in [2, -2]:
                yield "negative_or_zero_real_mass", kind, 1, [p], [first, 2-.03j]
        for p in ([-1, -2, -12], [-1, 3, 7]):
            yield "negative_or_zero_real_mass", 3, 1, p, [first, .7-.04j, 1.4-.06j]
        yield "negative_or_zero_real_mass", 4, 1, [-1, -2, -3, -4, 12, -5], [first, .5-.02j, .7-.03j, 1.4-.08j]
    p, m = [-1, -2, -12], [.5-.001j, .7-.8j, 1.4-.002j]
    for shift in range(3):
        yield "triangle_large_middle_width", 3, 1, p[shift:]+p[:shift], m[shift:]+m[:shift]
    for s in [-5, 8]:
        for p2, p3, t in [(-2, -3, -6), (7, -3, -6), (-2, 9, 12), (7, 9, 12)]:
            yield "box_ir16_continuation", 4, 1, [2, p2, p3, 4, s, t], [0, 2, 3-.3j, 4]
        for p1, p4, t in [(-1, -4, -6), (7, -4, -6), (-1, 8, 12), (7, 8, 12)]:
            yield "box_ir15_continuation", 4, 1, [p1, 2, 3, p4, s, t], [0, 2, 0, 3]


def main():
    accepted = rejected = 0
    print("# exploratory original-Fortran scalar parity scan; named categories precede rows")
    for name, kind, mu2, p, masses in cases():
        request = [f"{kind} {math.sqrt(mu2):.17e}"]
        request += [f"({complex(z).real:.17e},{complex(z).imag:.17e})" for z in p+masses]
        request.append("stop")
        result = subprocess.run([sys.argv[1]], input="\n".join(request)+"\n",
                                capture_output=True, text=True, check=True)
        pairs = re.findall(r"^\s*olo:\s*(\S+)\s+(\S+)", result.stdout, re.M)
        values = [float(x.replace("D", "E")) for pair in pairs for x in pair]
        if "ERROR" in result.stdout or len(values) != 6 or not all(map(math.isfinite, values)):
            rejected += 1
            print(f"ORACLE_REJECTED {name} kind={kind} mu2={mu2} p={p} masses={masses}\n{result.stdout}", file=sys.stderr)
            continue
        print(f"# {accepted+1} {name}")
        inputs = [mu2] + p + [v for m in masses for v in (complex(m).real, complex(m).imag)]
        print(kind, *(f"{x:.17e}" for x in inputs+values))
        accepted += 1
    print(f"Fortran scan: {accepted} accepted, {rejected} rejected/nonfinite", file=sys.stderr)


if __name__ == "__main__":
    main()
