#!/usr/bin/env python3
"""Regenerate native-expression fixtures using the original OneLOop library.

Usage: python3 tests/generate_fixtures.py /path/to/compiled/oracle
Compile tests/support/oracle.f90 against the unchanged original Fortran library.
The table is written to stdout. No Python packages or wrapper library are needed.
"""
import math
import random
import re
import subprocess
import sys


def cases():
    rng = random.Random(0x5EED)
    for legs in (3, 4):
        for sample in range(96):
            mask = sample % (1 << legs)
            scale = 10 ** (-4 + 8 * rng.random())
            p = [-scale * (0.2 + 3 * rng.random()) for _ in range(3 if legs == 3 else 6)]
            if legs == 4:
                p[4:] = [-12 * scale, -15 * scale]
            m = [complex(scale * (0.1 + 2 * rng.random()), -scale * (0.001 + .03 * rng.random()) if legs == 3 else 0) if mask & (1 << i) else 0j for i in range(legs)]
            yield legs, scale, p, m
        for sample in range(16):
            mask = sample % (1 << legs)
            m = [complex(.3+2*rng.random(),-.01-.08*rng.random()) if mask & (1<<j) else 0j for j in range(legs)]
            p = ([1,2,3] if sample%2==0 else [-1,2,3]) if legs==3 else [-1,-2,-3,-4,12,-5]
            yield legs,1,p,m
    yield 3, 1, [1, 2, 3], [.5, .7, 1.4]
    # Exact finite boundaries and both orientations of opposite masses.
    for p in ([-1,-2,3,-4,-5,-6], [-1,-2,-3,3,-5,-6]):
        yield 4, 1, p, [0,2,0,3]
    yield 4, 1, [-1,-2,-3,-4,-5,-6], [1.4,0,.7,0]
    # Widths and physical scattering kinematics are deliberately included.
    for m in ([.3-.01j,.5-.02j,.7-.03j,1.4-.08j], [0,.5-.02j,.7-.03j,1.4-.08j]):
        for p in ([-1,-2,-3,-4,-12,-15], [0,0,0,0,12,-5]):
            yield 4, 1, p, m
    # Independent infrared references, including both mixed-sign orderings.
    for p in ([0, 0, -3], [0, -2, -3], [0, 1, -1], [0, -1, 1]):
        yield 3, 1, p, [0, 0, 0]
    for p in ([0, 2, 2], [0, -3, 2], [0, -3, -4]):
        yield 3, 1, p, [0, 0, 2]
    yield 3, 1, [2, -3, 3], [0, 2, 3]
    for p in ([0, 1, -2], [0, 1, 2], [0, 2, 3]):
        yield 3, 1, p, [0, 0, 1]
    yield 3, 1, [2, 12, 3], [0, 2, 3]
    for p in (
        [0, 0, 0, 0, -5, -6],
        [0, 0, 0, -4, -5, -6],
        [0, 0, -3, -4, -5, -6],
        [0, -2, 0, -4, -5, -6],
        [0, -2, -3, -4, -5, -6],
        [0, 0, 0, 0, -2, 1],
        [0, 0, 0, 0, 1, -2],
    ):
        yield 4, 1, p, [0, 0, 0, 0]
    for p, m in (
        ([0, 0, 2, 2, -5, -6], [0, 0, 0, 2]),
        ([0, 0, -3, 2, -5, -6], [0, 0, 0, 2]),
        ([0, 0, -3, -4, -5, -6], [0, 0, 0, 2]),
        ([-1, 0, 2, -4, -5, -6], [0, 0, 0, 2]),
        ([-1, 0, -3, -4, -5, -6], [0, 0, 0, 2]),
        ([0, 2, -3, 3, -5, -6], [0, 0, 2, 3]),
        ([0, 2, -3, -4, -5, -6], [0, 0, 2, 3]),
        ([0, -2, -3, -4, -5, -6], [0, 0, 2, 3]),
        ([2, 2, 3, 3, -5, -6], [0, 2, 0, 3]),
        ([-1, 2, 3, 3, -5, -6], [0, 2, 0, 3]),
        ([2, -2, -3, 4, -5, -6], [0, 2, 3, 4]),
        ([0, 0, 2, 2, 5, -6], [0, 0, 0, 2]),
        ([0, 2, -3, 3, 5, -6], [0, 0, 2, 3]),
    ):
        yield 4, 1, p, m
    for mu2, mass in ((1, 0), (1, 2), (4, 2), (1, 100 - 1.4j)):
        yield 1, mu2, [], [mass]
    bubbles = (
        (5, [1, 1]),
        (7, [1, 2]),
        (9, [2, 1]),
        (3, [1, 1]),
        (-3, [1, 1]),
        (0, [1, 1]),
        (0, [2, 5]),
        (3, [0, 0]),
        (-1, [0, 2]),
        (2, [0, 2]),
        (3, [.7 - .03j, 1.4 - .08j]),
        (3, [1.4 - .08j, .7 - .03j]),
        (0, [0, 0]),
    )
    for p, m in bubbles:
        yield 2, 1, [p], m
        # The scaleless derivative is indeterminate in the original Fortran,
        # so it is not part of this finite-value comparison table.
        if p != 0 or any(m):
            yield -2, 1, [p], m  # dB0/dp², with the same input convention.


def main():
    samples = list(cases())
    lines = []
    for n, mu2, p, m in samples:
        lines.append(f"{n} {math.sqrt(mu2):.17e}")
        lines.extend(f"({complex(z).real:.17e},{complex(z).imag:.17e})" for z in p + m)
        lines.append("next")
    lines[-1] = "stop"
    output = subprocess.run([sys.argv[1]], input="\n".join(lines)+"\n", text=True, capture_output=True, check=True).stdout
    if "ERROR" in output:
        raise RuntimeError("Fortran rejected a fixture:\n" + output)
    coefficients = re.findall(r"^\s*olo:\s*(\S+)\s+(\S+)", output, re.M)
    if len(coefficients) != 3 * len(samples):
        raise RuntimeError("Incomplete Fortran output")
    rows = ["# kind(1=A0,2=B0,-2=dB0,3=C0,4=D0) mu_squared momenta... squared_mass_real_imag... finite_real_imag pole1_real_imag pole2_real_imag"]
    for i, (n,mu2,p,m) in enumerate(samples):
        result = [float(x.replace('D','E')) for pair in coefficients[3*i:3*i+3] for x in pair]
        if not all(math.isfinite(x) for x in result):
            raise RuntimeError(f"Nonfinite Fortran fixture {i}")
        values = [mu2] + p + [x for z in m for x in (complex(z).real, complex(z).imag)] + result
        rows.append(" ".join([str(n)] + [f"{x:.17e}" for x in values]))
    print("\n".join(rows))


if __name__ == '__main__':
    main()
