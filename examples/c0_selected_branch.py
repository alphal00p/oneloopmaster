"""Inspect the reported massive C0 via repo APIs; --full prints its full body.

Run in the shared Symbolica/OneLOop Python host:
    python examples/c0_selected_branch.py [--full]
"""
import sys
import time
from decimal import Decimal

from symbolica import N, Replacement, S
from symbolica.community import oneloop

mass2, mass2B, mu2 = S("mass2", "mass2B", "mu2")
C0 = S("oneloopmaster::C0")
rules = [Replacement(mass2, N(2)), Replacement(mass2B, N(1)), Replacement(mu2, N(1))]
master = C0(0, -mass2, mass2, mass2, mass2, mass2B, mu2)
started = time.perf_counter()
all_branches = oneloop.get_expression(master, max_nodes=100_000_000)
print("Complete expression characters:", len(str(all_branches)))
assert "if(" in str(all_branches)
selected = oneloop.select_branch(all_branches, rules)
print("Expansion and selection (seconds):", time.perf_counter() - started)
assert "if(" not in str(selected) and "__olo_" not in str(selected)
assert mass2 in selected[0].get_all_symbols(False)
assert mass2B in selected[0].get_all_symbols(False)
assert selected[1:] == (N(0), N(0))
print("Finite expression characters:", len(str(selected[0])))
if "--full" in sys.argv:
    print("Repo expressions (finite, simple pole, double pole):", selected)
values = {mass2: 2, mass2B: 1, mu2: 1}
evaluated = selected[0].evaluate(values, decimal_digit_precision=60)
actual = evaluated.to_decimal_tuple()
print("Repo expression evaluated at the probe:", evaluated.to_decimal_tuple(60))
print("Pole coefficients:", selected[1:])

expected = (Decimal("-0.313413985805891556410851729326115326931925200907803748216227"),Decimal(0))
native = oneloop.C0(0, -2, 2, 2, 2, 1, mu_squared=1, prec=60)[0]
for got, want, reference in zip(actual, expected, [native.real, native.imag]):
    assert abs(got - want) < Decimal("1e-45")
    assert abs(got - reference) < Decimal("1e-45")
print("PASS: complete repo expression, late selection, known value, and native backend agree.")
