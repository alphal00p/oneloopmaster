"""Expression-native inspection and condition-only branch selection."""
from decimal import Decimal
import unittest

from test_api import on_symbolica_thread, tearDownModule  # noqa: F401


class ExpressionInspection(unittest.TestCase):
    def in_host(self, check):
        def run(module):
            if not module.EXPRESSION_INTEROP:
                self.skipTest("symbolic APIs require the shared-kernel Python host")
            check(module)
        on_symbolica_thread(run)

    def test_readme_parametric_equal_mass_bubble_expression(self):
        def check(module):
            from symbolica import E, Expression, N, Replacement, S
            psq = S("inspection_readme_psq", is_real=True)
            mass = S("inspection_readme_m2", is_positive=True)
            master = S("oneloopmaster::B0")(psq, mass, mass, 1)
            result = module.get_expression(master)
            self.assertEqual(module.get_expression(master, coefficient=0), result[0])
            self.assertTrue(all(isinstance(c, Expression) for c in result))
            self.assertIn("if(", str(result[0]))
            self.assertNotIn("__olo_", str(result))
            self.assertNotIn("B0(", str(result))

            for sample, numeric in [(N("3.23"), 3.23), (N(3), 3),
                                    (E("1/2"), 0.5), (N(-2), -2),
                                    (N(5), 5), (N(0), 0), (N(4), 4)]:
                with self.subTest(sample=str(sample)):
                    rules = [Replacement(mass, N(1)), Replacement(psq, sample)]
                    selected = module.select_branch(result, rules)
                    self.assertIsInstance(selected, tuple)
                    self.assertNotIn("if(", str(selected))
                    self.assertEqual(selected[1:], (N(1), N(0)))
                    symbols = selected[0].get_all_symbols(False)
                    self.assertIn(mass, symbols)
                    if numeric != 0:
                        self.assertIn(psq, symbols)
                    value = selected[0].replace_multiple(rules).evaluate({})
                    reference = module.B0(numeric, 1, 1, mu_squared=1)[0]
                    self.assertLess(abs(value - reference), 2e-12)
                    self.assertEqual(module.select_branch(result[0], rules), selected[0])
                    self.assertEqual(module.select_branch(list(result), rules), list(selected))

            # Native tags on a compound squared mass are used by Symbolica;
            # neither this API nor the caller needs an assumption-name list.
            compound = S("oneloop::B0")(psq, mass + E("12/3"), mass + 4, 1)
            selected = module.select_branch(module.get_expression(compound), [
                Replacement(psq, N("3.23")), Replacement(mass, N(1)),
            ])
            self.assertNotIn("if(", str(selected))
            self.assertIn(mass, selected[0].get_all_symbols(False))
            namespaced = E(
                "oneloop::B0(inspection_ns::{real}::psq,"
                "inspection_ns::{real,positive}::m2+12/3,"
                "inspection_ns::{real,positive}::m2+4,1)"
            )
            native_tags = module.select_branch(module.get_expression(namespaced), [
                Replacement(S("inspection_ns::psq"), N("3.23")),
                Replacement(S("inspection_ns::m2"), N(1)),
            ])
            self.assertNotIn("if(", str(native_tags))
            self.assertIn(S("inspection_ns::m2"), native_tags[0].get_all_symbols(False))
        self.in_host(check)

    def test_only_conditions_change_and_unresolved_conditions_stay_parametric(self):
        def check(module):
            from symbolica import E, N, Replacement, S
            x, y, z = S("inspection_scope_x", "inspection_scope_y", "inspection_scope_z")
            native_if = S("if")
            rules = [Replacement(x, N(2)), Replacement(y, N(3))]
            expr = native_if(native_if(x, y, z), native_if(y, x + y, z), x - z)
            self.assertEqual(module.select_branch(expr, rules), x + y)
            unresolved = native_if(x + z, native_if(y, x, z), y)
            self.assertEqual(module.select_branch(unresolved, rules), native_if(x + z, x, y))
            self.assertEqual(module.select_branch(x + y, rules), x + y)
            self.assertEqual(module.select_branch(native_if(E("1/2") * x, x, y), rules), x)
            # Native IF tests zero/nonzero. For real x, x+abs(x) tests x>0;
            # Symbolica's Python comparison objects are matcher conditions,
            # not expression-tree predicates in this pinned development API.
            comparison = native_if(x + S("abs")(x), x, y)
            self.assertEqual(module.select_branch(comparison, [Replacement(x, N("4.23"))]), x)
            tiny = native_if((1 + x) ** E("1/2") - 1, x, y)
            self.assertEqual(module.select_branch(tiny, [Replacement(x, N("1e-1000"))]), x)
        self.in_host(check)

    def test_all_nested_ifs_in_functions_powers_and_deep_branches(self):
        def check(module):
            from symbolica import N, Replacement, S
            x, y, f = S("inspection_deep_x", "inspection_deep_y", "inspection_deep_f")
            native_if = S("if")
            branch = native_if(x, y, 17)
            expr = f(branch) + (branch + 1) ** branch * branch
            rules = [Replacement(x, N(1))]
            self.assertEqual(module.select_branch(expr, rules), f(y) + (y + 1) ** y * y)
            deep = y
            for _ in range(1200):
                deep = native_if(x, deep, 0)
            self.assertEqual(module.select_branch(deep, rules), y)
            wildcard = S("inspection_rule_a_")
            expr = native_if(f(x + 1), x, y)
            self.assertEqual(module.select_branch(expr, [
                Replacement(f(wildcard), N(0)), Replacement(f(wildcard), N(1)),
            ]), y)
            conditional_rules = [
                Replacement(f(wildcard), N(0), cond=wildcard.req_gt(N(0))),
                Replacement(f(wildcard), N(1)),
            ]
            self.assertEqual(module.select_branch(native_if(f(N(2)), x, y), conditional_rules), y)
            self.assertEqual(module.select_branch(native_if(f(N(-2)), x, y), conditional_rules), x)
        self.in_host(check)

    def test_complete_bodies_composites_coefficient_selection_and_errors(self):
        def check(module):
            from symbolica import E, N, Replacement, S
            for name, momenta, masses, expected in [
                ("A0", 0, 1, N(1)), ("B0", 1, 2, N(0)),
                ("dB0", 1, 2, E("1/6")), ("C0", 3, 3, E("-1/2")),
                ("D0", 6, 4, E("1/6")),
            ]:
                master = S("oneloopmaster::" + name)(*([0] * momenta + [1] * (masses + 1)))
                result = module.get_expression(master)
                self.assertEqual(result[0], expected)
                self.assertNotIn("__olo_", str(result))
            mass = S("inspection_compound_mass", is_positive=True)
            scale = S("inspection_compound_scale", is_positive=True)
            master = S("oneloopmaster::A0")(mass + 4, scale ** 2)
            result = module.get_expression(master)
            self.assertEqual(result[1], mass + 4)
            self.assertEqual(module.get_expression(master, coefficient=-1), result[1])
            self.assertNotIn("conj(", str(result))
            for options in [{"coefficient": 2}, {"max_nodes": 1}, {"max_depth": 1}, {"max_nodes": 0}]:
                with self.assertRaises(ValueError):
                    module.get_expression(master, **options)
            for invalid in [mass, master + 1, S("other::A0")(mass, 1), S("oneloopmaster::B0")(mass, 1)]:
                with self.assertRaises(ValueError):
                    module.get_expression(invalid)
            for operation in [
                lambda: module.get_expression("B0", ["psq", "m2", "m2", "1"]),
                lambda: module.get_expression(str(master)),
                lambda: module.get_expression(master, real=["inspection_compound_mass"]),
                lambda: module.select_branch("if(x,x,0)", []),
                lambda: module.select_branch(master, [(mass, N(1))]),
                lambda: module.master_coefficients("A0", [mass, 1]),
                lambda: module.compile_native(["x"], [mass]),
                lambda: module.Evaluator("A0"),
                lambda: module.Evaluator(mass + 1),
            ]:
                with self.assertRaises(TypeError):
                    operation()
        self.in_host(check)

    def test_triangle_native_expression_and_numeric_symbolica_inputs(self):
        def check(module):
            from symbolica import E, N, S
            invariant = S("inspection_c0_s", is_real=True)
            mass = S("inspection_c0_m")
            scale = S("inspection_c0_mu", is_positive=True)
            master = S("oneloopmaster::C0")(0, 0, invariant, 0, mass, 0, scale)
            result = module.get_expression(master)
            self.assertIn("if(", str(result[0]))
            self.assertNotIn("__olo_", str(result))
            self.assertLess(len(str(result[0])), 40_000)
            # A variable Expression is accepted where Symbol is required.
            evaluator = module.Evaluator(E("oneloopmaster::A0"), prec=64)
            actual = evaluator.evaluate([E("1/3"), N(1)])[1]
            self.assertIsInstance(actual.real, Decimal)
            self.assertLess(abs(actual.real - Decimal("0." + "3" * 64)), Decimal("1e-63"))
            with self.assertRaises(ValueError):
                evaluator.evaluate([mass, N(1)])
        self.in_host(check)
