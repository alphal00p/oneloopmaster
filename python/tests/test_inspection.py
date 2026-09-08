"""Inspection API checks; use the same single Symbolica worker as test_api."""
import unittest

from test_api import on_symbolica_thread, tearDownModule  # noqa: F401


class ExpressionInspection(unittest.TestCase):
    def test_complete_compact_complex_one_mass_triangle(self):
        def check(module):
            if module.EXPRESSION_INTEROP:
                from symbolica import S
                invariant = S("inspection_python_c0_s", is_real=True)
                mass = S("inspection_python_c0_m")
                scale = S("inspection_python_c0_mu", is_positive=True)
                arguments = [0, 0, invariant, 0, mass, 0, scale]
                options = {}
            else:
                arguments = ["0", "0", "inspection_python_c0_s", "0", "inspection_python_c0_m", "0", "inspection_python_c0_mu"]
                options = {"real": [arguments[2]], "positive": [arguments[6]]}
            result = module.get_expression("C0", arguments, **options)
            for coefficient in result:
                text = str(coefficient)
                self.assertNotIn("__olo_", text)
                self.assertNotIn("C0(", text)
                self.assertLess(len(text), 40_000)
            self.assertIn("if(", str(result[0]))
        on_symbolica_thread(check)

    def test_complete_bodies_and_coefficient_selection(self):
        def check(module):
            if module.EXPRESSION_INTEROP:
                from symbolica import S
                mass = S("inspection_python_mass", is_positive=True)
                scale = S("inspection_python_scale", is_positive=True)
                arguments = [mass, scale]
                options = {}
            else:
                arguments = ["inspection_python_mass", "inspection_python_scale"]
                options = {"positive": arguments}
            result = module.get_expression("A0", arguments, **options)
            self.assertIsInstance(result, tuple)
            self.assertEqual(len(result), 3)
            self.assertEqual(module.get_expression("A0", arguments, coefficient=-1, **options), result[1])
            self.assertNotIn("__olo_", str(result[0]))
            self.assertNotIn("A0(", str(result[0]))
            self.assertNotIn("conj(", str(result[0]))
            if module.EXPRESSION_INTEROP:
                self.assertEqual(result[1], mass)
                self.assertFalse(isinstance(result[0], str))
            else:
                self.assertIsInstance(result[0], str)
                self.assertIn("positive", result[1])
                self.assertIn("oneloop_input", result[1])
        on_symbolica_thread(check)

    def test_all_families_and_explicit_errors(self):
        def check(module):
            for family, momenta, masses, expected in [
                ("A0", 0, 1, "1"), ("B0", 1, 2, "0"),
                ("dB0", 1, 2, "1/6"), ("C0", 3, 3, "-1/2"),
                ("D0", 6, 4, "1/6"),
            ]:
                arguments = [0] * momenta + [1] * (masses + 1)
                if not module.EXPRESSION_INTEROP:
                    arguments = list(map(str, arguments))
                result = module.get_expression(family, arguments)
                self.assertEqual(str(result[0]), expected)
                self.assertNotIn("__olo_", str(result))
            arguments = [1, 1] if module.EXPRESSION_INTEROP else ["1", "1"]
            for options in [{"coefficient": 2}, {"max_nodes": 1}, {"max_depth": 1}, {"max_nodes": 0}]:
                with self.assertRaises(ValueError):
                    module.get_expression("A0", arguments, **options)
            with self.assertRaises(ValueError):
                module.get_expression("D0", arguments)
            if not module.EXPRESSION_INTEROP:
                module.get_expression("A0", ["inspection_python_untyped", "1"])
                with self.assertRaisesRegex(ValueError, "already exists"):
                    module.get_expression("A0", ["inspection_python_untyped", "1"], real=["inspection_python_untyped"])
                native = module.get_expression("A0", ["inspection_attributes::{positive}::m", "1"], coefficient=-1)
                self.assertIn("positive", native)
        on_symbolica_thread(check)
