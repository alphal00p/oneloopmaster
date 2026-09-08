"""Explicit numeric backends; finite controls, not universal parity evidence."""
from decimal import Decimal, localcontext
import unittest

from test_api import FAMILIES, on_symbolica_thread, read_fixture_groups, tearDownModule  # noqa: F401
from test_precision import vacuum_coefficients


class BackendSelection(unittest.TestCase):
    def assert_machine_close(self, actual, expected):
        for a, b in zip(actual, expected):
            self.assertIsInstance(a, complex)
            self.assertLessEqual(abs(a - b), 5e-10 * max(abs(b), 1.0))

    def assert_decimal_close(self, actual, expected, digits):
        with localcontext() as context:
            context.prec = digits + 40
            for a, b in zip(actual, expected):
                for component, reference in zip((a.real, a.imag), b):
                    self.assertIsInstance(component, Decimal)
                    self.assertLessEqual(abs(component - reference), Decimal(10) ** (4 - digits) * max(abs(reference), Decimal(1)))

    def test_machine_native_and_symjit_all_benchmark_points(self):
        def check(module):
            for kind, rows in read_fixture_groups().items():
                family = FAMILIES[kind]
                native = module.Evaluator(family, backend="native")
                jit = module.Evaluator(family, backend="symjit")
                expression = module.Evaluator(family, backend="expression")
                inputs = [arguments for _, arguments, _, _ in rows]
                native_batch = native.evaluate_batch(inputs)
                jit_batch = jit.evaluate_batch(inputs)
                for index, arguments in enumerate(inputs):
                    with self.subTest(family=family, row=index):
                        self.assert_machine_close(native_batch[index], jit_batch[index])
                        self.assert_machine_close(native.evaluate(arguments), jit_batch[index])
                        self.assert_machine_close(getattr(module, family)(*arguments, backend="native"), jit_batch[index])
                self.assert_machine_close(expression.evaluate(inputs[0]), jit_batch[0])
                self.assertEqual(native.backend, "native")
        on_symbolica_thread(check)

    def test_arbitrary_native_and_expression_all_families(self):
        def check(module):
            for digits in (32, 1000):
                for family, momenta, masses in [("A0", 0, 1), ("B0", 1, 2), ("dB0", 1, 2), ("C0", 3, 3), ("D0", 6, 4)]:
                    for mass_pair, scale in [
                        ((Decimal("2.125"), Decimal(0)), Decimal("1e-30")),
                        ((Decimal("1.23456789"), Decimal("-0.0123456789")), Decimal("4.75")),
                    ]:
                        arguments = [Decimal(0)] * momenta + [module.DecimalComplex(*mass_pair)] * masses + [scale]
                        reference = vacuum_coefficients(family, mass_pair, scale, digits)
                        for backend in ("native", "expression"):
                            with self.subTest(family=family, prec=digits, backend=backend):
                                evaluator = module.Evaluator(family, prec=digits, backend=backend)
                                self.assert_decimal_close(evaluator.evaluate(arguments), reference, digits)
                                self.assert_decimal_close(evaluator.evaluate_batch([arguments])[0], reference, digits)
                                self.assert_decimal_close(getattr(module, family)(*arguments, prec=digits, backend=backend), reference, digits)
        on_symbolica_thread(check)

    def test_genuinely_complex_triangle_backend_agreement(self):
        def check(module):
            arguments = [Decimal(-2), Decimal(-3), Decimal(-5)] + [
                module.DecimalComplex("1", "-0.1"), module.DecimalComplex("2", "-0.2"),
                module.DecimalComplex("3", "-0.3"),
            ] + [Decimal("4.75")]
            for digits in (32, 1000):
                reference = module.C0(*arguments, prec=digits, backend="expression")
                expected = [(z.real, z.imag) for z in reference]
                native = module.Evaluator("C0", prec=digits, backend="native")
                self.assert_decimal_close(native.evaluate(arguments), expected, digits)
                self.assert_decimal_close(native.evaluate_batch([arguments])[0], expected, digits)
        on_symbolica_thread(check)

    def test_batched_1024_native_machine_and_decimal(self):
        def check(module):
            machine_rows = [[1.125 + index / 1024, 4.75] for index in range(1024)]
            native = module.Evaluator("A0", backend="native")
            expected = module.Evaluator("A0", backend="symjit").evaluate_batch(machine_rows)
            actual = native.evaluate_batch(machine_rows)
            self.assertEqual(len(actual), 1024)
            for a, b in zip(actual, expected):
                self.assert_machine_close(a, b)
            decimal_rows = [[Decimal("1.125") + Decimal(index) / 1024, Decimal("4.75")] for index in range(1024)]
            native = module.Evaluator("A0", prec=32, backend="native")
            actual = native.evaluate_batch(decimal_rows)
            expected = module.Evaluator("A0", prec=32, backend="expression").evaluate_batch(decimal_rows)
            for a, b in zip(actual, expected):
                self.assert_decimal_close(a, [(z.real, z.imag) for z in b], 32)
            native.rebuild()
            self.assert_decimal_close(native.evaluate(decimal_rows[0]), [(z.real, z.imag) for z in expected[0]], 32)
        on_symbolica_thread(check)

    def test_selection_overrides_and_arbitrary_symjit_rejection(self):
        def check(module):
            self.assertIn(module.DEFAULT_BACKEND, ("native", "symjit", "expression"))
            evaluator = module.Evaluator("A0", backend="native")
            for backend in ("native", "symjit", "expression", "symbolica", "auto"):
                expected = module.A0(2.125, 4.75, backend=backend)
                self.assert_machine_close(evaluator.evaluate([2.125, 4.75], backend=backend), expected)
            self.assertEqual(evaluator.backend, "native")
            self.assertEqual(module.Evaluator("A0", backend="symbolica").backend, "expression")
            self.assertEqual(module.Evaluator("A0").backend, "auto")
            for value in (Decimal("2"), 2**53 + 1, module.DecimalComplex("2", "-1e-1000")):
                with self.assertRaisesRegex(ValueError, "binary64"):
                    module.A0(value, backend="symjit")
                with self.assertRaisesRegex(ValueError, "binary64"):
                    evaluator.evaluate([value, 1], backend="symjit")
            with self.assertRaisesRegex(ValueError, "binary64"):
                module.A0(2, mu_squared=Decimal(1), backend="symjit")
            with self.assertRaisesRegex(ValueError, "binary64"):
                module.Evaluator("A0", prec=32, backend="symjit")
            with self.assertRaisesRegex(ValueError, "binary64"):
                evaluator.evaluate_batch([], prec=1000, backend="symjit")
            with self.assertRaisesRegex(ValueError, "binary64"):
                evaluator.evaluate_batch([[2.0, 1.0], [Decimal(2), 1]], backend="symjit")
            for backend in ("native", "expression", "auto"):
                values = evaluator.evaluate([Decimal("2"), 1], backend=backend)
                self.assertTrue(all(isinstance(value, module.DecimalComplex) for value in values))
            for bad in ("", "cuda", "jit"):
                with self.assertRaises(ValueError):
                    module.Evaluator("A0", backend=bad)
                with self.assertRaises(ValueError):
                    evaluator.evaluate_batch([], backend=bad)
            with self.assertRaises(TypeError):
                module.A0(2, backend=True)
        on_symbolica_thread(check)
