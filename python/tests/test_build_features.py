"""Numeric API coverage for wheels with and without generated evaluators."""
from decimal import Decimal
import unittest

from test_api import family_selector, on_symbolica_thread, tearDownModule  # noqa: F401


class BuildFeatures(unittest.TestCase):
    def test_auto_backend_and_explicit_native_availability(self):
        def check(module):
            def components(values):
                return [(z.real, z.imag) for z in values]

            self.assertIn(module.DEFAULT_BACKEND, ("native", "expression"))
            family = family_selector(module, "A0")
            for digits, point in [(16, [2, 1]), (80, [Decimal(2), Decimal(1)])]:
                evaluator = module.Evaluator(family, prec=digits)
                self.assertEqual(evaluator.backend, "auto")
                result = evaluator.evaluate(point)
                self.assertEqual(result[1].real, 2)
                self.assertEqual(result[1].imag, 0)
                self.assertEqual(components(evaluator.evaluate_batch([point])[0]), components(result))
                if module.DEFAULT_BACKEND == "native":
                    native = module.Evaluator(family, prec=digits, backend="native")
                    self.assertEqual(components(native.evaluate(point)), components(result))
                else:
                    with self.assertRaisesRegex(ValueError, "generated-evaluators"):
                        module.Evaluator(family, prec=digits, backend="native")
                    with self.assertRaisesRegex(ValueError, "generated-evaluators"):
                        module.a0(*point, prec=digits, backend="native")
            # The compact build still supports an explicitly requested JIT.
            self.assertEqual(module.a0(2, 1, backend="symjit")[1], 2)

        on_symbolica_thread(check)
