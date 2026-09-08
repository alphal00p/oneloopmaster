"""Decimal/native-Float boundaries; no Symbolica work occurs during discovery.

Run serially with the other adapter tests. This module reuses their explicit
large-stack worker, never a hidden implementation worker. Decimal references
use Python's independent decimal transcendental implementation.
"""

from decimal import Decimal, getcontext, localcontext
import unittest

from test_api import on_symbolica_thread, tearDownModule


def decimal_multiply(first, second):
    a, b = first
    c, d = second
    return a * c - b * d, a * d + b * c


def decimal_reciprocal(value):
    real, imag = value
    norm = real * real + imag * imag
    return real / norm, -imag / norm


def decimal_complex_log(value, digits):
    """Independent log for the positive-real, narrow-width controls below."""
    real, imag = value
    x = imag / real
    if abs(x) > Decimal("0.125"):
        raise AssertionError("atan series controls require a small argument")
    power = x
    angle = x
    tolerance = Decimal(10) ** (-digits - 40)
    for index in range(1, 4 * digits + 100):
        power *= -x * x
        term = power / (2 * index + 1)
        angle += term
        if abs(term) < tolerance:
            return (real * real + imag * imag).ln() / 2, angle
    raise AssertionError("independent Decimal atan did not converge")


def vacuum_coefficients(family, mass, scale, digits):
    """Elementary scalar vacuum identities, independent of Symbolica."""
    zero = Decimal(0), Decimal(0)
    one = Decimal(1), Decimal(0)
    with localcontext() as context:
        context.prec = digits + 60
        if family in ("A0", "B0"):
            logarithm = decimal_complex_log(mass, digits)
            logarithm = logarithm[0] - scale.ln(), logarithm[1]
            if family == "A0":
                finite = decimal_multiply(mass, (1 - logarithm[0], -logarithm[1]))
                return finite, mass, zero
            return (-logarithm[0], -logarithm[1]), one, zero
        inverse = decimal_reciprocal(mass)
        if family == "dB0":
            return (inverse[0] / 6, inverse[1] / 6), zero, zero
        if family == "C0":
            return (-inverse[0] / 2, -inverse[1] / 2), zero, zero
        square = decimal_multiply(inverse, inverse)
        return (square[0] / 6, square[1] / 6), zero, zero


class PrecisionTests(unittest.TestCase):
    def assert_decimal_complex(self, value, module):
        self.assertIsInstance(value, module.DecimalComplex)
        self.assertIsInstance(value.real, Decimal)
        self.assertIsInstance(value.imag, Decimal)

    def assert_relative_decimal(self, actual, expected, digits):
        with localcontext() as context:
            context.prec = digits + 40
            tolerance = Decimal(10) ** (3 - digits) * max(abs(expected), Decimal(1))
            self.assertLessEqual(abs(actual - expected), tolerance)

    def test_decimal_context_does_not_truncate_inputs_or_outputs(self):
        def check(module):
            for digits in (32, 1000):
                with self.subTest(prec=digits), localcontext() as ambient:
                    ambient.prec = 28
                    mass = Decimal("1." + "1234567890" * 110)
                    scale = Decimal("3." + "9876543210" * 110)
                    result = module.A0(mass, mu_squared=scale, prec=digits)
                    self.assertEqual(getcontext().prec, 28)
                    for value in result:
                        self.assert_decimal_complex(value, module)
                        self.assertEqual(value.imag, 0)
                    with localcontext() as reference_context:
                        reference_context.prec = digits + 40
                        expected = mass * (1 - (mass / scale).ln())
                    self.assert_relative_decimal(result[0].real, expected, digits)
                    self.assert_relative_decimal(result[1].real, mass, digits)
                    self.assertEqual(result[2].real, 0)
                    self.assertGreaterEqual(len(result[0].real.as_tuple().digits), digits - 3)
        on_symbolica_thread(check)

    def test_default_machine_path_and_decimal_selection(self):
        def check(module):
            self.assertTrue(all(isinstance(z, complex) for z in module.A0(2)))
            self.assertTrue(all(isinstance(z, complex) for z in module.A0(2, prec=16)))
            self.assertTrue(all(isinstance(z, complex) for z in module.A0(2**53)))
            for arguments in ((Decimal("2"),), (2**53 + 1,)):
                values = module.A0(*arguments)
                self.assertTrue(all(isinstance(z, module.DecimalComplex) for z in values))
            self.assertEqual(module.A0(2**53 + 1)[1].real, Decimal(2**53 + 1))
            self.assertEqual(module.A0(2**100, prec=50)[1].real, Decimal(2**100))
            # Decimal(int) avoids CPython's default int->str digit limit.
            self.assert_relative_decimal(module.A0(2**20000, prec=32)[1].real, Decimal(2**20000), 32)
            # Raising precision preserves the supplied float's actual binary
            # value; it must not replace it with the shorter decimal str(float).
            result = module.A0(0.1, prec=80)
            self.assertEqual(result[1].real, Decimal.from_float(0.1))
            decimal_result = module.A0(Decimal("0.1"), prec=80)
            self.assertNotEqual(result[1].real, decimal_result[1].real)
            # Decimal scale alone selects AP, and cannot underflow through f64.
            scaled = module.A0(1, mu_squared=Decimal("1e-1000"), prec=32)
            with localcontext() as context:
                context.prec = 80
                expected = 1 - 1000 * Decimal(10).ln()
            self.assert_relative_decimal(scaled[0].real, expected, 32)
        on_symbolica_thread(check)

    def test_decimal_complex_preserves_tiny_width_and_large_precision(self):
        def check(module):
            with localcontext() as context:
                context.prec = 28
                mass = module.DecimalComplex(Decimal("2.125"), Decimal("-1e-1000"))
                self.assertEqual(mass.real, Decimal("2.125"))
                self.assertEqual(mass.imag, Decimal("-1e-1000"))
                values = module.A0(mass, prec=1050)
                for value in values:
                    self.assert_decimal_complex(value, module)
                self.assertEqual(values[1].real, mass.real)
                with localcontext() as reference_context:
                    reference_context.prec = 1100
                    self.assert_relative_decimal(values[1].imag / mass.imag, Decimal(1), 1000)
                    # F'(m)=-log(m) for A0(m,mu²=1). The omitted relative
                    # correction is O(width²), below 1e-2000 for this input.
                    self.assert_relative_decimal(values[0].imag / mass.imag, -mass.real.ln(), 1000)
                self.assertEqual(complex(mass), complex(2.125, -0.0))
                # Construction from strings is exact as well, even at context28.
                constructed = module.DecimalComplex("1.123456789012345678901234567890123456789", "-0")
                self.assertEqual(str(constructed.real), "1.123456789012345678901234567890123456789")
                self.assertTrue(constructed.imag.is_signed())
                self.assertIn("DecimalComplex", repr(constructed))
        on_symbolica_thread(check)

    def test_all_five_families_at_32_and_1000_digits(self):
        def check(module):
            for digits in (32, 1000):
                for family in ("A0", "B0", "dB0", "C0", "D0"):
                    momenta, masses = {"A0": (0, 1), "B0": (1, 2), "dB0": (1, 2),
                                       "C0": (3, 3), "D0": (6, 4)}[family]
                    evaluator = module.Evaluator(family, prec=digits)
                    rows, expected = [], []
                    for index in range(5):
                        real = Decimal("2.125") + index
                        imag = Decimal("-0.125") if index % 2 else Decimal(0)
                        scale = (Decimal("1e-30"), Decimal("4.75"), Decimal("1e30"))[index % 3]
                        mass = module.DecimalComplex(real, imag) if imag else real
                        row = [0] * momenta + [mass] * masses + [scale]
                        rows.append(row)
                        expected.append(vacuum_coefficients(family, (real, imag), scale, digits))
                    results = [getattr(module, family)(*row[:-1], mu_squared=row[-1], prec=digits)
                               for row in rows]
                    results += [evaluator.evaluate(row) for row in rows]
                    # Five rows retain the partial-four-lane-tail API shape,
                    # although this fixed arbitrary-precision path is scalar.
                    results += evaluator.evaluate_batch(rows)
                    self.assertEqual(len(results), 15)
                    for index, actual in enumerate(results):
                        for coefficient, (value, target) in enumerate(zip(actual, expected[index % 5])):
                            with self.subTest(prec=digits, family=family, row=index % 5,
                                              route=index // 5, coefficient=coefficient):
                                self.assert_decimal_complex(value, module)
                                self.assert_relative_decimal(value.real, target[0], digits)
                                self.assert_relative_decimal(value.imag, target[1], digits)
        on_symbolica_thread(check)

    def test_reusable_precision_overrides_and_mixed_batch_tails(self):
        def check(module):
            evaluator = module.Evaluator("A0", prec=32)
            self.assertEqual(evaluator.prec, 32)
            rows = [[Decimal("1.25") + index, Decimal("2.75")] for index in range(5)]
            singles = [evaluator.evaluate(row) for row in rows]
            for length in (1, 3, 4, 5):
                values = evaluator.evaluate_batch(rows[:length])
                for actual, expected in zip(values, singles):
                    for a, e in zip(actual, expected):
                        self.assertEqual((a.real, a.imag), (e.real, e.imag))
            self.assertEqual(evaluator.evaluate_batch([]), [])
            bulk = evaluator.evaluate_batch([rows[index % len(rows)] for index in range(1024)])
            self.assertEqual(len(bulk), 1024)
            for index, actual in enumerate(bulk):
                for value, expected in zip(actual, singles[index % len(rows)]):
                    self.assertEqual((value.real, value.imag), (expected.real, expected.imag))
            higher = evaluator.evaluate(rows[0], prec=1000)
            self.assertGreaterEqual(len(higher[0].real.as_tuple().digits), 997)
            self.assertEqual(evaluator.prec, 32)  # Overrides do not mutate default.
            self.assertTrue(all(isinstance(z, complex) for z in evaluator.evaluate([1, 1], prec=16)))
            # One Decimal anywhere promotes the complete mixed batch, not only
            # that row; all numerical work has one requested precision.
            default = module.Evaluator("A0")
            promoted = default.evaluate_batch([[2.0, 1.0], [Decimal("2.1"), 1.0], [3, 1]])
            self.assertTrue(all(isinstance(z, module.DecimalComplex) for row in promoted for z in row))
            evaluator.rebuild()
            rebuilt = evaluator.evaluate(rows[0])
            self.assertEqual(rebuilt[0].real, singles[0][0].real)
        on_symbolica_thread(check)

    def test_decimal_validation_never_uses_float_predicates(self):
        def check(module):
            for digits in (0, -1, 2**40):
                with self.subTest(prec=digits), self.assertRaises(ValueError):
                    module.A0(1, prec=digits)
            for digits in (True, False, 2.5, "32", Decimal(32)):
                with self.subTest(prec=digits), self.assertRaises(TypeError):
                    module.A0(1, prec=digits)
                with self.assertRaises(TypeError):
                    module.Evaluator("A0", prec=digits)
            for value in (Decimal("NaN"), Decimal("sNaN"), Decimal("Infinity"), Decimal("-Infinity")):
                with self.subTest(value=str(value)), self.assertRaises(ValueError):
                    module.A0(value)
            with self.assertRaises(ValueError):
                module.A0(module.DecimalComplex(1, Decimal("1e-10000")))
            with self.assertRaises(ValueError):
                module.A0(module.DecimalComplex(1, Decimal("1e-1000000000000")))
            with self.assertRaises(ValueError):
                module.B0(module.DecimalComplex(-1, Decimal("1e-10000")), 1, 1)
            for scale in (Decimal(0), Decimal("-1e-10000"), module.DecimalComplex(1, "1e-10000")):
                with self.subTest(scale=str(scale)), self.assertRaises(ValueError):
                    module.A0(1, mu_squared=scale)
            for value in (True, False):
                with self.assertRaises(TypeError):
                    module.A0(value)
            evaluator = module.Evaluator("A0")
            with self.assertRaises(TypeError):
                evaluator.evaluate_batch([[Decimal(2), 1], [True, 1]])
            with self.assertRaises(ValueError):
                evaluator.evaluate_batch([], prec=0)
            with self.assertRaises(ValueError):
                evaluator.evaluate([Decimal(2)])
            with self.assertRaises(ValueError):
                evaluator.evaluate_batch([[Decimal(2), 1], [Decimal(3)]])
        on_symbolica_thread(check)

    def test_undefined_integral_outputs_remain_decimal_nonfinite_values(self):
        def check(module):
            # The fully scaleless derivative is undefined, unlike scaleless B0.
            # Preserve its numerical result rather than failing during MPFR to
            # Decimal text conversion or silently inventing a finite value.
            result = module.dB0(Decimal(0), 0, 0, prec=32)
            self.assertTrue(all(isinstance(z, module.DecimalComplex) for z in result))
            self.assertTrue(any(not z.real.is_finite() or not z.imag.is_finite() for z in result))
        on_symbolica_thread(check)
