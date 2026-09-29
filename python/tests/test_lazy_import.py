"""Fresh-process coverage: importing the adapter must not prepare SymJIT."""

import os
import subprocess
import sys
import textwrap
import unittest

from test_api import MODULE


class LazyImport(unittest.TestCase):
    def test_import_and_native_calls_do_not_initialize_symjit(self):
        # An invalid SymJIT configuration must affect only a requested JIT
        # backend, including when this adapter is part of the community host.
        program = textwrap.dedent("""
            import importlib, math, os
            from decimal import Decimal
            module = importlib.import_module(os.environ['ONELOOP_PYTHON_MODULE'])
            assert not module.is_initialized()
            if module.EXPRESSION_INTEROP:
                from symbolica import Expression, S
                mass = S('lazy_exports::mass_squared')
                assert callable(module.A0)
                assert isinstance(module.A0(mass + 1, 1), Expression)
                finite = module.A0(0, mass + 1, 1).evaluate({mass: 1})
                assert abs(finite - 2 * (1 - math.log(2))) < 1e-14
                assert not module.is_initialized()
            result = module.a0(2, mu_squared=1)
            assert abs(result[0] - 2 * (1 - math.log(2))) < 1e-14
            assert result[1:] == (2+0j, 0j)
            precise = module.a0(Decimal(2), mu_squared=Decimal(1), prec=50)
            assert all(v.real.is_finite() and v.imag.is_finite() for v in precise)
            try:
                module.a0(2, mu_squared=1, backend='symjit')
            except RuntimeError as error:
                assert 'unset SYMJIT_TOML' in str(error), str(error)
            else:
                raise AssertionError('SymJIT configuration guard was bypassed')
            del os.environ['SYMJIT_TOML']
            compiled = module.a0(2, mu_squared=1, backend='symjit')
            assert all(abs(a-b) < 1e-14 for a,b in zip(compiled, result))
            assert not module.is_initialized()
        """)
        env = dict(os.environ, ONELOOP_PYTHON_MODULE=MODULE,
                   SYMJIT_TOML='/nonexistent/oneloop-lazy-import.toml')
        result = subprocess.run([sys.executable, '-c', program], env=env,
                                capture_output=True, text=True, timeout=120)
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, output)
        self.assertNotIn('Created infinity', output)
        self.assertNotIn('Created indeterminate', output)
