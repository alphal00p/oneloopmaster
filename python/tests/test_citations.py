"""Global citation metadata remains lazy until OneLoopMaster is used."""
import unittest
from test_api import on_symbolica_thread


class Citations(unittest.TestCase):
    def test_authors_and_requested_papers(self):
        def check(module):
            if not module.EXPRESSION_INTEROP:
                self.skipTest("citations use the shared Symbolica host")
            from symbolica import S, get_citations
            self.assertFalse(hasattr(module, "get_citations"))
            initialized = module.is_initialized()
            mass = S("citation_mass")
            module.A0(mass, 1)
            citations = get_citations()
            self.assertEqual({c.id for c in citations}, {
                "doi:10.5281/zenodo.17054381",
                "https://github.com/alphal00p/oneloopmaster",
                "arXiv:1007.4716", "arXiv:0903.4665",
            })
            self.assertTrue(all(c.reference and c.reasons for c in citations))
            self.assertTrue(all(c.to_bibtex().startswith("@") for c in citations))
            self.assertEqual(module.is_initialized(), initialized)
        on_symbolica_thread(check)
