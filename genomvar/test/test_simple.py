from unittest import TestCase
from genomvar import Reference
from pkg_resources import resource_filename as pkg_file

class test_Reference_class(TestCase):
    def test_cached_ref_retrieval(self):
        REF = Reference(pkg_file('genomvar.test','data/chr24.fna'))

        self.assertEqual(REF.get('chr24',16,19),'TCT')
        self.assertEqual(REF.get('chr24',17,20),'CTT')
        self.assertEqual(REF.get('chr24',860,861),'C')
        self.assertTrue(REF.get('chr24',860,860+1213).endswith('ACA'))
        self.assertTrue(REF.get('chr24',860,860+100).endswith('GCATTTT'))
        with self.assertRaises(IndexError):
            s = REF.get('chr24',-1,0)

        REF.close()
