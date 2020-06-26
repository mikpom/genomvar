import unittest
from pkg_resources import resource_filename as pkg_file
from genomvar.varset import IndexedVariantFileSet,VariantFileSet

class TestPerformanceCase(unittest.TestCase):
    def test_vcf_vs_itself(self):
        vs1 = IndexedVariantFileSet(pkg_file('genomvar.test',
                                      'data/test_vcf_1kg_2.vcf.gz'),
                             parse_samples=True)
        vs2 = IndexedVariantFileSet(pkg_file('genomvar.test',
                                      'data/test_vcf_1kg_2.vcf.gz'),
                             parse_samples=True)
        comm = []
        for vrt in vs1.comm_vrt(vs2).region(rgn='7:152134922-152436005'):
            comm.append(vrt)
        nof_vrt = len(list(vs1.iter_vrt()))
        self.assertEqual(len(comm),nof_vrt)

        

if __name__ == "__main__":
    unittest.main()        
