import unittest
from pkg_resources import resource_filename as pkg_file
from genomvar import Reference
from genomvar.varset import VariantBase
from genomvar import variant

CHR15RGN = pkg_file('genomvar.test','data/chr15rgn.fna')
# Factory normalizing indels
nvf = variant.VariantFactory(reference=CHR15RGN,normindel=True)
# Factory not normalizing indels
svf = variant.VariantFactory()

class TestVariantsCase(unittest.TestCase):
    def test_indel_equality(self):
        del1 = variant.AmbigDel('chr1',(1,1),(5,3),'CA','')
        del2 = variant.AmbigDel('chr1',(1,2),(5,4),'AC','')
        self.assertTrue(del1.ambig_equal(del2))
        self.assertFalse(del1.edit_equal(del2))

    def test_instantiation_from_edit(self):
        # **Simple insertion**
        #                15
        # TTCACTTAGCATAATG|TCTT
        #                 C
        vb = nvf.from_edit('chr15rgn',15,'G','GC')
        self.assertEqual([vb.start,vb.end],[16,17])
        # Now deletion
        vb = nvf.from_edit(chrom='chr15rgn',start=15,ref='GT',alt='G')
        self.assertEqual([vb.start,vb.end],[16,17])
        vb = nvf.from_edit('chr15rgn',15,'GT','G')
        self.assertEqual([vb.start,vb.end],[16,17])

        #    575
        # ATTTAATA
        #    T-AT   v1
        vb = nvf.from_edit('chr15rgn',575,'TA','T')
        self.assertTrue(isinstance(vb,variant.AmbigIndel))

        #       65
        # TAAGG CTG     AATACTAT
        #       CTC
        vb = svf.from_edit('chr15rgn',start=65,ref='CTG',alt='CTC')
        self.assertEqual([vb.start,vb.end,vb.ref,vb.alt,type(vb).__name__],
                         [67,68,'G','C','SNP'])

        #       65
        # TAAGG CTG     AATACTAT
        #       CCGTCGTG
        vb = svf.from_edit('chr15rgn',start=65,ref='CTG',alt='CCGTCGTG')
        self.assertEqual([vb.start,vb.end,bool(vb.ref),vb.alt,type(vb).__name__],
                         [66,67,False,'CGTCG','Ins'])

        #     2343
        # CTG TTTCCA    ACATACATCATGAGACTTCTG
        #     TTCCATTCCA
        vb = nvf.from_edit('chr15rgn',2343,'TTTCCA','TTCCATTCCA')
        self.assertEqual([vb.start,vb.end,vb.alt,type(vb).__name__],
                         [2344,2346,'CCAT','AmbigIns'])

        #      3300
        #   TATC TTTTTGAC TGG
        #        --------
        vb = nvf.from_edit('chr15rgn',3300,'CTTTTTGAC','C')
        self.assertEqual([vb.start,vb.end,vb.alt,type(vb).__name__],
                         [3300,3310,'','AmbigDel'])

        #   0
        #   TTCACTTAGCA
        vb = nvf.from_edit('chr15rgn',0,'T','TT')
        self.assertEqual([vb.start,vb.end,vb.alt,vb.vtp],
                         [0,3,'T',variant.AmbigIns])
    def test_instantiation_from_hgvs(self):
        # test SNP
        vb = svf.from_hgvs('chr1:g.15C>A')
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['chr1',14,15,'C','A'])
        vb = svf.from_hgvs('NG_012232.1:g.19_21del')
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['NG_012232.1',18,21,None,None])
        vb = svf.from_hgvs('NG_012232.1:g.19del')
        self.assertTrue(vb.is_instance(variant.Del))
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['NG_012232.1',18,19,None,None])
        vb = svf.from_hgvs('NC_000023.10:g.10_11insCCT')
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['NC_000023.10',10,11,None,'CCT'])
        vb = svf.from_hgvs('NC_000023.11:g.10delinsGA')
        self.assertTrue(vb.is_instance(variant.Mixed))
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['NC_000023.11',9,10,None,'GA'])
        vb = svf.from_hgvs('LRG_199:g.145_147delinsTGG')
        self.assertTrue(vb.is_instance(variant.MNP))
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['LRG_199',144,147,None,'TGG'])

        # Test Ambig
        #      3300
        #   TATC TTTTTGAC TGG
        #        --------
        ve = nvf.from_edit('chr15rgn',3299,'TCTTTTTGA','T')
        vb = nvf.from_hgvs('chr15rgn:g.3302_3309del')
        self.assertEqual([vb.start,vb.end,vb.alt,vb.vtp],
                         [3300,3310,'',variant.AmbigDel])
        self.assertFalse(vb.edit_equal(ve))
        self.assertTrue(vb.ambig_equal(ve))
