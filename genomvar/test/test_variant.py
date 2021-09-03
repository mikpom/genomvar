import copy
from pkg_resources import resource_filename as pkg_file
from genomvar import Reference
from genomvar.varset import VariantBase
from genomvar import variant
from genomvar.vcf import VCFReader, VCFWriter
from genomvar.vcf_utils import VCF_FIELDS, VCFRow
from genomvar.variant import GenomVariant
from genomvar.test import MyTestCase

# Factory normalizing indels

class TestVariantsCase(MyTestCase):
    writer = VCFWriter()
        
    def test_indel_equality(self):
        #       1
        # R     TCACAG
        # del1  T--CAG
        # del2  TCA--G
        adel1 = variant.AmbigDel('chrom',(1,1),(5,3),'CA','')
        adel2 = variant.AmbigDel('chrom',(1,2),(5,4),'AC','')
        self.assertTrue(adel1.ambig_equal(adel2))
        self.assertFalse(adel1.edit_equal(adel2))

        del1 = variant.Del('chrom',1,3)
        self.assertTrue(del1.edit_equal(adel1))
        self.assertTrue(del1.ambig_equal(adel1))
        self.assertTrue(adel1.edit_equal(del1))
        self.assertTrue(adel1.ambig_equal(del1))

    def test_instantiation_from_edit(self):
        # **Simple insertion**
        #                15
        # TTCACTTAGCATAATG|TCTT
        #                 C
        vb = self.nvf.from_edit('chr24',15,'G','GC')
        self.assertEqual([vb.start,vb.end],[16,17])
        # Now deletion
        vb = self.nvf.from_edit(chrom='chr24',start=15,ref='GT',alt='G')
        self.assertEqual([vb.start,vb.end],[16,17])
        vb = self.nvf.from_edit('chr24',15,'GT','G')
        self.assertEqual([vb.start,vb.end],[16,17])

        #    575
        # ATTTAATA
        #    T-AT   v1
        vb = self.nvf.from_edit('chr24',575,'TA','T')
        self.assertTrue(isinstance(vb,variant.AmbigIndel))

        #       65
        # TAAGG CTG     AATACTAT
        #       CTC
        vb = self.svf.from_edit('chr24',start=65,ref='CTG',alt='CTC')
        self.assertEqual([vb.start,vb.end,vb.ref,vb.alt,type(vb).__name__],
                         [67,68,'G','C','SNP'])

        #       65
        # TAAGG CTG     AATACTAT
        #       CCGTCGTG
        vb = self.svf.from_edit('chr24',start=65,ref='CTG',alt='CCGTCGTG')
        self.assertEqual([vb.start,vb.end,bool(vb.ref),vb.alt,type(vb).__name__],
                         [66,67,False,'CGTCG','Ins'])

        #     2343
        # CTG TTTCCA    ACATACATCATGAGACTTCTG
        #     TTCCATTCCA
        vb = self.nvf.from_edit('chr24',2343,'TTTCCA','TTCCATTCCA')
        self.assertEqual([vb.start,vb.end,vb.alt,type(vb).__name__],
                         [2344,2346,'CCAT','AmbigIns'])

        #      3300
        #   TATC TTTTTGAC TGG
        #        --------
        vb = self.nvf.from_edit('chr24',3300,'CTTTTTGAC','C')
        self.assertEqual([vb.start,vb.end,vb.alt,type(vb).__name__],
                         [3300,3310,'','AmbigDel'])

        #   0
        #   TTCACTTAGCA
        vb = self.nvf.from_edit('chr24',0,'T','TT')
        self.assertEqual([vb.start,vb.end,vb.alt,vb.vtp],
                         [0,3,'T',variant.AmbigIns])
    def test_instantiation_from_hgvs(self):
        # test SNP
        vb = self.svf.from_hgvs('chr1:g.15C>A')
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['chr1',14,15,'C','A'])
        vb = self.svf.from_hgvs('NG_012232.1:g.19_21del')
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['NG_012232.1',18,21,None,None])
        vb = self.svf.from_hgvs('NG_012232.1:g.19del')
        self.assertTrue(vb.is_variant_instance(variant.Del))
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['NG_012232.1',18,19,None,None])
        vb = self.svf.from_hgvs('NC_000023.10:g.10_11insCCT')
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['NC_000023.10',10,11,None,'CCT'])
        vb = self.svf.from_hgvs('NC_000023.11:g.10delinsGA')
        self.assertTrue(vb.is_variant_instance(variant.Mixed))
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['NC_000023.11',9,10,None,'GA'])
        vb = self.svf.from_hgvs('LRG_199:g.145_147delinsTGG')
        self.assertTrue(vb.is_variant_instance(variant.MNP))
        self.assertEqual([vb.chrom,vb.start,vb.end,vb.ref,vb.alt],
                         ['LRG_199',144,147,None,'TGG'])

        # Test Ambig
        #      3300
        #   TATC TTTTTGAC TGG
        #        --------
        ve = self.nvf.from_edit('chr24',3299,'TCTTTTTGA','T')
        vb = self.nvf.from_hgvs('chr24:g.3302_3309del')
        self.assertEqual([vb.start,vb.end,vb.alt,vb.vtp],
                         [3300,3310,'',variant.AmbigDel])
        self.assertFalse(vb.edit_equal(ve))
        self.assertTrue(vb.ambig_equal(ve))

    def test_haplotype_edit_equality(self):
        factory = variant.VariantFactory()
        v1 = factory.from_edit('chr24',2093,'TGG','CCC')
        v2 = factory.from_edit('chr24',2098,'TT','GG')
        v3 = factory.from_edit('chr24',2098,'TT','CC')
        h1 = variant.Haplotype.from_variants([v1,v2])
        h1_ = variant.Haplotype.from_variants([v1,v2])
        h2 = variant.Haplotype.from_variants([v1,v3])
        self.assertTrue(h1.edit_equal(h1_))
        self.assertFalse(h1.edit_equal(h2))

    def test_get_vcf_row_instantiated_variant(self):
        factory = variant.VariantFactory()
        v1 = factory.from_edit('chr24',2093,'TGG','CCC')

        row = self.writer.get_row(v1)
        self.assertEqual(row.REF, 'TGG')
        self.assertEqual(row.POS, 2094)
        self.assertEqual(str(row), 'chr24\t2094\t.\tTGG\tCCC\t.\t.\t.')

        gv = GenomVariant(v1, attrib={'id':'vrtid', 'filter':'LOWQUAL',
                                      'qual':100})
        row = self.writer.get_row(gv)
        self.assertEqual(
            str(row),'chr24\t2094\tvrtid\tTGG\tCCC\t100\tLOWQUAL\t.')
        
        vrt = factory.from_edit('chr20', 1253922,'TGT','G')
        # test vcf notations TODO
        
        row = self.writer.get_row(vrt)
        self.assertEqual(str(row), 'chr20\t1253923\t.\tTGT\tG\t.\t.\t.')

        # vf = VariantFactory(reference=ref,
        #                     normindel=True)
        vrt = factory.from_edit('chr1',13957,'TCCCCCA','TCCCCA')
        with self.assertRaises(ValueError) as cm:
            row = self.writer.get_row(vrt)
        self.assertIn('Reference is required',cm.exception.args[0])

    def test_change_of_attributes(self):
        reader = VCFReader(
            pkg_file('genomvar.test','data/example1.vcf'))
        vrt = list(reader.iter_vrt())[0]
        self.assertEqual(str(self.writer.get_row(vrt)),
                         'chr24\t23\t1\tAG\tA\t100\tPASS\t.')
        vrt2 = copy.deepcopy(vrt)
        vrt2.attrib['id'] = '.'
        vrt2.attrib['qual'] = '.'
        vrt2.attrib['filter'] = '.'
        self.assertEqual(str(self.writer.get_row(vrt2)),
                         'chr24\t23\t.\tAG\tA\t.\t.\t.')
        self.assertEqual(str(self.writer.get_row(vrt2, id='.', qual='.', filter='.')),
                         'chr24\t23\t.\tAG\tA\t.\t.\t.')
        
        vrt3 = copy.deepcopy(vrt)
        vrt3.attrib['id'] = None
        vrt3.attrib['qual'] = None
        vrt3.attrib['filter'] = None
        self.assertEqual(str(self.writer.get_row(vrt3)),
                         'chr24\t23\t.\tAG\tA\t.\t.\t.')
        reader.close()
        reader.close()

    def test_to_vcf_row_from_file(self):
        def _split_multiallelic(rows):
            for row in rows:
                for alt in row.ALT.split(','):
                    kwds = {f:getattr(row,f) for f in VCF_FIELDS}
                    kwds['ALT'] = alt
                    kwds['INFO'] = '.'
                    kwds['FORMAT'] = None
                    kwds['SAMPLES'] = None
                    yield str(VCFRow(**kwds))

        reader = VCFReader(pkg_file('genomvar.test','data/example1.vcf'))
        variants = list(reader.iter_vrt(
            parse_info=False,parse_samples=False))
        rows = [str(self.writer.get_row(v)) for v in variants]
        
        for r1, r2 in zip(
                _split_multiallelic(reader.iter_rows()), rows):
            if 'AG\tAGG' in r1: # stripping 
                continue
            self.assertEqual(r1,r2)
        reader.close()

    def test_to_vcf_row_instantiated_variant_numeric_chrom(self):
        factory = variant.VariantFactory()
        v1 = factory.from_edit(1,2093,'TGG','CCC')
        row = str(self.writer.get_row(v1))
        row2 = str(v1)
        self.assertIn('TGG', row)
        self.assertIn('TGG', row2)
