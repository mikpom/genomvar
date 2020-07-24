from unittest import TestCase
import tempfile
import io
import operator
import itertools
import re
from pkg_resources import resource_filename as pkg_file
import numpy as np
from genomvar.varset import MutableVariantSet,\
    VariantFileSet,IndexedVariantFileSet,VariantSet
from genomvar.variant import VariantBase,AmbigIndel,Haplotype,VariantFactory
from genomvar import varset,Reference
from genomvar import OverlappingHaplotypeVars,\
    UnsortedVariantFileError,VCFSampleMismatch,\
    DuplicateVariants,DifferentlySortedChromsError
from genomvar.vcf import VCFRow,VCFReader
from genomvar import variant

CHR1RGN = pkg_file('genomvar.test','data/chr1rgn.fasta')
CHR15RGN = pkg_file('genomvar.test','data/chr15rgn.fna')
test_vcf1 = pkg_file('genomvar.test','data/test_vcf.vcf')
test_vcf8 = pkg_file('genomvar.test','data/test_vcf8.vcf')
test_vcf_1kg = pkg_file('genomvar.test','data/test_vcf_1kg.vcf.gz')

# REF      AG    T|   T|  C    G     TGG   TT    G|      T    T      CACAGTTCCAC
#          22    154  165 453  1206  2093  2099  3200    4754 6145   10044
# varset1  A     TT   TG  CT   C     CCC   GG    GG      TCG  C      T----------
#          AGG   TT   TG       T     CCC   GG    GG      T    G      G----------
#          AG                        --------            ------ phased
# FILT           h    h;n                  n     n

factory = variant.VariantFactory()

class TestVariantSetCase(TestCase):
    def test_empty_vcf(self):
        vs = VariantSet.from_vcf(pkg_file('genomvar.test',
                                          'data/empty_vcf.vcf'))
        self.assertEqual(vs.nof_unit_vrt(),0)
    
    def test_sort_chroms(self):
        vs = VariantSet.from_vcf(pkg_file('genomvar.test',
                                          'data/test_vcf2.vcf'))
        vs.sort_chroms()
        self.assertEqual(list(vs.get_chroms()),['chr11rgn','chr15rgn'])

        vs.sort_chroms(key=lambda c: 1 if c=='chr15rgn' else 2)
        self.assertEqual(list(vs.get_chroms()),['chr15rgn','chr11rgn'])

    def test_diff_del(self):
        vrt = variant.Del(chrom="chr1",start=6751613,end=6751627)
        vs1 = varset.VariantSet.from_variants([vrt])
        vs2 = varset.VariantSet.from_variants([vrt])
        self.assertEqual(len(list(vs1.diff(vs2).iter_vrt())),0)
        self.assertEqual(len(list(vs1.comm(vs2).iter_vrt())),1)

    def test_diff_ins(self):
        vrt = variant.Ins(chrom="chr1",start=6751613,alt='AGTC')
        vs1 = varset.VariantSet.from_variants([vrt])
        vs2 = varset.VariantSet.from_variants([vrt])
        self.assertEqual(len(list(vs1.diff(vs2).iter_vrt())),0)
        self.assertEqual(len(list(vs1.comm(vs2).iter_vrt())),1)
        
    def test_ambig_difference_different_ambig(self):
        #         10043
        # R       T--CA--CA--G
        # v1 VCF  T--CA--CACAG ins CA right
        # v2 VCF  T------CA--G del CA left
        # v1      TCACA--CA--G ins CA left
        # v2      T--CA------G del CA right
        s1 = VariantSet.from_vcf(pkg_file('genomvar.test',
                                          'data/test_vcf1_ambig2.vcf'),
                                 reference=CHR15RGN,normindel=True)
        fac = variant.VariantFactory(reference=CHR15RGN,normindel=True)
        s2 = VariantSet.from_variants(
            [fac.from_edit('chr15rgn',10043,'T','TCA'),
             fac.from_edit('chr15rgn',10045,'ACA','A')])
        self.assertEqual(len(list(s1.diff(s2,match_ambig=True).iter_vrt())),0)

    def test_strip_order_dependent_Ambig(self):
        #    10043
        # R  T--CA--CAG
        # v1 TCACA--CAG
        # v2 T--CACACAG
        factory = variant.VariantFactory(reference=CHR15RGN,
                                         normindel=True)
        v1 = factory.from_edit('chr15rgn',10043,'T','TCA')
        v2 = factory.from_edit('chr15rgn',10045,'A','ACA')
        s1 = VariantSet.from_variants([v1])
        s2 = VariantSet.from_variants([v2])

        diff = s1.diff(s2,match_ambig=True)
        self.assertEqual(len(list(diff.iter_vrt())),0)
        diff = s1.diff(s2,match_ambig=False)
        self.assertEqual(len(list(diff.iter_vrt())),1)

    def test_ovlp(self):
        #                        23
        # TTCACTTAGCATAATGTCTTCAAG|ATT
        #                         G
        #                          C <- not interfering

        ins = factory.from_edit(chrom='chr15rgn',start=23,ref='G',
                                    alt='GG')
        s1 = VariantSet.from_variants([ins],)
        vb = factory.from_edit(chrom='chr15rgn',start=24,ref='A',alt='C')
        self.assertEqual(len(s1.ovlp(vb)),0)
        vb2 = factory.from_edit(chrom='chr15rgn',start=23,ref='G',alt='GG')
        self.assertEqual(len(s1.ovlp(vb2,match_ambig=False)),1)

    def test_from_variants(self):
        vfset = VariantFileSet(test_vcf1)
        vset = VariantSet.from_variants(list(vfset.iter_vrt()))
        vrt = list(vset.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(vrt),2)

        # Test error on out of reference bounds
        with self.assertRaises(ValueError):
            VariantSet.from_variants(
                list(vfset.iter_vrt())+[variant.SNP('chr15rgn',10000000,'C')],
                reference=Reference(CHR15RGN))

        # Test error on chromosome not in reference
        with self.assertRaises(ValueError):
            vs = VariantSet.from_variants(
                list(vfset.iter_vrt())+[variant.SNP('chr2',10,'C')],
                reference=Reference(CHR15RGN))
            print(vs.chroms)

    def test_from_vcf(self):
        vset = VariantSet.from_vcf(test_vcf1)

        self.assertEqual(vset.nof_unit_vrt(),18)
        self.assertEqual(list(vset.chroms),['chr15rgn'])
        vrt = list(vset.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertTrue(v1.is_instance(variant.SNP))
        self.assertEqual(v1.attrib['vcf_notation']['ref'],'G')
        self.assertEqual(v1.attrib['vcf_notation']['start'],1206)
        self.assertEqual(v1.attrib['vcf_notation']['id'],'5')
        self.assertEqual(v1.attrib['vcf_notation']['filt'],'PASS')
        vrt = list(vset.find_vrt('chr15rgn',154,156))
        self.assertEqual(len(vrt),1)
        v1 = vrt[0]
        self.assertTrue(v1.is_instance(variant.Ins))
        vrt = list(vset.find_vrt('chr15rgn',20,25))
        self.assertEqual(len(vrt),2)
        self.assertEqual(len(list(vset.iter_vrt())),16)

    def test_from_vcf_with_info(self):
        vset = VariantSet.from_vcf(test_vcf1,parse_info=True)

        vrt = list(vset.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.alt,'C')
        self.assertEqual(v1.attrib['info']['NSV'],1)

    def test_from_vcf_with_samples(self):
        vset = VariantSet.from_vcf(test_vcf1,
                                   parse_samples=True)

        vrt = list(vset.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(1,0))
        
    def test_asterisk_variant(self):
        vset = VariantSet.from_vcf(pkg_file('genomvar.test',
                            'data/test_vcf11_asterisk.vcf'),
                                   parse_info=True)

        vrt = list(vset.find_vrt('chr1',995507,995515))
        self.assertEqual(len(vrt),3)
        # v1,v2 = vrt
        # self.assertEqual(v1.alt,'TTTTT')
        # self.assertEqual(v2.alt,'C')

    def test_find_vrt_chrom_only(self):
        s1 = VariantSet.from_vcf(pkg_file('genomvar.test','data/test_vcf2.vcf'),
                                 parse_info=True,parse_samples=True)
        self.assertEqual(len(list(s1.find_vrt('chr15rgn'))),4)


    def test_match(self):
        # REF      TGG   TT    
        #          2093  2099  
        # vs1      CCC   GG
        # vs2            CG    
        #          r1   r2,r3
        vs1 = VariantSet.from_vcf(test_vcf1)
        vrt_CG = factory.from_edit('chr15rgn',2098,'TT','CG')
        vrt_CC = factory.from_edit('chr15rgn',2098,'TT','CC')
        self.assertEqual(len(vs1.match(vrt_CG)),1)
        self.assertEqual(len(vs1.match(vrt_CG,match_partial=False)),0)
        self.assertEqual(len(vs1.match(vrt_CC)),0)
        self.assertEqual(len(vs1.match(vrt_CC,match_partial=False)),0)

    def test_drop_duplicates2(self):
        with open(pkg_file('genomvar.test','data/test_vcf.vcf')) as fh:
            header = []
            variants = []
            for line in fh:
                if line.startswith('#'):
                    header.append(line)
                else:
                    variants.append(line)
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')
        with open(tf.name,'wt') as fh:
            for line in header:
                fh.write(line)
            for i in range(2):
                for line in variants:
                    fh.write(line)
        vs = VariantSet.from_vcf(tf.name)
        vs2,dropped = vs.drop_duplicates(return_dropped=True)
        self.assertEqual(vs.nof_unit_vrt()/2,vs2.nof_unit_vrt())
        self.assertEqual(vs2.nof_unit_vrt(),
                         sum([v.nof_unit_vrt() for v in dropped]))

class TestIndexedVariantFileCase(TestCase):
    vcf1_bgz = pkg_file('genomvar.test','data/test_vcf.vcf.gz')
    vcf4_bgz = pkg_file('genomvar.test','data/test_vcf4.vcf.gz')

    def test_complex_INFO_example(self):
        vset = IndexedVariantFileSet(self.vcf4_bgz,parse_info=True)
        checked = False
        for vrt in vset.find_vrt(rgn='chr1:69090-69091'):
            if vrt.alt!='C':
                continue
            if not vrt.attrib['info']['MutPred_Top5features'] is None:
                checked = True
                self.assertTrue(vrt.attrib['info']['MutPred_Top5features']\
                                .startswith('Loss of sheet (P = 0.0817)| L'))
                
        self.assertTrue(checked)

    def test_find_vrt(self):
        ivfs = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/test_vcf2.vcf.gz'))
        vs = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/test_vcf2.vcf.gz'))

        self.assertEqual(
            sum([v.nof_unit_vrt() for v in ivfs.find_vrt('chr15rgn')]),
            sum([v.nof_unit_vrt() for v in vs.find_vrt('chr15rgn')]))

        self.assertEqual(
            sum([v.nof_unit_vrt() for v in ivfs.iter_vrt()]),
            sum([v.nof_unit_vrt() for v in vs.iter_vrt()]))
        
    def test_many_samples(self):
        vcf  = pkg_file('genomvar.test','data/test_vcf_1kg.vcf.gz')
        vset = IndexedVariantFileSet(vcf)
        checked = False
        for vrt in vset.find_vrt(rgn='17:120339-120340'):
            checked = True
                
        self.assertTrue(checked)

    def test_find_nonexistent_chrom(self):
        vcf  = pkg_file('genomvar.test','data/test_vcf_1kg.vcf.gz')
        vset = IndexedVariantFileSet(vcf)
        self.assertEqual(list(vset.find_vrt('chrZ')),[])
        
    def test_match(self):
        # REF      TGG   TT    
        #          2093  2099  
        # vs1      CCC   GG
        # vrt            CG    
        #          r1   r2,r3
        for ref in CHR15RGN,None: # should work with or without ref
            vs1 = IndexedVariantFileSet(self.vcf1_bgz)
            vrt = factory.from_edit('chr15rgn',2098,'TT','CG')
            self.assertEqual(len(vs1.match(vrt)),1)

            # Test insertion
            vrt = factory.from_edit('chr15rgn',22,'AG','AGG')
            match = vs1.match(vrt)
            self.assertEqual(len(match),1)

    def  test_error_on_format(self):
        with self.assertRaises(OSError):
            vset = IndexedVariantFileSet(test_vcf1,reference=CHR15RGN)

    def test_some_variants(self):
        vset = IndexedVariantFileSet(self.vcf1_bgz,reference=CHR15RGN)
        self.assertEqual(len(list(vset.find_vrt(rgn='chr15rgn:1200-1210'))),2)
        v1,v2 = list(vset.find_vrt(rgn='chr15rgn:1200-1210'))
        self.assertEqual([v1.start,v1.end],[1206,1207])
        self.assertEqual([v2.start,v2.end],[1206,1207])

        self.assertEqual(len(list(vset.find_vrt(rgn='chr15rgn:3200-3205'))),1)
        v1 = list(vset.find_vrt(rgn='chr15rgn:3200-3205'))[0]
        self.assertEqual([v1.start,v1.end],[3201,3202])
        
        self.assertEqual(len(list(vset.find_vrt(rgn='chr15rgn:20-30'))),2)
        v1,v2 = list(vset.find_vrt(rgn='chr15rgn:20-30'))
        self.assertEqual([v1.start,v1.end,v1.vtp],
                         [23,24,variant.Del])
        self.assertEqual([v2.start,v2.end,v2.vtp],
                         [24,25,variant.Ins])

    def test_wrong_chrom_name_in_ref(self):
        vset = IndexedVariantFileSet(self.vcf1_bgz,reference=CHR1RGN)
        self.assertEqual(len(list(vset.find_vrt(rgn='chr15rgn:1200-1210'))),2)

    def test_class(self):
        vset = IndexedVariantFileSet(self.vcf1_bgz,parse_info=True,
                                     reference=CHR15RGN,parse_samples='SAMP1')

        # Test find_vrt and returned INFO
        vrt = list(vset.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.attrib['info']['NSV'],1)
        self.assertEqual(v2.attrib['info']['RECN'],19)

        # Test multiallelic
        vrt = list(vset.find_vrt('chr15rgn',20,30))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.attrib['info']['AF'],0.5)
        self.assertEqual(v2.attrib['info']['AF'],0.5)

        # Test find_vrt precision
        vrt = list(vset.find_vrt('chr15rgn',2095,2096))
        self.assertEqual(len(vrt),1)
        vrt = list(vset.find_vrt('chr15rgn',2098,2100))
        self.assertEqual(len(vrt),1)

        # Test find all variants
        self.assertEqual(len(list(vset.find_vrt())),16)

        # Test finding all variants
        self.assertEqual(len(list(vset.find_vrt())),16)
        
    def test_ctg_len_without_ref(self):
        vset = IndexedVariantFileSet(self.vcf1_bgz,
                        parse_samples='SAMP1')
        self.assertEqual(vset.chroms,{'chr15rgn'})
        
class TestVariantFileCase(TestCase):
    def test_unsorted_VCF_input(self):
        header = []
        lines = []
        with open(pkg_file('genomvar.test','data/test_vcf.vcf'),'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    header.append(line)
                else:
                    lines.append(line)
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')
        with open(tf.name,'wt') as fh:
            fh.writelines(header)
            fh.writelines(reversed(lines))

        vs1 = VariantFileSet(pkg_file('genomvar.test','data/test_vcf.vcf'))
        vs2 = VariantFileSet(tf.name)
        with self.assertRaises(UnsortedVariantFileError):
            list(vs1.diff_vrt(vs2).iter_vrt())

class MutableVariantSetTestCase(TestCase):
    def test_sort_chroms(self):
        vs = MutableVariantSet.from_vcf(pkg_file('genomvar.test',
                                                 'data/test_vcf2.vcf'),max_ploidy=4)
        vs.sort_chroms()
        self.assertEqual(list(vs.get_chroms()),['chr11rgn','chr15rgn'])

        vs.sort_chroms(key=lambda c: 1 if c=='chr15rgn' else 2)
        self.assertEqual(list(vs.get_chroms()),['chr15rgn','chr11rgn'])
        
        
    def test_variant_addition(self):
        # **Insertion deletion nearby**
        #                        23
        # TTCACTTAGCATAATGTCTTCAAG|ATT
        # v1                   AA-|ATT
        # v2                   AAGGATT
        s1 = MutableVariantSet(reference=CHR15RGN)
        vfac = s1.get_factory(normindel=True)
        vb = vfac.from_edit('chr15rgn',22,'AG','A')
        v1 = s1.add_vrt(vb,GT=(1,0))
        self.assertEqual(type(v1.base),variant.Del)
        self.assertEqual([v1.start,v1.end,v1.ref,v1.alt],
                         [23,24,'G',''])
        # Adding insertion of G
        vb = vfac.from_edit('chr15rgn',23,'G','GG')
        v2 = s1.add_vrt(vb,GT=(0,1))
        self.assertEqual(type(v2.base),variant.AmbigIns)
        self.assertEqual([v2.start,v2.end],[23,25])
        self.assertEqual([v2.act_start,v2.act_end],[24,25])
        self.assertEqual(len(list(s1.find_vrt('chr15rgn',20,30))),2)

    def test_nof_unit_vrt(self):
        # REF      TGG   TT    G|
        #          2093  2099  3200
        # varset1  CCC   GG
        #                CG    GG
        #          v1   v2,v3   r4
        variants = [0]*4
        variants[0] = factory.from_edit('chr15rgn',2093,'TGG','CCC')
        variants[1] = factory.from_edit('chr15rgn',2098,'TT','GG')
        variants[2] = factory.from_edit('chr15rgn',2098,'TT','CG')
        variants[3] = factory.from_edit('chr15rgn',3200,'G','GG')
        s1 = MutableVariantSet(reference=CHR15RGN)
        for vrt in variants:
            try:
                s1.add_vrt(vrt,GT=(1,0))
            except OverlappingHaplotypeVars:
                s1.add_vrt(vrt,GT=(0,1))
        self.assertEqual(s1.nof_unit_vrt(),7)

    def test_nof_unit_vrt2(self):
        vs = MutableVariantSet.from_vcf(
            pkg_file('genomvar.test','data/test_vcf7.vcf'))
        self.assertEqual(vs.nof_unit_vrt(),1000)
        
    def test_iter_vrt_by_chrom(self):
        vset = MutableVariantSet.from_vcf(test_vcf1,
                        reference=CHR15RGN,sample='SAMP1')

        chroms = {}
        for chrom,it in vset.iter_vrt_by_chrom():
            chroms[chrom] = len(list(it))

        self.assertEqual(set(chroms),{'chr15rgn'})
    def test_mixed_variant(self):
        #                        23
        # TTCACTTAGCATAATGTCTTCAAGATT
        #                       AT-TT
        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = factory.from_edit('chr15rgn',23,'GA','T')
        v1 = s1.add_vrt(vb,GT=(0,1),attrib={'info':{'f1':'1'}})
        self.assertTrue(type(v1),variant.Mixed)
        self.assertEqual([v1.start,v1.end,v1.ref,v1.alt,type(v1.base)],
                         [23,25,'GA','T',variant.Mixed])
        self.assertEqual(v1.attrib['info']['f1'],'1')
    
    def test_ovlp(self):
        #                        23
        # TTCACTTAGCATAATGTCTTCAAG|ATT
        #                         G
        #          interfering-> C C <- not interfering
        s1 = MutableVariantSet(reference=CHR15RGN)
        vfac = s1.get_factory(normindel=True)
        s1.add_vrt(vfac.from_edit('chr15rgn',23,'G','GG'),
                   GT=(0,1))
        vb = vfac.from_edit(chrom='chr15rgn',start=24,ref='A',alt='C')
        self.assertEqual(len(s1.ovlp(vb)),0)
        vb2 = vfac.from_edit(chrom='chr15rgn',start=23,ref='G',alt='C')
        self.assertEqual(len(s1.ovlp(vb2,match_ambig=False)),0)
        self.assertEqual(len(s1.ovlp(vb2,match_ambig=True)),1)

        # Same reversed
        s1 = MutableVariantSet(reference=CHR15RGN)
        s1.add_vrt(vfac.from_edit('chr15rgn',23,'G','C'),
                   GT=(0,1))
        s1.add_vrt(vfac.from_edit('chr15rgn',24,'A','C'),
                   GT=(0,1))
        vb = vfac.from_edit(chrom='chr15rgn',start=23,ref='G',alt='GG')
        self.assertEqual(len(s1.ovlp(vb)),0)
        self.assertEqual(len(s1.ovlp(vb,match_ambig=True)),1)

        # ** Test insertion interference **
        #                        23
        # TTCACTTAGCATAATGTCTTCAAG|ATT
        #                         G
        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = vfac.from_edit('chr15rgn',23,'G','GG')
        s1.add_vrt(vb,GT=(0,1))
        vb = vfac.from_edit(chrom='chr15rgn',start=23,
                                   ref='G',alt='GG')
        self.assertEqual(len(s1.ovlp(vb,match_ambig=False)),1)
        self.assertEqual(len(s1.ovlp(vb,match_ambig=True)),1)
        vb2 = vfac.from_edit(chrom='chr15rgn',start=23,ref='G',alt='GT')
        self.assertEqual(len(s1.ovlp(vb,match_ambig=False)),1)
        self.assertEqual(len(s1.ovlp(vb,match_ambig=True)),1)

        # ** Test ambig insertion and deletion non-interference **
        #                        23
        # TTCACTTAGCATAATGTCTTCAAG|ATT
        #                         G
        #                       AG -T
        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = vfac.from_edit('chr15rgn',23,'G','GG')
        s1.add_vrt(vb,GT=(0,1))
        vb = vfac.from_edit(chrom='chr15rgn',start=23,ref='GA',alt='G')
        self.assertEqual(len(s1.ovlp(vb,match_ambig=False)),0)
        self.assertEqual(len(s1.ovlp(vb,match_ambig=True)),0)

        # ** Test insertion and deletion non-interference **
        #                        23
        # TTCACTTAGCATAATGTCTTCAAG|ATT
        #                         C
        #                       AG -T
        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = vfac.from_edit('chr15rgn',23,'G','GC')
        s1.add_vrt(vb,GT=(0,1))
        vb = vfac.from_edit(chrom='chr15rgn',start=23,ref='GA',alt='G')
        self.assertEqual(len(s1.ovlp(vb,match_ambig=False)),0)
        self.assertEqual(len(s1.ovlp(vb,match_ambig=True)),0)

        #            11
        # TTCACTTAGCATAATGTC
        #            T-AT
        #              C
        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = vfac.from_edit('chr15rgn',11,'TA','T')
        s1.add_vrt(vb,GT=(0,1))
        vb = vfac.from_edit(chrom='chr15rgn',start=13,ref='A',alt='C')
        self.assertEqual(len(s1.ovlp(vb,match_ambig=False)),0)
        self.assertEqual(len(s1.ovlp(vb,match_ambig=True)),1)

    def test_match_func(self):
        # REF      TGG   TT    G|
        #          2093  2099  3200
        # varset1  CCC   GG
        #                CC    GG
        #          r1   r2,r3   r4
        r1 = factory.from_edit('chr15rgn',2093,'TGG','CCC')
        r2 = factory.from_edit('chr15rgn',2098,'TT','GG')
        r3 = factory.from_edit('chr15rgn',2098,'TT','CC')
        r4 = factory.from_edit('chr15rgn',3200,'GG','G')

        s1 = MutableVariantSet(reference=CHR15RGN)
        hap1 = s1.add_hap_variants([r1,r2],GT=(0,1))
        ccc1 = sorted(hap1.variants,key=lambda o: o.start)[0]
        hap2 = s1.add_hap_variants([r3,r4],GT=(0,1),
                                   allow_adjust_genotype=True)

        s2 = MutableVariantSet(reference=CHR15RGN)
        hap = s2.add_hap_variants([r1,r4],GT=(0,1))
        ccc2 = sorted(hap.variants,key=lambda o: o.start)[0]

        match = s1.match(hap)
        self.assertEqual(len(match),1)
        self.assertEqual(match,{ccc2.key:[ccc1]})

    def test_copy_with_attr(self):
        kwargs = {'vcf':test_vcf1,
                  'reference':CHR15RGN,
                  'sample':'SAMP1',
                  'normindel':True,
                  'parse_info':True}
        s = MutableVariantSet.from_vcf(**kwargs)

        s2 = s.copy()
        _vrt = list(s2.find_vrt('chr15rgn',150,160))
        self.assertEqual(len(_vrt),1)
        v = _vrt[0]
        self.assertEqual(v.attrib['info']['AF'],1.0)


        
    def test_attributes(self):
        vset = MutableVariantSet()
        vb = factory.from_edit('chr15rgn',2092,'TT','AC')
        v = vset.add_vrt(vb,GT=(0,1),
                         attrib={'info':{'origin':'manual'}})
        self.assertEqual(v.attrib['info']['origin'],'manual')

class TestSetComparisonCase(TestCase):
    vfac = VariantFactory(reference=CHR15RGN)
    vfac_norm = VariantFactory(reference=CHR15RGN,normindel=True)
    vcf1_bgz = pkg_file('genomvar.test','data/test_vcf.vcf.gz')
    vcf2_bgz = pkg_file('genomvar.test','data/test_vcf2.vcf.gz')
    vcf4_bgz = pkg_file('genomvar.test','data/test_vcf4.vcf.gz')

    def test_diff_vrt(self):
        # REF    G     TTGG         C
        #        1206   2093        10044
        # s1     C      CCC         T
        #        G      CCC         G
        #                           C
        #              2092         10044
        # s2           TT           T
        #              AC           T

        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = self.vfac.from_edit('chr15rgn',1206,'G','C')
        s1.add_vrt(vb,GT=(1,0))
        vb = self.vfac.from_edit('chr15rgn',2093,'TG','CC')
        s1.add_vrt(vb,GT=(1,1))
        vb = self.vfac.from_edit('chr15rgn',2095,'G','C')
        s1.add_vrt(vb,GT=(1,1))
        vb = self.vfac.from_edit('chr15rgn',10044,'C','T')
        s1.add_vrt(vb,GT=(1,0))
        vb = self.vfac.from_edit('chr15rgn',10044,'C','G')
        s1.add_vrt(vb,GT=(0,1))

        s2 = MutableVariantSet(reference=CHR15RGN)
        vb = self.vfac.from_edit('chr15rgn',2092,'TT','AC')
        s2.add_vrt(vb,GT=(0,1),allow_adjust_genotype=True)
        vb = self.vfac.from_edit('chr15rgn',10044,'C','T')
        s2.add_vrt(vb,GT=(1,1),allow_adjust_genotype=True)

        diff = s1.diff(s2)
        # vset = diff.all_variants()
        self.assertEqual(diff.nof_unit_vrt(),4)
        
        v1,v2,v3,v4 = list(diff.iter_vrt())

        self.assertEqual([v1.start,v1.ref,v1.alt,v1.GT],
                         [1206,'G','C',(1,0)])
        # self.assertEqual([v2.start,v2.ref,v2.alt,v2.GT],
        #                  [2094,'GG','CC',(1,1)])
        self.assertEqual([v4.start,v4.ref,v4.alt,v4.GT],
                         [10044,'C','G',(0,1)])

    def test_diff_with_attr_filters(self):
        s1 = MutableVariantSet.from_vcf(test_vcf1,reference=CHR15RGN,
                                        normindel=True,parse_info=True,
                                        sample='SAMP1')
        _vrt = list(s1.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(_vrt),2)
        v1,v2 = _vrt
        self.assertEqual(v1.attrib['info']['AF'],1.0)

        s2 = MutableVariantSet(reference=CHR15RGN)
        s2.add_vrt(factory.from_edit('chr15rgn',1206,'G','C'),
                   GT=(0,1))

        diff = s1.diff(s2)
        _vrt = list(diff.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(_vrt),1)
        v = _vrt[0]
        self.assertEqual(v.attrib['filters'],['PASS'])
        self.assertEqual(v.attrib['info']['AF'],1.0)

    def test_strip_order_dependent_Ambig(self):
        #    10043
        #    TCACAG
        # s1 TC--AG
        # s2 T--CAG
        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = self.vfac_norm.from_edit('chr15rgn',10044,'CAC','C')
        s1.add_vrt(vb,GT=(0,1))
        s2 = MutableVariantSet(reference=CHR15RGN)
        vb = self.vfac_norm.from_edit('chr15rgn',10043,'TCA','T')
        s2.add_vrt(vb,GT=(0,1))

        diff = s1.diff(s2,match_ambig=True)
        self.assertEqual(len(list(diff.iter_vrt())),0)

        diff = s1.diff(s2,match_ambig=False)
        mset = list(diff.iter_vrt())
        self.assertEqual(len(mset),1)
        
    def test_diff_of_ambig_indel(self):
        #     9945
        #     CTTTTTCAT
        # s1  CTT--TCAT
        # s2  C--TTTCAT
        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = factory.from_edit('chr15rgn',9947,'TT','T')
        s1.add_vrt(vb,GT=(0,1))
        s2 = MutableVariantSet(reference=CHR15RGN)
        vb = factory.from_edit('chr15rgn',9945,'CTT','C')
        s2.add_vrt(vb,GT=(1,0))

        diff = s1.diff(s2)
        mset = list(diff.iter_vrt())
        self.assertEqual(len(mset),1) #TODO fix/reconsider

    def test_consistency_of_diff_and_com(self):
        # REF    G     TTGG         CACAGTTC---CA-C
        #        1206   2093        10044
        # s1     C      CCC         T       CC  G
        #        C      CCC         G       CC
        #              2092         10044
        # s2           TT           T        CC
        #              AC           T           G

        s1 = MutableVariantSet(reference=CHR15RGN)
        vfac = s1.get_factory(normindel=True)
        vb = vfac.from_edit('chr15rgn',1206,'G','C')
        s1.add_vrt(vb,GT=(1,1))
        vb = vfac.from_edit('chr15rgn',2093,'TG','CC')
        s1.add_vrt(vb,GT=(1,1))
        vb = vfac.from_edit('chr15rgn',2095,'G','C')
        s1.add_vrt(vb,GT=(1,1))
        vb = vfac.from_edit('chr15rgn',10044,'C','T')
        s1.add_vrt(vb,GT=(1,0))
        vb = vfac.from_edit('chr15rgn',10044,'C','G')
        s1.add_vrt(vb,GT=(0,1))
        vb = vfac.from_edit('chr15rgn',10051,'C','CCC')
        s1.add_vrt(vb,GT=(1,1))
        vb = vfac.from_edit('chr15rgn',10053,'A','AG')
        s1.add_vrt(vb,GT=(1,0))

        s2 = MutableVariantSet(reference=CHR15RGN)
        vfac = s2.get_factory(normindel=True)
        vb = vfac.from_edit('chr15rgn',2092,'TT','AC')
        s2.add_vrt(vb,GT=(0,1))
        vb = vfac.from_edit('chr15rgn',10044,'C','T')
        s2.add_vrt(vb,GT=(1,1))
        vb = vfac.from_edit('chr15rgn',10052,'C','CCC')
        s2.add_vrt(vb,GT=(1,0))
        vb = vfac.from_edit('chr15rgn',10053,'A','AG')
        s2.add_vrt(vb,GT=(0,1))

        com1 = s1.comm(s2,match_ambig=True)
        com2 = s2.comm(s1,match_ambig=True)
        s1_s2 = s1.diff(s2,match_ambig=True)
        s2_s1 = s2.diff(s1,match_ambig=True)

        nof = {}
        nof['com1'] = com1.nof_unit_vrt()
        nof['com2'] = com1.nof_unit_vrt()
        nof['s1'] = s1.nof_unit_vrt()
        nof['s2'] = s2.nof_unit_vrt()
        nof['s1_s2'] = s1_s2.nof_unit_vrt()
        nof['s2_s1'] = s2_s1.nof_unit_vrt()

        self.assertEqual(nof['com1'],4)
        self.assertEqual(nof['com2'],nof['com1'])
        self.assertEqual(nof['com1']+nof['s1_s2'],nof['s1'])
        self.assertEqual(nof['com1']+nof['s2_s1'],nof['s2'])

    def test_different_instatiation_of_variants(self):
        vset1 = MutableVariantSet()
        vset2 = MutableVariantSet()
        fac = VariantFactory()
        ins1 = fac.from_hgvs('chr15rgn:g.10053_10054insG')
        ins2 = fac.from_edit('chr15rgn',10052,'A','AG')

        del1 = fac.from_hgvs('NG_012232.1:g.19_21del')
        del2 = fac.from_edit('NG_012232.1',17,'GTTT','G')

        snp1 = fac.from_hgvs('chr1:g.15C>A')
        snp2 = fac.from_edit('chr1',14,'C','A')

        mixed1 = fac.from_hgvs('NC_000023.11:g.10delinsGA')
        mixed2 = fac.from_edit('NC_000023.11',9,'C','GA')

        for vrt in (ins1,del1,snp1,mixed1):
            vset1.add_vrt(vrt,GT=(1,0))
        for vrt in (ins2,del2,snp2,mixed2):
            vset2.add_vrt(vrt,GT=(1,0))
        self.assertEqual(len(list(vset1.iter_vrt())),4)
        self.assertEqual(len(list(vset2.iter_vrt())),4)

        self.assertEqual(len(vset2.match(ins1)),1)
        com = vset1.comm(vset2)
        self.assertEqual(len(list(com.iter_vrt())),4)

        diff = vset1.diff(vset2)
        self.assertEqual(len(list(diff.iter_vrt())),0)

    def test_variant_cluster(self):
        vs1 = varset.VariantSet.from_variants(
            [variant.SNP('chr1',13366968,'A'),
             variant.SNP('chr1',13366969,'T'),
             variant.SNP('chr1',13366970,'G')])

        vs2 = varset.VariantSet.from_variants(
            [variant.Del('chr1',13366967,13366969),
             variant.Ins('chr1',13366971,'TG')])

        diff = list(vs1.diff(vs2).iter_vrt())
        self.assertEqual(len(diff), 3)

        vs3 = varset.VariantSet.from_variants(
            [variant.MNP("chr1",13366968,'ATG')])
        self.assertEqual(len(list(vs1.diff(vs3,match_partial=True)\
                                     .iter_vrt())), 0)
        self.assertEqual(len(list(vs3.diff(vs1,match_partial=True)\
                                     .iter_vrt())), 0)
        self.assertEqual(len(list(vs1.diff(vs3,match_partial=False)\
                                     .iter_vrt())), 3)
        self.assertEqual(len(list(vs3.diff(vs1,match_partial=False)\
                                     .iter_vrt())), 1)
        

    def test_mnp_com_split(self):
        #                           23
        #    TTCACTTAGCATAATGTCTTCAAGATT
        # v1                       TT -single
        # v2                       TT -splitted
        vset1 = MutableVariantSet(reference=CHR15RGN)
        vb = factory.from_edit('chr15rgn',22,'AG','TT')
        vset1.add_vrt(vb,GT=(1,0))
        vset2 = MutableVariantSet(reference=CHR15RGN)
        vb = factory.from_edit('chr15rgn',22,'A','T')
        vset2.add_vrt(vb,GT=(0,1))
        vb = factory.from_edit('chr15rgn',23,'G','T')
        vset2.add_vrt(vb,GT=(0,1))

        com = vset1.comm(vset2)
        v = list(com.iter_vrt())[0]

        self.assertEqual([v.start,v.ref],[22,'AG'])
        self.assertFalse(list(vset1.diff(vset2).iter_vrt()))


    def test_in_the_middle(self):
        #      2
        #    TTCACTTAGCAT
        # v1   GGG
        # v2    G
        s1 = MutableVariantSet(reference=CHR15RGN)
        vb = factory.from_edit('chr15rgn',2,'CAC','GGG')
        s1.add_vrt(vb,GT=(1,0))
        s2 = MutableVariantSet(reference=CHR15RGN)
        vb = factory.from_edit('chr15rgn',3,'A','G')
        s2.add_vrt(vb,GT=(1,0))

        diff = s1.diff(s2)
        vv = list(diff.iter_vrt())
        self.assertEqual(len(vv),2)
        v1,v2 = vv
        self.assertEqual([v1.start,v1.ref,v1.alt],[2,'C','G'])
        self.assertEqual([v2.start,v2.ref,v2.alt],[4,'C','G'])

        self.assertEqual(len(list(s1.diff(s2,match_partial=False)
                                  .iter_vrt())),1)

    def test_diff_without_ref(self):
        vset1 = MutableVariantSet()
        vb = factory.from_edit('chr15rgn',10053,'A','AG')
        vset1.add_vrt(vb,GT=(1,0))

        vset2 = MutableVariantSet()
        vb = factory.from_edit('chr15rgn',10053,'A','AG')
        vset2.add_vrt(vb,GT=(0,1))

        com = vset1.comm(vset2)
        self.assertEqual(len(list(com.iter_vrt())),1)

        diff = vset1.diff(vset2)
        self.assertEqual(len(list(diff.iter_vrt())),0)

    def test_comm_vs_self_different_types(self):
        vcf = pkg_file('genomvar.test','data/test_vcf.vcf')
        vs1 = VariantSet.from_vcf(vcf)
        vs2 = MutableVariantSet.from_vcf(vcf)
        vs3 = VariantFileSet(vcf)

        for meth in 'comm','comm_vrt':
            for _vs1,_vs2 in itertools.product([vs1,vs2,vs3],repeat=2):
                method = operator.methodcaller(meth,_vs2)
                if meth=='comm':
                    if isinstance(_vs1,VariantFileSet): # method doesn't exist
                        continue
                    elif isinstance(_vs2,VariantFileSet):
                        with self.assertRaises(TypeError):
                            _vs1.comm(_vs2)
                        continue
                    elif isinstance(_vs1,VariantSet) and \
                              isinstance(_vs2,MutableVariantSet):
                        with self.assertRaises(TypeError):
                            _vs1.comm(_vs2)
                        continue
                
                comm = sorted(method(_vs1).iter_vrt(),key=lambda v: v.key)
                known = sorted(vs1.iter_vrt(),key=lambda v: v.key)
                self.assertEqual(len(comm),len(known))
                self.assertEqual(len(comm),len(known))
                for ind,vrt in enumerate(comm):
                    self.assertTrue(vrt.edit_equal(known[ind]))

    def test_cmp_stream(self):
        s1 = VariantFileSet(pkg_file('genomvar.test','data/test_vcf.vcf'))
        s2 = VariantFileSet(pkg_file('genomvar.test','data/test_vcf2.vcf'))
        nofv = 0
        for vrt in s1.diff_vrt(s2).iter_vrt():
            nofv += vrt.nof_unit_vrt()
        self.assertEqual(nofv,14)

    def test_diff_callback(self):
        s1 = VariantFileSet(pkg_file('genomvar.test','data/test_vcf5.vcf'))
        s2 = VariantFileSet(pkg_file('genomvar.test','data/test_vcf5.vcf'))
        cb = lambda m: [v.attrib['vcf_notation']['row'] for v in m]
        for N,vrt in enumerate(s1.comm_vrt(s2).iter_vrt(callback=cb)):
            self.assertEqual(vrt.attrib['vcf_notation']['row'],
                             vrt.attrib['cmp'][0])
        self.assertEqual(N,7)

    def test_differently_sorted_chroms(self):
        s1 = VariantFileSet(pkg_file(
            'genomvar.test','data/test_vcf5.vcf'))
        s2 = VariantFileSet(pkg_file(
            'genomvar.test','data/test_vcf6.vcf.gz'))
        with self.assertRaises(DifferentlySortedChromsError):
            list(s1.diff_vrt(s2).iter_vrt())

    def test_cmp_vrt_iter_same(self):
        vs = VariantFileSet(self.vcf2_bgz)
        tot = list(vs.find_vrt())
        # print(tot)
        comm = list(vs.comm_vrt(vs).iter_vrt())
        self.assertEqual(len(comm),len(tot))
    
    def test_cmp_vrt_iter_vrt(self):
        vs1 = VariantFileSet(self.vcf1_bgz,parse_samples=True)
        vs2 = VariantFileSet(self.vcf2_bgz,parse_samples=True)
        comm = list()
        for vrt in vs1.comm_vrt(vs2).iter_vrt():
            comm.append(vrt)
            self.assertTrue(vrt.attrib['samples'],
                            msg='Vrt {} has no samples'.format(vrt))
        self.assertEqual(len(comm),4)
        diff = vs1.diff_vrt(vs2).iter_vrt()
        self.assertEqual(len(list(diff)),12)

    def test_cmp_vrt_iter_vrt2(self):
        vs1 = VariantFileSet(test_vcf1)
        vs2 = VariantFileSet(test_vcf8)
        self.assertEqual(len(list(vs1.diff_vrt(vs2).iter_vrt())),
                         len(list(vs1.iter_vrt())))
        
    def test_cmp_vrt_region(self):
        vs1 = IndexedVariantFileSet(self.vcf1_bgz,parse_samples=True,
                                    parse_info=True)
        vs2 = IndexedVariantFileSet(self.vcf2_bgz,parse_samples='SAMP1',
                                    parse_info=True)
        comm = list(vs1.comm_vrt(vs2).region(rgn='chr15rgn:10040-10050'))
        self.assertEqual(len(comm),2)
        v1,v2 = comm
        self.assertEqual(v1.attrib['info']['AF'],1.0)
        self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(0,1))

    def test_ValueError_on_nonindexed(self):
        vs1 = IndexedVariantFileSet(pkg_file('genomvar.test',
                                             'data/test_vcf9.vcf.gz'),
                                    parse_samples=True)
        vs2 = VariantFileSet(pkg_file('genomvar.test',
                                             'data/test_vcf10.vcf.gz'),
                                    parse_samples=True)
        with self.assertRaises(ValueError) as cm:
            list(vs1.comm_vrt(vs2).region(rgn='7:152134922-152436005'))
        self.assertIn('region access',cm.exception.args[0])


    def test_cmp_vrt_region_multisample(self):
        vs1 = IndexedVariantFileSet(test_vcf_1kg,parse_samples=True,
                                    parse_info=True)
        vs2 = MutableVariantSet()
        vs2.add_vrt(factory.from_edit('17',120359,'A','C'),GT=(1,0))
        comm = list(vs1.comm_vrt(vs2).region(rgn='17:120350-120360'))
        self.assertEqual(len(comm),1)
        v1 = comm[0]
        # self.assertEqual(v1.attrib['info']['AF'],1.0)
        # self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(0,1))

    def test_cmp_vrt_region_multisample2(self):
        vs1 = IndexedVariantFileSet(pkg_file('genomvar.test',
                                             'data/test_vcf9.vcf.gz'),
                                    parse_samples=True)
        vs2 = IndexedVariantFileSet(pkg_file('genomvar.test',
                                             'data/test_vcf10.vcf.gz'),
                                    parse_samples=True)
        comm = []
        for vrt in vs2.comm_vrt(vs1).region(rgn='7:152134922-152436005'):
            comm.append(vrt)
            self.assertTrue(hasattr(vrt,'attrib'),msg='False for'+str(vrt))
        comm = list(vs2.comm_vrt(vs1).region(rgn='7:152134922-152436005'))
        self.assertGreater(len(comm),0)
            
    def test_haplotypes(self):
        # REF      TGG   TT    G-
        #          2093  2099  3200
        # varset1  CCC   GG
        #                CC    GG
        #          v1   v2,v3   v4
        v1 = factory.from_edit('chr15rgn',2093,'TGG','CCC')
        v2 = factory.from_edit('chr15rgn',2098,'TT','GG')
        v3 = factory.from_edit('chr15rgn',2098,'TT','CC')
        v4 = factory.from_edit('chr15rgn',3200,'G','GG')
        s1 = MutableVariantSet(reference=CHR15RGN)
        s1.add_hap_variants([v1,v2],GT=(1,0))
        s1.add_vrt(v3,GT=(0,1))
        s1.add_vrt(v4,GT=(0,1))
        s2 = VariantSet.from_variants(s1.iter_vrt())
        vrt = list(s2.find_vrt('chr15rgn',2090,2095,expand=True))
        self.assertEqual(len(vrt),2)
        self.assertEqual(s2.nof_unit_vrt(),8)
        self.assertEqual(s1.diff(s2).nof_unit_vrt(), 0)
        self.assertEqual(len(list(s2.diff_vrt(s1).iter_vrt())), 0)

        h1 = Haplotype.from_variants([v1,v2])
        h2 = Haplotype.from_variants([v3,v4])
        s3 = varset.VariantSet.from_variants([h1,h2])
        
        self.assertEqual(s3.diff(s2).nof_unit_vrt(),0)
        self.assertEqual(s2.diff(s3).nof_unit_vrt(),0)

    def test_VariantSet_cmp(self):
        vcf1 = pkg_file('genomvar.test','data/test_vcf.vcf')
        vcf2 = pkg_file('genomvar.test','data/test_vcf2.vcf')
        s1 = VariantSet.from_vcf(vcf1,parse_info=True,parse_samples=True)
        s2 = VariantSet.from_vcf(vcf2)
        diff = s1.diff(s2)
        self.assertEqual(diff.nof_unit_vrt(),14)
        # Now same diff but without loading in memory
        N = 0
        for vrt in s1.diff_vrt(s2).iter_vrt():
            N += vrt.nof_unit_vrt()
        self.assertEqual(N,14)
        
        comm = s1.comm(s2)
        self.assertEqual(len(list(comm.iter_vrt())),4)
        v1,v2 = sorted(comm.iter_vrt(),
                       key=lambda v: v.key)[:2]
        self.assertEqual(v1.attrib['info']['NSV'],1)
        self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(0,1))

        self.assertEqual(s2.comm(s1).nof_unit_vrt(), comm.nof_unit_vrt())

    def test_VariantSet_cmp2(self):
        vcf1 = pkg_file('genomvar.test','data/test_vcf5.vcf')
        vcf2 = pkg_file('genomvar.test','data/test_vcf6.vcf.gz')
        s1 = VariantSet.from_vcf(vcf1)
        s2 = VariantSet.from_vcf(vcf2)
        diff = s1.diff(s2)
        self.assertEqual(diff.nof_unit_vrt(),0)
        
        comm = s1.comm(s2)
        self.assertEqual(comm.nof_unit_vrt(),s1.nof_unit_vrt())
        self.assertEqual(s2.comm(s1).nof_unit_vrt(), comm.nof_unit_vrt())
        
class TestIO(TestCase):
    def test_from_vcf(self):
        vset = MutableVariantSet.from_vcf(test_vcf1,
                        reference=CHR15RGN,sample='SAMP1',normindel=True)

        vrt = list(vset.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(vrt),2)
        self.assertEqual(len({o.GT for o in vrt}),2)

        # Test presence of null operation
        vrt = list(vset.find_vrt('chr15rgn',20,25))
        self.assertEqual(len(vrt),2)
        self.assertEqual({o.GT for o in vrt},{(0,0,1),(0,1,0)})
        for _v in vrt:
            if _v.vtp!=variant.Null:
                self.assertEqual(_v.attrib['id'],'1')

    def test_from_variants_to_records(self):
        fac = variant.VariantFactory(reference=CHR15RGN,normindel=True)
        hap = Haplotype.from_variants([fac.from_edit('chr15rgn',1207,'G','C'),
                                      fac.from_edit('chr15rgn',1207,'G','T')])
        vs = VariantSet.from_variants(
            [fac.from_edit('chr15rgn',10043,'T','TCA'),
             fac.from_edit('chr15rgn',10045,'ACA','A'),
             hap])
        recs = vs.to_records()
        self.assertEqual(recs.shape,(4,))
        self.assertEqual(list(recs.dtype.fields),
                         ['chrom','start','end','ref','alt',
                          'vartype','phase_group'])
        
    def test_from_vcf_to_records(self):
        vs = VariantSet.from_vcf(test_vcf1,parse_info=True,
                                 parse_samples=True)

        # Test nested dtype
        recs = vs.to_records(nested=True)
        self.assertEqual(list(recs.dtype.fields),
                         ['chrom','start','end','ref','alt',
                          'vartype','phase_group',
                          'INFO','SAMPLES'])
        self.assertEqual(list(recs['INFO'].dtype.fields),
                         ['NSV','AF','DP4','ECNT','pl','mt','RECN','STR'])
        self.assertEqual(list(recs['SAMPLES'].dtype.fields),
                         ['SAMP1'])
        self.assertEqual(list(recs['SAMPLES']['SAMP1'].dtype.fields),
                         ['GT'])

        # Test not nested
        recs = vs.to_records(nested=False)
        self.assertEqual(list(recs.dtype.fields),
                         ['chrom','start','end','ref','alt',
                          'vartype','phase_group',
                          'INFO_NSV', 'INFO_AF', 'INFO_DP4', 'INFO_ECNT',
                          'INFO_pl', 'INFO_mt', 'INFO_RECN', 'INFO_STR',
                          'SAMPLES_SAMP1_GT'])
    # ΤΟDO
    # def test_reading_empty_ALT(self):
    #     fl = pkg_file('genomvar.test','data/test_vcf3.vcf')
    #     vset = from_vcf(fl,reference=CHR15RGN,parse_samples='SAMP1',
    #                                     normindel=True,
    #                                     allow_adjust_ploidy=True,
    #                                     debug=False,parse_null=False)

    #     vrt = list(vset.find_vrt('chr15rgn',20,30))
    #     self.assertEqual(len(vrt),1)
    #     v = vrt[0]
    #     self.assertEqual(type(v).__name__,'Null')

    def test_from_vcf_with_attr(self):
        s = MutableVariantSet.from_vcf(test_vcf1,parse_info=True,
                                       sample='SAMP1')
        _vrt = list(s.find_vrt('chr15rgn',150,160))
        self.assertEqual(len(_vrt),1)
        vrt = _vrt[0]
        self.assertEqual(vrt.attrib['info']['AF'],1.0)

        # Check multiallelic locus
        _vrt = list(s.find_vrt('chr15rgn',20,30))
        self.assertEqual(len(_vrt),2)
        for vrt in _vrt:
            if not vrt.is_instance(variant.Null):
                self.assertEqual(vrt.attrib['info']['AF'],0.5)

        # Check None/KeyError cases (".",field absent...)
        _vrt = list(filter(lambda o: not o.is_instance(variant.Null),
                           s.find_vrt('chr15rgn',450,460)))
        self.assertEqual(len(_vrt),1)
        vrt = _vrt[0]
        with self.assertRaises(KeyError):
            vrt.attrib['info']['Randomfields']

        _vrt = list(filter(lambda o: not o.is_instance(variant.Null),
                           s.find_vrt('chr15rgn',4750,4760)))
        self.assertEqual(len(_vrt),1)
        vrt = _vrt[0]
        self.assertEqual(vrt.attrib['info']['STR'],True)

    def test_from_vcf_no_sample_no_ref(self):
        _vcf = test_vcf1
        vset = MutableVariantSet.from_vcf(_vcf)
        v1,v2 = list(vset.find_vrt('chr15rgn',1200,1210))

        self.assertTrue(True if v1.GT in [(0,1),(1,0)] else False)
        self.assertTrue(True if v2.GT in [(0,1),(1,0)] else False)
        self.assertNotEqual(v1,v2)


    def test_from_vcf2(self):
        _vcf = pkg_file('genomvar.test','data/test_vcf7.vcf')
        vset = MutableVariantSet.from_vcf(_vcf)

    def test_wrong_sample(self):
        with self.assertRaises(VCFSampleMismatch):
            MutableVariantSet.from_vcf(test_vcf1,
                     reference=CHR15RGN,sample='ZAMP1',normindel=True)

    def test_to_vcf(self):
        s1 = MutableVariantSet(reference=CHR15RGN)
        s1.add_vrt(variant.Del("chr15rgn",23,24),GT=(1,0))
        s1.add_vrt(variant.SNP("chr15rgn",1206,"C"),GT=(1,0))
        s1.add_vrt(variant.MNP("chr15rgn",2093,"CCC"),GT=(1,0))

        # Ambigous indels
        vf = s1.get_factory(normindel=True)
        vb = vf.from_edit('chr15rgn',2343,'TTTCCA','TTCCATTCCA')
        s1.add_vrt(vb,GT=(0,1))
        vb = vf.from_edit('chr15rgn',9947,'TTT','T')
        s1.add_vrt(vb,GT=(0,1))

        stream = io.StringIO()
        s1.to_vcf(stream)
        stream.seek(0)

        rows = []
        for line in stream:
            if line.startswith('##'):
                continue
            else:
                vals = line.split('\t')
                self.assertEqual(len(vals),8)

                if vals[0]=='#CHROM':
                    continue

                rows.append(VCFRow(*vals))

        row0 = rows[0]
        self.assertEqual([row0.CHROM,row0.POS,row0.REF,row0.ALT],
                         ['chr15rgn',23,'AG','A'])
        row1 = rows[1]
        self.assertEqual([row1.CHROM,row1.POS,row1.REF,row1.ALT],
                         ['chr15rgn',1207,'G','C'])

        row3 = rows[3]
        self.assertEqual([row3.CHROM,row3.POS,row3.REF,row3.ALT],
                         ['chr15rgn',2344,'T','TTCCA'])

        row4 = rows[4]
        self.assertEqual([row4.CHROM,row4.POS,row4.REF,row4.ALT],
                         ['chr15rgn',9946,'CTT','C'])
        
    def test_from_to_vcf(self):
        variants1 = sorted(VCFReader(test_vcf1).iter_vrt(),
                           key=lambda v: v.key)
        vs = VariantSet.from_vcf(test_vcf1)
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')
        with open(tf.name,'wt') as fh:
            vs.to_vcf(fh)
        variants2 = sorted(VCFReader(tf.name).iter_vrt(),
                           key=lambda v: v.key)
        self.assertEqual(len(variants1),len(variants2))
        cnt = 0
        for v1,v2 in zip(variants1,variants2):
            self.assertTrue(v1.edit_equal(v2))

        
if __name__ == '__main__':
    unittest.main()
