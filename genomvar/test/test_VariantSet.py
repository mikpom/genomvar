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
    UnsortedVariantFileError,VCFSampleMismatch
from genomvar.vcf import VCFRow,VCFReader,header_simple
from genomvar import variant

# The following is a scheme for example used here a lot (example1.vcf)
# 
# REF      AG    T|   T|  C    G     TGG   TT    G|      T    T      CACAGTTCCAC
#          22    154  165 453  1206  2093  2099  3200    4754 6145   10044
# varset1  A     TT   TG  CT   C     CCC   GG    GG      TCG  C      T----------
#          AGG   TT   TG       T     CCC   GG    GG      T    G      G----------
#          AG                        --------            ------ phased
# FILT           h    h;n                  n     n

factory = variant.VariantFactory()
CHR15RGN = pkg_file('genomvar.test','data/chr15rgn.fna')

class TestVariantSetCase(TestCase):
    def test_empty_vcf(self):
        buf = io.StringIO()
        with open(pkg_file('genomvar.test','data/example1.vcf')) as fh:
            for line in itertools.takewhile(
                    lambda l: l.startswith('#'), fh):
                buf.write(line)
        buf.seek(0)
        vs = VariantSet.from_vcf(buf)
        self.assertEqual(vs.nof_unit_vrt(),0)
    
    def test_sort_chroms(self):
        vs = VariantSet.from_vcf(pkg_file('genomvar.test',
                                          'data/example2.vcf.gz'))
        vs.sort_chroms()
        self.assertEqual(list(vs.get_chroms()),['chr11rgn','chr15rgn'])

        vs.sort_chroms(key=lambda c: 1 if c=='chr15rgn' else 2)
        self.assertEqual(list(vs.get_chroms()),['chr15rgn','chr11rgn'])

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
        vfset = VariantFileSet(pkg_file('genomvar.test','data/example1.vcf'))
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
        vset = VariantSet.from_vcf(pkg_file('genomvar.test','data/example1.vcf'))

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
        self.assertEqual(v1.attrib['vcf_notation']['row'],4)
        vrt = list(vset.find_vrt('chr15rgn',154,156))
        self.assertEqual(len(vrt),1)
        v1 = vrt[0]
        self.assertTrue(v1.is_instance(variant.Ins))
        vrt = list(vset.find_vrt('chr15rgn',20,25))
        self.assertEqual(len(vrt),2)
        self.assertEqual(len(list(vset.iter_vrt())),16)

    def test_from_vcf_with_info(self):
        vset = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            parse_info=True)

        vrt = list(vset.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.alt,'C')
        self.assertEqual(v1.attrib['info']['NSV'],1)

    def test_from_vcf_with_samples(self):
        vset = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            parse_samples=True)

        vrt = list(vset.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(1,0))
        
    def test_asterisk_variant(self):
        vset = VariantSet.from_vcf(pkg_file('genomvar.test',
                            'data/example_with_asterisk.vcf'),
                                   parse_info=True)

        vrt = list(vset.find_vrt('chr1',995507,995515))
        self.assertEqual(len(vrt),3)
        # v1,v2 = vrt
        # self.assertEqual(v1.alt,'TTTTT')
        # self.assertEqual(v2.alt,'C')

    def test_find_vrt_chrom_only(self):
        s1 = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example2.vcf.gz'),
            parse_info=True,parse_samples=True)
        self.assertEqual(len(list(s1.find_vrt('chr15rgn'))),4)


    def test_match(self):
        # REF      TGG   TT    
        #          2093  2099  
        # vs1      CCC   GG
        # vs2            CG    
        #          r1   r2,r3
        vs1 = VariantSet.from_vcf(pkg_file('genomvar.test','data/example1.vcf'))
        vrt_CG = factory.from_edit('chr15rgn',2098,'TT','CG')
        vrt_CC = factory.from_edit('chr15rgn',2098,'TT','CC')
        self.assertEqual(len(vs1.match(vrt_CG)),1)
        self.assertEqual(len(vs1.match(vrt_CG,match_partial=False)),0)
        self.assertEqual(len(vs1.match(vrt_CC)),0)
        self.assertEqual(len(vs1.match(vrt_CC,match_partial=False)),0)

    def test_match2(self):
        #       1
        # R     TCACAG
        # del1  T--CAG
        # del2  TCA--G
        del1 = variant.Del('chrom',1,3)
        adel1 = variant.AmbigDel('chrom',(1,1),(5,3),'CA','')
        del2 = variant.Del('chrom',2,4)
        adel2 = variant.AmbigDel('chrom',(1,2),(5,4),'AC','')
        

        vs = VariantSet.from_variants([del1])
        self.assertEqual(len(vs.match(adel1)),1)
        self.assertEqual(len(vs.match(adel2)),0)
        self.assertEqual(len(vs.match(adel2,match_ambig=True)),1)

        vs2 = VariantSet.from_variants([adel1])
        self.assertEqual(len(vs2.match(adel1)),1)
        self.assertEqual(len(vs2.match(del1)),1)
        self.assertEqual(len(vs2.match(adel2)),0)
        self.assertEqual(len(vs2.match(adel2,match_ambig=True)),1)
        self.assertEqual(len(vs2.match(del2,match_ambig=True)),1)
        self.assertEqual(len(vs2.match(del2,match_ambig=False)),0)

    def test_drop_duplicates(self):
        buf = io.StringIO()
        with open(pkg_file('genomvar.test','data/example1.vcf')) as fh:
            for line in fh:
                if line.startswith('#'):
                    buf.write(line)
                else:
                    for i in range(2):
                        buf.write(line)
        buf.seek(0)
        vs = VariantSet.from_vcf(buf,parse_info=True,parse_samples=True)
        vs2,dropped = vs.drop_duplicates(return_dropped=True)
        self.assertEqual(vs.nof_unit_vrt()/2,vs2.nof_unit_vrt())
        self.assertEqual(vs2.nof_unit_vrt(),
                         sum([v.nof_unit_vrt() for v in dropped]))

class TestIndexedVariantFileCase(TestCase):
    def test_complex_INFO_example(self):
        vset = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example_gnomad_1.vcf.gz'),
            parse_info=True)
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
            pkg_file('genomvar.test','data/example2.vcf.gz'))
        vs = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example2.vcf.gz'))

        self.assertEqual(
            sum([v.nof_unit_vrt() for v in ivfs.find_vrt('chr15rgn')]),
            sum([v.nof_unit_vrt() for v in vs.find_vrt('chr15rgn')]))

        self.assertEqual(
            sum([v.nof_unit_vrt() for v in ivfs.iter_vrt()]),
            sum([v.nof_unit_vrt() for v in vs.iter_vrt()]))
        
    def test_find_vrt2(self):
        vset = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            reference=CHR15RGN)
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

    def test_find_nonexistent_chrom(self):
        vcf  = pkg_file('genomvar.test','data/example_1000genomes_1.vcf.gz')
        vset = IndexedVariantFileSet(vcf)
        self.assertEqual(list(vset.find_vrt('chrZ')),[])
        
    def test_match(self):
        # REF      TGG   TT    
        #          2093  2099  
        # vs1      CCC   GG
        # vrt            CG    
        #          r1   r2,r3
        for ref in CHR15RGN,None: # should work with or without ref
            vs1 = IndexedVariantFileSet(
                pkg_file('genomvar.test','data/example1.vcf.gz'))
            vrt = factory.from_edit('chr15rgn',2098,'TT','CG')
            self.assertEqual(len(vs1.match(vrt)),1)

            # Test insertion
            vrt = factory.from_edit('chr15rgn',22,'AG','AGG')
            match = vs1.match(vrt)
            self.assertEqual(len(match),1)

    def test_error_on_format(self):
        with self.assertRaises(OSError):
            vset = IndexedVariantFileSet(
                pkg_file('genomvar.test','data/example1.vcf'),
                reference=CHR15RGN)

    def test_wrong_chrom_name_in_ref(self):
        vset = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            reference=pkg_file('genomvar.test','data/chr1rgn.fasta'))
        self.assertEqual(len(list(vset.find_vrt(rgn='chr15rgn:1200-1210'))),2)

    def test_class(self):
        vset = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            parse_info=True,reference=CHR15RGN,parse_samples='SAMP1')

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
        vset = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            parse_samples='SAMP1')
        self.assertEqual(vset.chroms,{'chr15rgn'})
        
class TestVariantFileCase(TestCase):
    def test_unsorted_VCF_input(self):
        header = []
        lines = []
        with open(pkg_file('genomvar.test','data/example1.vcf'),'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    header.append(line)
                else:
                    lines.append(line)
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')
        with open(tf.name,'wt') as fh:
            fh.writelines(header)
            fh.writelines(reversed(lines))

        vs1 = VariantFileSet(pkg_file('genomvar.test','data/example1.vcf'))
        vs2 = VariantFileSet(tf.name)
        with self.assertRaises(UnsortedVariantFileError):
            list(vs1.diff_vrt(vs2).iter_vrt())

class MutableVariantSetTestCase(TestCase):
    def test_sort_chroms(self):
        vs = MutableVariantSet.from_vcf(
            pkg_file('genomvar.test','data/example2.vcf.gz'),
            max_ploidy=4)
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

    def test_iter_vrt_by_chrom(self):
        vset = MutableVariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
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
        kwargs = {'vcf':pkg_file('genomvar.test','data/example1.vcf'),
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

        
class TestIO(TestCase):
    def test_from_vcf(self):
        vset = MutableVariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
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
        vs = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            parse_info=True,parse_samples=True)

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
    # ....
    #     self.assertEqual(type(v).__name__,'Null')

    def test_from_vcf_with_attr(self):
        s = MutableVariantSet.from_vcf(pkg_file('genomvar.test','data/example1.vcf'),parse_info=True,
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
        _vcf = pkg_file('genomvar.test','data/example1.vcf')
        vset = MutableVariantSet.from_vcf(_vcf)
        v1,v2 = list(vset.find_vrt('chr15rgn',1200,1210))

        self.assertTrue(True if v1.GT in [(0,1),(1,0)] else False)
        self.assertTrue(True if v2.GT in [(0,1),(1,0)] else False)
        self.assertNotEqual(v1,v2)

    def test_wrong_sample(self):
        with self.assertRaises(VCFSampleMismatch):
            MutableVariantSet.from_vcf(pkg_file('genomvar.test','data/example1.vcf'),
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
        variants1 = sorted(VCFReader(
            pkg_file('genomvar.test','data/example1.vcf')).iter_vrt(),
                        key=lambda v: v.key)
        vs = VariantSet.from_vcf(pkg_file('genomvar.test','data/example1.vcf'))
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')
        with open(tf.name,'wt') as fh:
            vs.to_vcf(fh)
        variants2 = sorted(VCFReader(tf.name).iter_vrt(),
                           key=lambda v: v.key)
        self.assertEqual(len(variants1),len(variants2))
        cnt = 0
        for v1,v2 in zip(variants1,variants2):
            self.assertTrue(v1.edit_equal(v2))

    def test_from_string_buffer(self):
        buf = io.StringIO() 
        with open(pkg_file('genomvar.test','data/example1.vcf'),'rt') as fh:
            for line in fh:
                buf.write(line)
        buf.seek(0)
        vs = VariantSet.from_vcf(buf)
        self.assertEqual(len(list(vs.find_vrt('chr15rgn',150,160))), 1)
if __name__ == '__main__':
    unittest.main()
