from unittest import TestCase
import tempfile
import io
import operator
import itertools
import re
import warnings
from pkg_resources import resource_filename as pkg_file
import numpy as np
from genomvar.varset import VariantSetFromFile, VariantSetFromFile, VariantSet
from genomvar.variant import VariantBase,AmbigIndel,Haplotype,VariantFactory
from genomvar import varset, Reference
from genomvar import OverlappingHaplotypeVars,\
    UnsortedVariantFileError,VCFSampleMismatch,NoIndexFoundError
from genomvar.vcf import VCFRow,VCFReader, RESERVED_FORMAT
from genomvar import variant
from genomvar.test import MyTestCase

# Representation of an example used here a lot (example1.vcf)
# 
# REF      AG    T|   T|  C    G     TGG   TT    G|      T    T      CACAGTTCCAC
#          22    154  165 453  1206  2093  2099  3200    4754 6145   10044
# varset1  A     TT   TG  CT   C     CCC   GG    GG      TCG  C      T----------
#          AGG   TT   TG       T     CCC   GG    GG      T    G      G----------
#          AG                        --------            ------ phased
# FILT           h    h;n                  n     n

factory = variant.VariantFactory()

class TestVariantSetCase(MyTestCase):
    def test_empty_vcf(self):
        buf = io.StringIO()
        with open(pkg_file('genomvar.test','data/example1.vcf')) as fh:
            for line in itertools.takewhile(
                    lambda l: l.startswith('#'), fh):
                buf.write(line)
        buf.seek(0)
        vs = VariantSet.from_vcf(buf)
        self.assertEqual(vs.nof_unit_vrt(),0)

    def test_random_sample(self):
        vs = VariantSet.from_vcf(pkg_file('genomvar.test',
                                          'data/example1.vcf'))
        sample = vs.sample(5)
        self.assertEqual(len(sample), 5)
        
    
    def test_sort_chroms(self):
        vs = VariantSet.from_vcf(pkg_file('genomvar.test',
                                          'data/example2.vcf.gz'))
        vs.sort_chroms()
        self.assertEqual(list(vs.get_chroms()),['chr23','chr24'])

        vs.sort_chroms(key=lambda c: 1 if c=='chr24' else 2)
        self.assertEqual(list(vs.get_chroms()),['chr24','chr23'])

    def test_ovlp(self):
        #                        23
        # TTCACTTAGCATAATGTCTTCAAG|ATT
        #                         G
        #                          C <- not interfering

        ins = factory.from_edit(chrom='chr24',start=23,ref='G',
                                    alt='GG')
        s1 = VariantSet.from_variants([ins],)
        vb = factory.from_edit(chrom='chr24',start=24,ref='A',alt='C')
        self.assertEqual(len(s1.ovlp(vb)),0)
        vb2 = factory.from_edit(chrom='chr24',start=23,ref='G',alt='GG')
        self.assertEqual(len(s1.ovlp(vb2,match_ambig=False)),1)

    def test_from_variants(self):
        vfset = VariantSetFromFile(pkg_file('genomvar.test','data/example1.vcf'))
        vset = VariantSet.from_variants(list(vfset.iter_vrt()))
        vrt = list(vset.find_vrt('chr24',1200,1210))
        self.assertEqual(len(vrt),2)

        # Test error on out of reference bounds
        with self.assertRaises(ValueError):
            VariantSet.from_variants(
                list(vfset.iter_vrt())+[variant.SNP('chr24',10000000,'C')],
                reference=self.chr24)

        # Test error on chromosome not in reference
        with self.assertRaises(ValueError):
            vs = VariantSet.from_variants(
                list(vfset.iter_vrt())+[variant.SNP('chr2',10,'C')],
                reference=self.chr24)

    def test_from_variants_with_attributes(self):
        reader = VCFReader(
            pkg_file('genomvar.test','data/example1.vcf'))
        vset = VariantSet.from_variants(list(reader.iter_vrt(parse_info=True)))
        vrt = list(vset.find_vrt('chr24',1200,1210))
        self.assertEqual(len(vrt),2)

        v1 = vrt[0]
        self.assertEqual(v1.attrib['info']['NSV'], 1)
        self.assertEqual(v1.attrib['id'], '5')

        v2 = vrt[1]
        self.assertEqual(v2.attrib['id'], None)

        recs = vset.to_records()
        self.assertEqual(recs[0]['attrib']['info']['NSV'], 2)
        
    def test_from_vcf(self):
        vset = VariantSet.from_vcf(pkg_file('genomvar.test','data/example1.vcf'))

        self.assertEqual(vset.nof_unit_vrt(),18)
        self.assertEqual(list(vset.chroms),['chr24'])
        vrt = list(vset.find_vrt('chr24',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertTrue(v1.is_variant_instance(variant.SNP))
        self.assertEqual(v1.attrib['vcf_notation']['ref'],'G')
        self.assertEqual(v1.attrib['vcf_notation']['start'],1206)
        self.assertEqual(v1.attrib['vcf_notation']['row'],4)
        self.assertEqual(v1.attrib['id'],'5')
        self.assertEqual(v1.attrib['filter'],'PASS')
        self.assertEqual(v1.attrib['qual'],100)
        vrt = list(vset.find_vrt('chr24',154,156))
        self.assertEqual(len(vrt),1)
        v1 = vrt[0]
        self.assertTrue(v1.is_variant_instance(variant.Ins))
        vrt = list(vset.find_vrt('chr24',20,25))
        self.assertEqual(len(vrt),2)
        self.assertEqual(len(list(vset.iter_vrt())),16)

    def test_from_vcf_with_info(self):
        vset = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            parse_info=True)

        vrt = list(vset.find_vrt('chr24',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.alt,'C')
        self.assertEqual(v1.attrib['info']['NSV'],1)

    def test_from_vcf_with_samples(self):
        vset = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            parse_samples=True)

        vrt = list(vset.find_vrt('chr24',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(1,0))
        
    def test_asterisk_variant(self):
        vset = VariantSet.from_vcf(pkg_file('genomvar.test',
                            'data/example_with_asterisk.vcf.gz'),
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
        self.assertEqual(len(list(s1.find_vrt('chr24'))),4)


    def test_match(self):
        # REF      TGG   TT    
        #          2093  2099  
        # vs1      CCC   GG
        # vs2            CG    
        #          r1   r2,r3
        vs1 = VariantSet.from_vcf(pkg_file('genomvar.test','data/example1.vcf'))
        vrt_CG = factory.from_edit('chr24',2098,'TT','CG')
        vrt_CC = factory.from_edit('chr24',2098,'TT','CC')
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

    def test_match_with_haplotypes(self):
        # REF      TGG   TT    G|
        #          2093  2099  3200
        # varset1  CCC   GG
        #                CC    GG
        #          r1   r2,r3   r4
        r1 = factory.from_edit('chr24',2093,'TGG','CCC')
        r2 = factory.from_edit('chr24',2098,'TT','GG')
        r3 = factory.from_edit('chr24',2098,'TT','CC')
        r4 = factory.from_edit('chr24',3200,'GG','G')

        hap1 = Haplotype.from_variants([r1,r2])
        ccc1 = sorted(hap1.variants,key=lambda o: o.start)[0]
        hap2 = Haplotype.from_variants([r3,r4])
        s1 = VariantSet.from_variants([hap1, hap2])

        hap = Haplotype.from_variants([r1,r4])
        ccc2 = sorted(hap.variants,key=lambda o: o.start)[0]

        match = s1.match(hap)
        self.assertEqual(len(match),1)
        d1 = {k:[v2.base for v2 in v] for k,v in match.items()}
        d2 = {ccc2.key:[ccc1]}
        self.assertEqual(len(d1), len(d2))
        self.assertEqual(list(d1.keys()), list(d2.keys()))
        for k in d1:
            self.assertTrue(
                all(
                    [v1.edit_equal(v2) for v1,v2 in zip(d1[k], d2[k])]))

class TestVariantFileSetWithIndexCase(MyTestCase):
    def test_complex_info_example(self):
        vset = VariantSetFromFile(
            pkg_file('genomvar.test','data/example_gnomad_1.vcf.gz'),
            parse_info=True,
            index=True)
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
        ivfs = VariantSetFromFile(
            pkg_file('genomvar.test','data/example2.vcf.gz'),
            index=True)
        vs = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example2.vcf.gz'))

        self.assertEqual(
            sum([v.nof_unit_vrt() for v in ivfs.find_vrt('chr24')]),
            sum([v.nof_unit_vrt() for v in vs.find_vrt('chr24')]))

        self.assertEqual(
            sum([v.nof_unit_vrt() for v in ivfs.iter_vrt()]),
            sum([v.nof_unit_vrt() for v in vs.iter_vrt()]))
        
    def test_find_vrt2(self):
        vset = VariantSetFromFile(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            reference=self.chr24, index=True)
        self.assertEqual(len(list(vset.find_vrt(rgn='chr24:1200-1210'))),2)
        v1,v2 = list(vset.find_vrt(rgn='chr24:1200-1210'))
        self.assertEqual([v1.start,v1.end],[1206,1207])
        self.assertEqual([v2.start,v2.end],[1206,1207])

        self.assertEqual(len(list(vset.find_vrt(rgn='chr24:3200-3205'))),1)
        v1 = list(vset.find_vrt(rgn='chr24:3200-3205'))[0]
        self.assertEqual([v1.start,v1.end],[3201,3202])
        
        self.assertEqual(len(list(vset.find_vrt(rgn='chr24:20-30'))),2)
        v1,v2 = list(vset.find_vrt(rgn='chr24:20-30'))
        self.assertEqual([v1.start,v1.end,type(v1.base)],
                         [23,24,variant.Del])
        self.assertEqual([v2.start,v2.end,type(v2.base)],
                         [24,25,variant.Ins])

    def test_find_nonexistent_chrom(self):
        vcf  = pkg_file('genomvar.test','data/example_1000genomes_1.vcf.gz')
        vset = VariantSetFromFile(vcf, index=True)
        self.assertEqual(list(vset.find_vrt('chr24')),[])
        
    def test_match(self):
        # REF      TGG   TT    
        #          2093  2099  
        # vs1      CCC   GG
        # vrt            CG    
        #          r1   r2,r3
        vs1 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            index=True)
        vrt = factory.from_edit('chr24',2098,'TT','CG')
        self.assertEqual(len(vs1.match(vrt)),1)

        # Test insertion
        vrt = factory.from_edit('chr24',22,'AG','AGG')
        match = vs1.match(vrt)
        self.assertEqual(len(match),1)

    def test_index_errors(self):
        file = pkg_file('genomvar.test','data/example1.vcf')
        with self.assertRaises(NoIndexFoundError):
            vset = VariantSetFromFile(
                file, reference=self.chr24, index=True)

        vset = VariantSetFromFile(
            file, reference=self.chr24)

        with self.assertRaises(ValueError) as cm:
            list(vset.find_vrt('chr1', 1, 100))
        error = cm.exception
        self.assertIn('index is required', error.args[0].lower())

    def test_wrong_chrom_name_in_ref(self):
        ref = Reference(pkg_file(__name__,'data/chr25.fasta'))
        vset = VariantSetFromFile(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            reference=ref, index=True)
        self.assertEqual(len(list(vset.find_vrt(rgn='chr24:1200-1210'))),2)
        ref.close()
        

    def test_class(self):
        vset = VariantSetFromFile(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            parse_info=True,
            reference=self.chr24,
            parse_samples='SAMP1')

        # Test find_vrt and returned INFO
        vrt = list(vset.find_vrt('chr24',1200,1210))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.attrib['info']['NSV'],1)
        self.assertEqual(v2.attrib['info']['RECN'],19)

        # Test multiallelic
        vrt = list(vset.find_vrt('chr24',20,30))
        self.assertEqual(len(vrt),2)
        v1,v2 = vrt
        self.assertEqual(v1.attrib['info']['AF'],0.5)
        self.assertEqual(v2.attrib['info']['AF'],0.5)

        # Test find_vrt precision
        vrt = list(vset.find_vrt('chr24',2095,2096))
        self.assertEqual(len(vrt),1)
        vrt = list(vset.find_vrt('chr24',2098,2100))
        self.assertEqual(len(vrt),1)

        # Test find all variants
        self.assertEqual(len(list(vset.find_vrt())),16)

        # Test finding all variants
        self.assertEqual(len(list(vset.find_vrt())),16)
        
    def test_no_index(self):
        with self.assertRaises(NoIndexFoundError):
            vset = VariantSetFromFile(
                pkg_file('genomvar.test','data/example3.vcf'),
                index=True)

    def test_ctg_len_without_ref(self):
        vset = VariantSetFromFile(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            parse_samples='SAMP1', index=True)
        self.assertEqual(vset.chroms,{'chr24'})
        
class TestVariantFileCase(MyTestCase):
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

        vs1 = VariantSetFromFile(pkg_file('genomvar.test','data/example1.vcf'))
        vs2 = VariantSetFromFile(tf.name)
        with self.assertRaises(UnsortedVariantFileError):
            list(vs1.diff_vrt(vs2).iter_vrt())

class TestIO(MyTestCase):
    def test_from_vcf(self):
        vset = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            reference=self.chr24, normindel=True)

        vrt = list(vset.find_vrt('chr24',1200,1210))
        self.assertEqual(len(vrt),2)

        # Test presence of null operation
        vrt = list(vset.find_vrt('chr24',20,25))
        self.assertEqual(len(vrt),2)
        for _v in vrt:
            if not _v.is_variant_instance(variant.Null):
                self.assertEqual(_v.attrib['id'],'1')

    def test_from_vcf_problematic(self):
        vset = VariantSet.from_vcf(
            pkg_file('genomvar.test', 'data/example5.vcf.gz'),
            parse_info=True)

        v = list(vset.find_vrt('1', 160109406, 160109409))
        
        self.assertEqual(len(v), 1)
        v = v[0]
        self.assertEqual(v.attrib['info']['PH'], ('.',))

    def test_from_bcf(self):
        vset = VariantSet.from_vcf(
            pkg_file('genomvar.test', 'data/example3.bcf'))
        self.assertEqual(len(list(vset.iter_vrt())), 8)

    def test_from_variants_to_records(self):
        fac = variant.VariantFactory(reference=self.chr24,normindel=True)
        hap = Haplotype.from_variants([fac.from_edit('chr24',1207,'G','C'),
                                      fac.from_edit('chr24',1207,'G','T')])
        vs = VariantSet.from_variants(
            [fac.from_edit('chr24',10043,'T','TCA'),
             fac.from_edit('chr24',10045,'ACA','A'),
             hap])
        recs = vs.to_records()
        self.assertEqual(recs.shape, (4,))
        self.assertEqual(list(recs.dtype.fields),
                         ['chrom','start','end','ref','alt',
                          'vartype','phase_group','attrib'])
        
    def test_from_vcf_to_records(self):
        vs = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            parse_info=True,parse_samples=True)

        self.assertEqual(vs._samples, ['SAMP1'])

        # Test nested dtype
        recs = vs.to_records(nested=True)
        self.assertEqual(list(recs.dtype.fields),
                         ['chrom','start','end','ref','alt',
                          'vartype','phase_group',
                          'info','SAMPLES'])
        self.assertEqual(list(recs['info'].dtype.fields),
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
                          'info_NSV', 'info_AF', 'info_DP4', 'info_ECNT',
                          'info_pl', 'info_mt', 'info_RECN', 'info_STR',
                          'SAMPLES_SAMP1_GT'])
    # ΤΟDO
    # def test_reading_empty_ALT(self):
    # ....
    #     self.assertEqual(type(v).__name__,'Null')

    def test_from_vcf_with_attr(self):
        s = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'), parse_info=True)
        _vrt = list(s.find_vrt('chr24',150,160))
        self.assertEqual(len(_vrt),1)
        vrt = _vrt[0]
        self.assertEqual(vrt.attrib['info']['AF'],1.0)

        # Check multiallelic locus
        _vrt = list(s.find_vrt('chr24',20,30))
        self.assertEqual(len(_vrt),2)
        for vrt in _vrt:
            if not vrt.is_variant_instance(variant.Null):
                self.assertEqual(vrt.attrib['info']['AF'],0.5)

        # Check None/KeyError cases (".",field absent...)
        _vrt = list(filter(lambda o: not o.is_variant_instance(variant.Null),
                           s.find_vrt('chr24',450,460)))
        self.assertEqual(len(_vrt),1)
        vrt = _vrt[0]
        with self.assertRaises(ValueError):
            vrt.attrib['info']['Randomfields']

        _vrt = list(filter(lambda o: not o.is_variant_instance(variant.Null),
                           s.find_vrt('chr24',4750,4760)))
        self.assertEqual(len(_vrt),1)
        vrt = _vrt[0]
        self.assertEqual(vrt.attrib['info']['STR'],True)

    def test_to_vcf_no_ref(self):
        vs1 = VariantSet.from_variants(
            [variant.Del("chr24",23,24),
             variant.SNP("chr24",1206,"C")]
        )

        buf = io.StringIO()
        vs1.to_vcf(buf, reference=self.chr24)
        buf.seek(0)

        vs2 = VariantSet.from_vcf(buf)
        self.assertEqual(vs1.comm(vs2).nof_unit_vrt(), 2)

    def test_from_to_vcf(self):
        fl = pkg_file('genomvar.test','data/example1.vcf')
        variants1 = sorted(
            VCFReader(fl).iter_vrt(parse_info=True),
            key=lambda v: v.key)
        vs = VariantSet.from_vcf(fl, parse_info=True)
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')
        with open(tf.name,'wt') as fh:
            vs.to_vcf(fh)

            
        with open(tf.name,'rt') as fh:
            fh.seek(0)
            self.assertIn(
                '##INFO=<ID=DP4,Number=4,Type=Integer,Description="Test for multinumber field">',
                fh.read().splitlines())
        variants2 = sorted(VCFReader(tf.name).iter_vrt(parse_info=True),
                           key=lambda v: v.key)
        self.assertEqual(len(variants1),len(variants2))
        cnt = 0
        for v1,v2 in zip(variants1,variants2):
            self.assertTrue(v1.edit_equal(v2))
            self.assertEqual(v1.attrib['info']['NSV'], v2.attrib['info']['NSV'])

        
    def test_from_variants_vcf(self):
        vs0 = varset.VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            parse_info=True)
        variants1 = sorted(
            vs0.iter_vrt(),
            key=lambda v: v.key)
        vs = VariantSet.from_variants(variants1)
        _desc = 'Test for multinumber field'
        info_spec_tuples = [('DP4', 4, 'Integer', _desc),
                            ('NSV', 1, 'Integer')]
        info_spec_dict = vs0.dtype['info']
        for info_spec in (info_spec_tuples, info_spec_dict):
            tf = tempfile.NamedTemporaryFile(suffix='.vcf')
            with open(tf.name,'wt') as fh:
                vs.to_vcf(fh, info_spec=info_spec)

            with open(tf.name,'rt') as fh:
                self.assertIn(
                    '##INFO=<ID=DP4,Number=4,Type=Integer,Description="{}">'\
                          .format(_desc),
                    fh.read().splitlines())
                fh.seek(0)
                # print(fh.read())
            variants2 = sorted(VCFReader(tf.name).iter_vrt(parse_info=True),
                               key=lambda v: v.key)
            self.assertEqual(len(variants1),len(variants2))
            cnt = 0
            for v1,v2 in zip(variants1,variants2):
                self.assertTrue(v1.edit_equal(v2))
                self.assertEqual(v1.attrib['info']['NSV'], v2.attrib['info']['NSV'])

    def test_from_variants_to_vcf_with_info(self):
        variants1 = sorted(VCFReader(
            pkg_file('genomvar.test','data/example1.vcf')).iter_vrt(parse_info=True),
                        key=lambda v: v.key)
        vs = VariantSet.from_variants(variants1)
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')

        # Test Invalid specs
        invalid_specs = [('NSV',),
                         ('NSV', 1, 'Integedr'),
                         ('NSV', 'C', 'Integer', 'Number of Simple Variants')]
        _buf = io.StringIO()
        for spec in invalid_specs:
            with self.assertRaises(ValueError) as cm:
                vs.to_vcf(_buf, info_spec=[spec])
            exc = cm.exception
            self.assertTrue('INFO spec' in exc.args[0])
        
        with open(tf.name,'wt') as fh:
            vs.to_vcf(fh, info_spec=[('NSV', 1, 'Integer', 'Number of Simple Variants'),
                                     ('AF', 'A', 'Float', '', 'source', 'version')])
            
        with open(tf.name,'rt') as fh:
            self.assertIn(
                '##INFO=<ID=NSV,Number=1,Type=Integer,'\
                +'Description="Number of Simple Variants">',
                fh.read().splitlines())
        variants2 = sorted(VCFReader(tf.name).iter_vrt(parse_info=True),
                           key=lambda v: v.key)
        self.assertEqual(len(variants1),len(variants2))
        cnt = 0
        for v1,v2 in zip(variants1,variants2):
            self.assertTrue(v1.edit_equal(v2))
            self.assertEqual(v1.attrib['info']['NSV'], v2.attrib['info']['NSV'])

    def test_from_variants_to_vcf_with_sampdata(self):
        file = pkg_file('genomvar.test', 'data/example3.vcf')
        variants1 = sorted(
            VCFReader(file).iter_vrt(parse_samples=True),
            key=lambda v: v.key)
        vs = VariantSet.from_variants(variants1)
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')

        with open(tf.name,'wt') as fh:
            vs.to_vcf(fh,
                      format_spec=[RESERVED_FORMAT.GT,
                                   ('AD', 'R', 'Integer', '')],
                      samples=['SAMP1'])
            
        with open(tf.name,'rt') as fh:
            fh.seek(0)
            self.assertIn(
                '##FORMAT=<ID=AD,Number=R,Type=Integer,'\
                +'Description="">',
                fh.read().splitlines())
        variants2 = sorted(VCFReader(tf.name).iter_vrt(parse_samples=True),
                           key=lambda v: v.key)
        self.assertEqual(len(variants1),len(variants2))
        cnt = 0
        for v1, v2 in zip(variants1, variants2):
            self.assertTrue(v1.edit_equal(v2))
            self.assertEqual(v1.attrib['samples']['SAMP1']['AD'],
                             v2.attrib['samples']['SAMP1']['AD'])

    def test_from_string_buffer(self):
        buf = io.StringIO() 
        with open(pkg_file('genomvar.test','data/example1.vcf'),'rt') as fh:
            for line in fh:
                buf.write(line)
        buf.seek(0)
        vs = VariantSet.from_vcf(buf)
        self.assertEqual(len(list(vs.find_vrt('chr24',150,160))), 1)
        self.assertEqual(len(list(vs.find_vrt('chr24',20,30))), 2)

    def test_sv_types(self):
        with warnings.catch_warnings(record=True) as wrn:
            vs = VariantSet.from_vcf(pkg_file('genomvar.test', 'data/example4.vcf.gz'))
            warnings.simplefilter('always')
            self.assertEqual(vs.nof_unit_vrt(), 100)
            self.assertGreater(len(wrn), 1)
            self.assertIn('Structural', str(wrn[-1].message))

if __name__ == '__main__':
    unittest.main()
