import operator
from unittest import TestCase
from pkg_resources import resource_filename as pkg_file
import itertools
import tempfile
from genomvar import variant
from genomvar.variant import VariantFactory
from genomvar.varset import VariantSet,\
    VariantSetFromFile,VariantSetFromFile
from genomvar import DifferentlySortedChromsError, NoIndexFoundError
from genomvar.test import MyTestCase

factory = variant.VariantFactory()

class TestVariantSetComparison(MyTestCase):
    def test_many_chroms_shuffled(self):
        vs1 = VariantSet.from_vcf(
            pkg_file('genomvar.test','data/example3.vcf'))
        vrt = list(vs1.iter_vrt())
        vrt2 = [vrt[i] for i in (0,7,1,6,2,5,3,4)]
        vs2 = VariantSet.from_variants(vrt2)

        self.assertEqual(vs1.comm(vs2).nof_unit_vrt(),
                         vs1.nof_unit_vrt())
        self.assertEqual(vs1.diff(vs2).nof_unit_vrt(), 0)

        self.assertEqual(vs2.comm(vs1).nof_unit_vrt(),
                         vs2.nof_unit_vrt())
        self.assertEqual(vs2.diff(vs1).nof_unit_vrt(), 0)
        
    def test_diff_del(self):
        vrt = variant.Del(chrom="chr1",start=6751613,end=6751627)
        vs1 = VariantSet.from_variants([vrt])
        vs2 = VariantSet.from_variants([vrt])
        self.assertEqual(len(list(vs1.diff(vs2).iter_vrt())),0)
        self.assertEqual(len(list(vs1.comm(vs2).iter_vrt())),1)

    def test_diff_ins(self):
        vrt = variant.Ins(chrom="chr1",start=6751613,alt='AGTC')
        vs1 = VariantSet.from_variants([vrt])
        vs2 = VariantSet.from_variants([vrt])
        self.assertEqual(len(list(vs1.diff(vs2).iter_vrt())),0)
        self.assertEqual(len(list(vs1.comm(vs2).iter_vrt())),1)

    def test_strip_order_dependent_Ambig(self):
        #    10043
        # R  T--CA--CAG
        # v1 TCACA--CAG
        # v2 T--CACACAG
        factory = variant.VariantFactory(
            reference=pkg_file('genomvar.test','data/chr24.fna'),
            normindel=True)
        v1 = factory.from_edit('chr24',10043,'T','TCA')
        v2 = factory.from_edit('chr24',10045,'A','ACA')
        s1 = VariantSet.from_variants([v1])
        s2 = VariantSet.from_variants([v2])

        diff = s1.diff(s2,match_ambig=True)
        self.assertEqual(len(list(diff.iter_vrt())),0)
        diff = s1.diff(s2,match_ambig=False)
        self.assertEqual(len(list(diff.iter_vrt())),1)
    
class TestSetComparisonCase(MyTestCase):

    def test_diff_vrt(self):
        # REF    G     TTGG         C
        #        1206   2093        10044
        # s1     C      CCC         T
        #        G      CCC         G
        #                           C
        #              2092         10044
        # s2           TT           T
        #              AC           T

        variants1 = [
            ('chr24',1206,'G','C'),
            ('chr24',2093,'TG','CC'),
            ('chr24',2095,'G','C'),
            ('chr24',10044,'C','T'),
            ('chr24',10044,'C','G')
        ]
        variants2 = [
            ('chr24',2092,'TT','AC'),
            ('chr24',10044,'C','T')
        ]

        s1 = VariantSet.from_variants(
            [factory.from_edit(*v) for v  in variants1])
        s2 = VariantSet.from_variants(
            [factory.from_edit(*v) for v  in variants2])

        diff = s1.diff(s2)
        self.assertEqual(diff.nof_unit_vrt(),4)
        
        v1,v2,v3,v4 = list(diff.iter_vrt())

        self.assertEqual([v1.start,v1.ref,v1.alt],
                         [1206,'G','C'])
        self.assertEqual([v4.start,v4.ref,v4.alt],
                         [10044,'C','G'])

        
    def test_diff_of_ambig_indel(self):
        #     9945
        #     CTTTTTCAT
        # s1  CTT--TCAT
        # s2  C--TTTCAT
        s1 = VariantSet.from_variants(
            [factory.from_edit('chr24',9947,'TT','T')])
        s2 = VariantSet.from_variants(
            [factory.from_edit('chr24',9945,'CTT','C')])

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

        vfac = VariantFactory(reference=self.chr24, normindel=True)
        variants1 = [
            ('chr24',1206,'G','C'),
            ('chr24',2093,'TG','CC'),
            ('chr24',2095,'G','C'),
            ('chr24',10044,'C','T'),
            ('chr24',10044,'C','G'),
            ('chr24',10051,'C','CCC')
        ]
        variants2 = [
            ('chr24',2092,'TT','AC'),
            ('chr24',10044,'C','T'),
            ('chr24',10052,'C','CCC'),
            ('chr24',10053,'A','AG')
        ]

        s1 = VariantSet.from_variants(
            [vfac.from_edit(*v) for v  in variants1])
        s2 = VariantSet.from_variants(
            [vfac.from_edit(*v) for v  in variants2])
        
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

        self.assertEqual(nof['com1'],3)
        self.assertEqual(nof['com2'],nof['com1'])
        self.assertEqual(nof['com1']+nof['s1_s2'],nof['s1'])
        self.assertEqual(nof['com1']+nof['s2_s1'],nof['s2'])

    def test_variant_cluster(self):
        vs1 = VariantSet.from_variants(
            [variant.SNP('chr1',13366968,'A'),
             variant.SNP('chr1',13366969,'T'),
             variant.SNP('chr1',13366970,'G')])

        vs2 = VariantSet.from_variants(
            [variant.Del('chr1',13366967,13366969),
             variant.Ins('chr1',13366971,'TG')])

        diff = list(vs1.diff(vs2).iter_vrt())
        self.assertEqual(len(diff), 3)

        vs3 = VariantSet.from_variants(
            [variant.MNP("chr1",13366968,'ATG')])
        self.assertEqual(len(list(vs1.diff(vs3,match_partial=True)\
                                     .iter_vrt())), 0)
        self.assertEqual(len(list(vs3.diff(vs1,match_partial=True)\
                                     .iter_vrt())), 0)
        self.assertEqual(len(list(vs1.diff(vs3,match_partial=False)\
                                     .iter_vrt())), 3)
        self.assertEqual(len(list(vs3.diff(vs1,match_partial=False)\
                                     .iter_vrt())), 1)
        
    def test_variant_cluster2(self):
        vs1 = VariantSet.from_variants(
            [variant.SNP('chr1',13366968,'A'),
             variant.SNP('chr1',13366969,'T'),
             variant.SNP('chr1',13366970,'G')])

        vs2 = VariantSet.from_variants(
            [variant.SNP('chr1',13366968,'A'),
             variant.SNP('chr1',13366969,'T'),
             variant.SNP('chr1',13366970,'G'),
             variant.Ins('chr1',13366970,'TTT')])

        com = list(vs1.comm(vs2, match_partial=False).iter_vrt())
        self.assertEqual(len(com), 3)

    def test_mnp_com_split(self):
        #                           23
        #    TTCACTTAGCATAATGTCTTCAAGATT
        # v1                       TT -single
        # v2                       TT -splitted
        vset1 = VariantSet.from_variants(
            [factory.from_edit('chr24',22,'AG','TT')])
        vset2 = VariantSet.from_variants(
            [factory.from_edit('chr24',22,'A','T'),
             factory.from_edit('chr24',23,'G','T')])

        com = vset1.comm(vset2)
        v = list(com.iter_vrt())[0]

        self.assertEqual([v.start,v.ref],[22,'AG'])
        self.assertFalse(list(vset1.diff(vset2).iter_vrt()))


    def test_in_the_middle(self):
        #      2
        #    TTCACTTAGCAT
        # v1   GGG
        # v2    G
        s1 = VariantSet.from_variants(
            [factory.from_edit('chr24',2,'CAC','GGG')])
        s2 = VariantSet.from_variants(
            [factory.from_edit('chr24',3,'A','G')])

        diff = s1.diff(s2)
        vv = list(diff.iter_vrt())
        self.assertEqual(len(vv),2)
        v1,v2 = vv
        self.assertEqual([v1.start,v1.ref,v1.alt],[2,'C','G'])
        self.assertEqual([v2.start,v2.ref,v2.alt],[4,'C','G'])

        self.assertEqual(len(list(s1.diff(s2,match_partial=False)
                                  .iter_vrt())),1)



    def test_cmp_stream(self):
        s1 = VariantSetFromFile(pkg_file('genomvar.test','data/example1.vcf'))
        s2 = VariantSetFromFile(pkg_file('genomvar.test','data/example2.vcf.gz'))
        nofv = 0
        for vrt in s1.diff_vrt(s2).iter_vrt():
            nofv += vrt.nof_unit_vrt()
        self.assertEqual(nofv,14)

    def test_diff_callback(self):
        s1 = VariantSetFromFile(pkg_file('genomvar.test','data/example3.vcf'))
        s2 = VariantSetFromFile(pkg_file('genomvar.test','data/example3.vcf'))
        cb = lambda m: [v.attrib['vcf_notation']['row'] for v in m]
        for N,vrt in enumerate(s1.comm_vrt(s2).iter_vrt(callback=cb)):
            self.assertEqual(vrt.attrib['vcf_notation']['row'],
                             vrt.attrib['cmp'][0])
        self.assertEqual(N,7)

    def test_ambig_difference_different_ambig(self):
        #         10043
        # Ref     T--CA--CA--G
        # v1 s1   T--CA--CACAG ins CA right
        # v2 s1   T------CA--G del CA left
        # v1 s2   TCACA--CA--G ins CA left
        # v2 s2   T--CA------G del CA right
        fac = variant.VariantFactory(reference=self.chr24,normindel=True)
        s1 = VariantSet.from_variants(
            [fac.from_edit('chr24',10047,'A','ACA'),
             fac.from_edit('chr24',10043,'TCA','T')])
        s2 = VariantSet.from_variants(
            [fac.from_edit('chr24',10043,'T','TCA'),
             fac.from_edit('chr24',10045,'ACA','A')])
        self.assertEqual(len(list(s1.diff(s2,match_ambig=True).iter_vrt())),0)

    def test_ambig_difference_snp_in_locus(self):
        #         10043
        # Ref     TC-ACA--G
        # v1 s1         CA
        # v1 s2    G
        # v2 s2     T
        fac = variant.VariantFactory(reference=self.chr24,normindel=True)
        s1 = VariantSet.from_variants(
            [fac.from_edit('chr24',10047,'A','ACA')])
        s2 = VariantSet.from_variants(
            [fac.from_edit('chr24',10044,'C','G'),
             fac.from_edit('chr24',10044,'C','CT')])
        self.assertEqual(len(list(s1.comm(s2,match_ambig=True).iter_vrt())),0)

    def test_differently_sorted_chroms(self):
        s1 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example3.vcf'))
        header = []
        variants = {}
        with open(pkg_file('genomvar.test','data/example3.vcf')) as fh:
            for line in fh:
                if line.startswith('#'):
                    header.append(line)
                else:
                    variants.setdefault(line.split(maxsplit=1)[0],[])\
                            .append(line)
        tf = tempfile.NamedTemporaryFile(suffix='.vcf')
        with open(tf.name,'wt') as fh:
            fh.writelines(header)
            for chrom in ['chr1','chr10','chr2']:
                fh.writelines(variants[chrom])
            
        s2 = VariantSetFromFile(tf.name)
        with self.assertRaises(DifferentlySortedChromsError):
            list(s1.diff_vrt(s2).iter_vrt())
            

    def test_cmp_vrt_iter_same(self):
        vs = VariantSetFromFile(pkg_file('genomvar.test','data/example2.vcf.gz'))
        tot = list(vs.find_vrt())
        # print(tot)
        comm = list(vs.comm_vrt(vs).iter_vrt())
        self.assertEqual(len(comm),len(tot))
    
    def test_cmp_vrt_iter_vrt(self):
        vs1 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            parse_samples=True)
        vs2 = VariantSetFromFile(pkg_file('genomvar.test','data/example2.vcf.gz'),parse_samples=True)
        comm = list()
        for vrt in vs1.comm_vrt(vs2).iter_vrt():
            comm.append(vrt)
            self.assertTrue(vrt.attrib['samples'],
                            msg='Vrt {} has no samples'.format(vrt))
        self.assertEqual(len(comm),4)
        diff = vs1.diff_vrt(vs2).iter_vrt()
        self.assertEqual(len(list(diff)),12)

    def test_cmp_vrt_iter_vrt2(self):
        vs1 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example1.vcf'))
        vs2 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example_gnomad_1.vcf.gz'))
        vrt = list(vs1.iter_vrt())
        self.assertEqual(len(list(vs1.diff_vrt(vs2).iter_vrt())),
                         len(vrt))
        
    def test_cmp_vrt_region(self):
        vs1 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            parse_samples=True,parse_info=True, index=True)
        vs2 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example2.vcf.gz'),
            parse_samples='SAMP1',parse_info=True, index=True)
        comm = list(vs1.comm_vrt(vs2).region(rgn='chr24:10040-10050'))
        self.assertEqual(len(comm),2)
        v1,v2 = comm
        self.assertEqual(v1.attrib['info']['AF'],1.0)
        self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(0,1))

    def test_ValueError_on_nonindexed(self):
        vs1 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example_1000genomes_1.vcf.gz'),
            parse_samples=True)
        vs2 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example_1000genomes_2.vcf.gz'),
            parse_samples=True)
        with self.assertRaises(NoIndexFoundError) as cm:
            list(vs1.comm_vrt(vs2).region(rgn='7:152134922-152436005'))

    def test_cmp_vrt_region_multisample2(self):
        vs1 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example_1000genomes_1.vcf.gz'),
            parse_samples=True, index=True)
        vs2 = VariantSetFromFile(
            pkg_file('genomvar.test','data/example_1000genomes_2.vcf.gz'),
            parse_samples=True, index=True)
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
        v1 = factory.from_edit('chr24',2093,'TGG','CCC')
        v2 = factory.from_edit('chr24',2098,'TT','GG')
        v3 = factory.from_edit('chr24',2098,'TT','CC')
        v4 = factory.from_edit('chr24',3200,'G','GG')
        hap = variant.Haplotype.from_variants([v1, v2])
        s1 = VariantSet.from_variants([hap, v3, v4])
        s2 = VariantSet.from_variants(s1.iter_vrt())
        vrt = list(s2.find_vrt('chr24', 2090, 2095, expand=True))
        self.assertEqual(len(vrt),2)
        self.assertEqual(s2.nof_unit_vrt(),8)
        self.assertEqual(s1.diff(s2).nof_unit_vrt(), 0)
        self.assertEqual(len(list(s2.diff_vrt(s1).iter_vrt())), 0)

        h1 = variant.Haplotype.from_variants([v1,v2])
        h2 = variant.Haplotype.from_variants([v3,v4])
        s3 = VariantSet.from_variants([h1,h2])
        
        self.assertEqual(s3.diff(s2).nof_unit_vrt(),0)
        self.assertEqual(s2.diff(s3).nof_unit_vrt(),0)

    def test_VariantSet_cmp(self):
        vcf1 = pkg_file('genomvar.test','data/example1.vcf')
        vcf2 = pkg_file('genomvar.test','data/example2.vcf.gz')
        s1 = VariantSet.from_vcf(vcf1,parse_info=True,parse_samples=True)
        s2 = VariantSet.from_vcf(vcf2)
        diff = s1.diff(s2)
        self.assertEqual(diff.nof_unit_vrt(),14)
        # Now same diff but without loading in memory
        N = 0
        for vrt in s1.diff_vrt(s2).iter_vrt():
            N += vrt.nof_unit_vrt()
        self.assertEqual(N, 14)
        
        comm = s1.comm(s2)
        self.assertEqual(len(list(comm.iter_vrt())),4)
        v1,v2 = sorted(comm.iter_vrt(),
                       key=lambda v: v.key)[:2]
        self.assertEqual(v1.attrib['info']['NSV'],1)
        self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(0,1))

        self.assertEqual(s2.comm(s1).nof_unit_vrt(), comm.nof_unit_vrt())

    def test_VariantSet_cmp2(self):
        vcf1 = pkg_file('genomvar.test','data/example3.vcf')
        s1 = VariantSet.from_vcf(vcf1)
        s2 = VariantSet.from_variants(
            reversed(list(s1.iter_vrt())))
        diff = s1.diff(s2)
        self.assertEqual(diff.nof_unit_vrt(),0)
        
        comm = s1.comm(s2)
        self.assertEqual(comm.nof_unit_vrt(),s1.nof_unit_vrt())
        self.assertEqual(s2.comm(s1).nof_unit_vrt(), comm.nof_unit_vrt())

    def test_VariantSet_no_common(self):
        vcf1 = pkg_file('genomvar.test','data/example1.vcf')
        vcf2 = pkg_file('genomvar.test','data/example_gnomad_1.vcf.gz')
        s1 = VariantSet.from_vcf(vcf1)
        s2 = VariantSet.from_vcf(vcf2)
        self.assertEqual(s1.comm(s2).nof_unit_vrt(), 0)
        
