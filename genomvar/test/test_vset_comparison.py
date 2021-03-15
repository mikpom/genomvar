import operator
from unittest import TestCase
from pkg_resources import resource_filename as pkg_file
import itertools
import tempfile
from genomvar import variant
from genomvar.variant import VariantFactory
from genomvar.varset import VariantSet,MutableVariantSet,\
    VariantFileSet,IndexedVariantFileSet
from genomvar import DifferentlySortedChromsError
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
            reference=pkg_file('genomvar.test','data/chr15rgn.fna'),
            normindel=True)
        v1 = factory.from_edit('chr15rgn',10043,'T','TCA')
        v2 = factory.from_edit('chr15rgn',10045,'A','ACA')
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

        s1 = MutableVariantSet(reference=self.chr15rgn)
        vb = factory.from_edit('chr15rgn',1206,'G','C')
        s1.add_vrt(vb,GT=(1,0))
        vb = factory.from_edit('chr15rgn',2093,'TG','CC')
        s1.add_vrt(vb,GT=(1,1))
        vb = factory.from_edit('chr15rgn',2095,'G','C')
        s1.add_vrt(vb,GT=(1,1))
        vb = factory.from_edit('chr15rgn',10044,'C','T')
        s1.add_vrt(vb,GT=(1,0))
        vb = factory.from_edit('chr15rgn',10044,'C','G')
        s1.add_vrt(vb,GT=(0,1))

        s2 = MutableVariantSet(reference=self.chr15rgn)
        vb = factory.from_edit('chr15rgn',2092,'TT','AC')
        s2.add_vrt(vb,GT=(0,1),allow_adjust_genotype=True)
        vb = factory.from_edit('chr15rgn',10044,'C','T')
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
        s1 = MutableVariantSet.from_vcf(
            pkg_file('genomvar.test','data/example1.vcf'),
            reference=self.chr15rgn,normindel=True,parse_info=True,
            sample='SAMP1')
        _vrt = list(s1.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(_vrt),2)
        v1,v2 = _vrt
        self.assertEqual(v1.attrib['info']['AF'],1.0)

        s2 = MutableVariantSet(reference=self.chr15rgn)
        s2.add_vrt(factory.from_edit('chr15rgn',1206,'G','C'),
                   GT=(0,1))

        diff = s1.diff(s2)
        _vrt = list(diff.find_vrt('chr15rgn',1200,1210))
        self.assertEqual(len(_vrt),1)
        v = _vrt[0]
        self.assertEqual(v.attrib['filter'],['PASS'])
        self.assertEqual(v.attrib['info']['AF'],1.0)

    def test_strip_order_dependent_Ambig(self):
        #    10043
        #    TCACAG
        # s1 TC--AG
        # s2 T--CAG
        vfac = VariantFactory(reference=self.chr15rgn,normindel=True)
        
        s1 = MutableVariantSet(reference=self.chr15rgn)
        vb = vfac.from_edit('chr15rgn',10044,'CAC','C')
        s1.add_vrt(vb,GT=(0,1))
        s2 = MutableVariantSet(reference=self.chr15rgn)
        vb = vfac.from_edit('chr15rgn',10043,'TCA','T')
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
        s1 = MutableVariantSet(reference=self.chr15rgn)
        vb = factory.from_edit('chr15rgn',9947,'TT','T')
        s1.add_vrt(vb,GT=(0,1))
        s2 = MutableVariantSet(reference=self.chr15rgn)
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

        s1 = MutableVariantSet(reference=self.chr15rgn)
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

        s2 = MutableVariantSet(reference=self.chr15rgn)
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
        

    def test_mnp_com_split(self):
        #                           23
        #    TTCACTTAGCATAATGTCTTCAAGATT
        # v1                       TT -single
        # v2                       TT -splitted
        vset1 = MutableVariantSet(reference=self.chr15rgn)
        vb = factory.from_edit('chr15rgn',22,'AG','TT')
        vset1.add_vrt(vb,GT=(1,0))
        vset2 = MutableVariantSet(reference=self.chr15rgn)
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
        s1 = MutableVariantSet(reference=self.chr15rgn)
        vb = factory.from_edit('chr15rgn',2,'CAC','GGG')
        s1.add_vrt(vb,GT=(1,0))
        s2 = MutableVariantSet(reference=self.chr15rgn)
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
        vcf = pkg_file('genomvar.test','data/example1.vcf')
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
        s1 = VariantFileSet(pkg_file('genomvar.test','data/example1.vcf'))
        s2 = VariantFileSet(pkg_file('genomvar.test','data/example2.vcf.gz'))
        nofv = 0
        for vrt in s1.diff_vrt(s2).iter_vrt():
            nofv += vrt.nof_unit_vrt()
        self.assertEqual(nofv,14)

    def test_diff_callback(self):
        s1 = VariantFileSet(pkg_file('genomvar.test','data/example3.vcf'))
        s2 = VariantFileSet(pkg_file('genomvar.test','data/example3.vcf'))
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
        fac = variant.VariantFactory(reference=self.chr15rgn,normindel=True)
        s1 = VariantSet.from_variants(
            [fac.from_edit('chr15rgn',10047,'A','ACA'),
             fac.from_edit('chr15rgn',10043,'TCA','T')])
        s2 = VariantSet.from_variants(
            [fac.from_edit('chr15rgn',10043,'T','TCA'),
             fac.from_edit('chr15rgn',10045,'ACA','A')])
        self.assertEqual(len(list(s1.diff(s2,match_ambig=True).iter_vrt())),0)

    def test_ambig_difference_snp_in_locus(self):
        #         10043
        # Ref     TC-ACA--G
        # v1 s1         CA
        # v1 s2    G
        # v2 s2     T
        fac = variant.VariantFactory(reference=self.chr15rgn,normindel=True)
        s1 = VariantSet.from_variants(
            [fac.from_edit('chr15rgn',10047,'A','ACA')])
        s2 = VariantSet.from_variants(
            [fac.from_edit('chr15rgn',10044,'C','G'),
             fac.from_edit('chr15rgn',10044,'C','CT')])
        self.assertEqual(len(list(s1.comm(s2,match_ambig=True).iter_vrt())),0)

    def test_differently_sorted_chroms(self):
        s1 = VariantFileSet(
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
            
        s2 = VariantFileSet(tf.name)
        with self.assertRaises(DifferentlySortedChromsError):
            list(s1.diff_vrt(s2).iter_vrt())
            

    def test_cmp_vrt_iter_same(self):
        vs = VariantFileSet(pkg_file('genomvar.test','data/example2.vcf.gz'))
        tot = list(vs.find_vrt())
        # print(tot)
        comm = list(vs.comm_vrt(vs).iter_vrt())
        self.assertEqual(len(comm),len(tot))
    
    def test_cmp_vrt_iter_vrt(self):
        vs1 = VariantFileSet(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            parse_samples=True)
        vs2 = VariantFileSet(pkg_file('genomvar.test','data/example2.vcf.gz'),parse_samples=True)
        comm = list()
        for vrt in vs1.comm_vrt(vs2).iter_vrt():
            comm.append(vrt)
            self.assertTrue(vrt.attrib['samples'],
                            msg='Vrt {} has no samples'.format(vrt))
        self.assertEqual(len(comm),4)
        diff = vs1.diff_vrt(vs2).iter_vrt()
        self.assertEqual(len(list(diff)),12)

    def test_cmp_vrt_iter_vrt2(self):
        vs1 = VariantFileSet(
            pkg_file('genomvar.test','data/example1.vcf'))
        vs2 = VariantFileSet(
            pkg_file('genomvar.test','data/example_gnomad_1.vcf.gz'))
        vrt = list(vs1.iter_vrt())
        self.assertEqual(len(list(vs1.diff_vrt(vs2).iter_vrt())),
                         len(vrt))
        
    def test_cmp_vrt_region(self):
        vs1 = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example1.vcf.gz'),
            parse_samples=True,parse_info=True)
        vs2 = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example2.vcf.gz'),
            parse_samples='SAMP1',parse_info=True)
        comm = list(vs1.comm_vrt(vs2).region(rgn='chr15rgn:10040-10050'))
        self.assertEqual(len(comm),2)
        v1,v2 = comm
        self.assertEqual(v1.attrib['info']['AF'],1.0)
        self.assertEqual(v1.attrib['samples']['SAMP1']['GT'],(0,1))

    def test_ValueError_on_nonindexed(self):
        vs1 = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example_1000genomes_1.vcf.gz'),
            parse_samples=True)
        vs2 = VariantFileSet(
            pkg_file('genomvar.test','data/example_1000genomes_2.vcf.gz'),
            parse_samples=True)
        with self.assertRaises(ValueError) as cm:
            list(vs1.comm_vrt(vs2).region(rgn='7:152134922-152436005'))
        self.assertIn('region access',cm.exception.args[0])

    def test_cmp_vrt_region_multisample2(self):
        vs1 = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example_1000genomes_1.vcf.gz'),
            parse_samples=True)
        vs2 = IndexedVariantFileSet(
            pkg_file('genomvar.test','data/example_1000genomes_2.vcf.gz'),
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
        s1 = MutableVariantSet(reference=self.chr15rgn)
        s1.add_hap_variants([v1,v2],GT=(1,0))
        s1.add_vrt(v3,GT=(0,1))
        s1.add_vrt(v4,GT=(0,1))
        s2 = VariantSet.from_variants(s1.iter_vrt())
        vrt = list(s2.find_vrt('chr15rgn', 2090, 2095, expand=True))
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
        
