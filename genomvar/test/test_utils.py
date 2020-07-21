from pkg_resources import resource_filename as pkg_file
from unittest import TestCase
from genomvar.variantops import novlp_grp,novlp_cmbtns,_mnp_isect,\
    _cmp_mnps,nof_snp_vrt
from genomvar import variant
from genomvar.variant import VariantBase,GenomVariant,AmbigIndel
from genomvar.vcf import VCFReader
from genomvar.utils import zip_variants,_strip_ref_alt,no_ovlp

class test_MNPs(TestCase):
    # R         TTATATT
    #           -------
    # v0        AGTCTGA
    # v1,v3     AGT TG
    # v2,v4      G    A

    v0 = variant.MNP(chrom='1',start=1,alt='AGTCTGA')
    v1 = variant.MNP(chrom='1',start=1,alt='AGT')
    v2 = variant.SNP(chrom='1',start=2,alt='G')
    v3 = variant.MNP(chrom='1',start=5,alt='TG')
    v4 = variant.SNP(chrom='1',start=7,alt='A')

    def test_overlaps(self):
        cmb = novlp_cmbtns([self.v1,self.v2,self.v3,self.v4])
        self.assertEqual(len(cmb),2)

        grp = list(novlp_grp([self.v1,self.v2,self.v3,self.v4]))
        self.assertEqual(len(grp),3)
        self.assertEqual(len(grp[0]),2)

    def test_mnp_isect(self):
        isect = _mnp_isect(self.v0,[self.v1,self.v3])
        self.assertEqual(''.join(isect),'AGT-TG-')
        isect = _mnp_isect(self.v0,[self.v4])
        self.assertEqual(''.join(isect),'------A')

class test_utils(TestCase):
    def test_zip_variants(self):
        fl = pkg_file('genomvar.test','data/test_vcf5.vcf')
        it1 = VCFReader(fl).iter_vrt()
        it2 = VCFReader(fl).iter_vrt()
        matches = list(zip_variants(it1,it2))
        d1,v1,ovlp1 = matches[0]
        self.assertEqual(d1,0)
        self.assertEqual(v1.alt,'')
        self.assertEqual([v.ref for v in ovlp1],['G',''])

        dlast,vlast,ovlplast = matches[-1]
        self.assertEqual(vlast.alt,'CC')
        self.assertEqual([v.ref for v in ovlplast],['TG'])
        
    def test_cmp_mnps(self):
        v1 = variant.MNP(chrom='',start=1,end=4,alt='GGG')
        v2 = variant.SNP(chrom='',start=2,alt='G')
        d = _cmp_mnps(v1,[v2],action='diff')
        self.assertEqual(len(d),2)
        l,r = d
        self.assertEqual([l.start,l.alt],[1,'G'])
        self.assertEqual([r.start,r.alt],[3,'G'])

    def test_strip_ref_alt(self):
        self.assertEqual(_strip_ref_alt('GACTGC','G'),
                         ('ACTGC','',1))
        self.assertEqual(_strip_ref_alt('A','C'),
                         ('A','C',0))
        self.assertEqual(_strip_ref_alt('C','CTA'),
                         ('','TA',1))

    def test_min_mnps(self):
        factory = variant.VariantFactory()
        v1 = factory.from_edit('1',1,'AAA','GGG')
        v2 = factory.from_edit('1',1,'AAA','GGC')
        v3 = factory.from_edit('1',1,'AAA','CGG')
        self.assertEqual(nof_snp_vrt([v1,v2,v3]),5)

    def test_no_ovlp(self):
        # REF      TGG   TT   
        #          2093  2099 
        #          CCC   GG
        #                CG   
        factory = variant.VariantFactory()
        variants = [
            factory.from_edit('1',2093,'TGG','CCC'),
            factory.from_edit('1',2098,'TT','GG'),
            factory.from_edit('1',2098,'TT','CG'),
            factory.from_edit('2',3200,'G','GG')]

        for ind,chunk in enumerate(no_ovlp(variants)):
            if ind==0:
                self.assertEqual(len(chunk),1)
            elif ind==1:
                self.assertEqual(len(chunk),2)
            if ind==2:
                self.assertEqual(len(chunk),1)
        
        
        
