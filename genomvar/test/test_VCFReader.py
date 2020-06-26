import unittest
import numpy as np
from genomvar.vcf import VCFReader
from pkg_resources import resource_filename as pkg_file

class TestVCFReaderCase(unittest.TestCase):
    def test_infer_dtypes(self):
        reader = VCFReader(pkg_file('genomvar.test','data/test_vcf.vcf'))
        self.assertEqual(reader.header_len,15)
        dtype = reader._parse_header_dtypes()
        self.assertEqual(len(dtype['format']),1)
        self.assertTrue(issubclass(dtype['format']['GT']['type'],np.object_),
                        msg='Got type'+str(dtype['format']['GT']['type']))
    
    def test_iter_chrom_rows(self):
        reader = VCFReader(pkg_file('genomvar.test','data/test_vcf2.vcf'))
        chroms = set()
        rows = {}
        for chrom,it in reader.iter_rows_by_chrom():
            chroms.add(chrom)
            l = list(it)
            rows[chrom] = len(l)
        self.assertEqual(chroms,{'chr15rgn','chr11rgn'})
        self.assertEqual(rows['chr11rgn'],5)
        self.assertEqual(rows['chr15rgn'],4)

        # same but without reading the chroms
        chroms = set()
        for chrom,it in reader.iter_rows_by_chrom():
            chroms.add(chrom)
        self.assertEqual(chroms,{'chr15rgn','chr11rgn'})

    def test_get_chroms(self):
        reader = VCFReader(pkg_file(
            'genomvar.test','data/test_vcf6.vcf.gz'))
        self.assertEqual(list(reader.get_chroms(unindexed=True)),
                         ['chr1','chr10','chr2'])

    def test_check_getting_vrt_is_sorted(self):
        reader = VCFReader(pkg_file(
            'genomvar.test','data/test_vcf12_gnomad.vcf.gz'),index=True)
        starts = [v.start for v in reader.iter_vrt()]
        self.assertEqual(starts,sorted(starts))

        starts2 = [v.start for v in reader.find_vrt(
            'chr15',74719587,74824401)]
        self.assertEqual(starts2,sorted(starts2))
        
    def test_iter_vrt_example1(self):
        reader = VCFReader(pkg_file('genomvar.test','data/test_vcf.vcf'))
        vrts = list(reader.iter_vrt(parse_info=True,parse_samples=True))
        vrt1,vrt2 = vrts[:2]

        # Test Ref and Alt
        self.assertEqual([vrt1.start,vrt1.ref,vrt1.alt],
                         [23,'G',''])
        self.assertEqual([vrt2.start,vrt2.ref,vrt2.alt],
                         [24,'','G'])

        # Check row numbers
        self.assertEqual(vrt1.attrib['vcf_notation']['row'],0)
        self.assertEqual(vrt2.attrib['vcf_notation']['row'],0)
        # Test INFO
        self.assertEqual(vrt1.attrib['info']['AF'],0.5)
        self.assertEqual(vrt2.attrib['info']['AF'],0.5)

        # Test SAMPLES fields
        self.assertEqual(vrt1.attrib['samples']['SAMP1']['GT'],(0,1,0))
        self.assertEqual(vrt2.attrib['samples']['SAMP1']['GT'],(0,0,1))

    def test_iter_vrt_example2(self):
        reader = VCFReader(pkg_file('genomvar.test','data/test_vcf7.vcf'))
        for vrt in reader.iter_vrt(parse_info=True,parse_samples=True):
            self.assertTrue(bool(vrt.attrib['samples']['TUMOR']))

        recs = reader.get_records(parse_info=True,parse_samples=True)
        inforow2 = recs['info'][1]
        self.assertTrue(np.isnan(inforow2['MQ0']))

        samprow2 = recs['sampdata']['NORMAL'][1]
        self.assertTrue(np.isnan(samprow2['FDP']))
        self.assertTrue(np.isnan(samprow2['AU'][1]))

    def test_iter_vrt_gzipped(self):
        reader = VCFReader(pkg_file('genomvar.test','data/test_vcf2.vcf.gz'),
                           index=True)
        self.assertEqual(list(reader.chroms),['chr11rgn','chr15rgn'])
        for vrt in reader.iter_vrt():
            self.assertIn(vrt.chrom,['chr15rgn','chr11rgn'])

        self.assertEqual(len(list(reader.find_vrt(chrom='chr15rgn'))),4)
        self.assertEqual(len(list(
            reader.find_vrt('chr11rgn',7464,7465))),3)
