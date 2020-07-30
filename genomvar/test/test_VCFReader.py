import unittest
import numpy as np
import io
from genomvar.vcf import VCFReader, header as vcf_header
from pkg_resources import resource_filename as pkg_file

class TestVCFReaderCase(unittest.TestCase):
    def test_init(self):
        reader = VCFReader(pkg_file('genomvar.test','data/example1.vcf'))
        self.assertEqual(reader.header_len,15)
        dtype = reader._dtype
        self.assertEqual(len(dtype['format']),1)
        self.assertTrue(issubclass(dtype['format']['GT']['type'],np.object_),
                        msg='Got type'+str(dtype['format']['GT']['type']))
    
    def test_iter_chrom_rows(self):
        reader = VCFReader(pkg_file('genomvar.test','data/example2.vcf.gz'))
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

    def test_example3(self):
        reader = VCFReader(pkg_file(
            'genomvar.test','data/example3.vcf'))
        self.assertEqual(list(reader.get_chroms(unindexed=True)),
                         ['chr1','chr2','chr10'])
        vrt = list(reader.iter_vrt(parse_info=True,parse_samples=True))
        self.assertGreater(len(vrt),0)

    def test_check_getting_vrt_is_sorted(self):
        reader = VCFReader(pkg_file(
            'genomvar.test','data/example_gnomad_2.vcf.gz'),index=True)
        starts = [v.start for v in reader.iter_vrt()]
        self.assertEqual(starts,sorted(starts))

        starts2 = [v.start for v in reader.find_vrt(
            'chr15',74719587,74824401)]
        self.assertEqual(starts2,sorted(starts2))
        
    def test_iter_vrt_example1(self):
        reader = VCFReader(pkg_file('genomvar.test','data/example1.vcf'))
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

    def test_iter_vrt_gzipped(self):
        reader = VCFReader(pkg_file('genomvar.test','data/example2.vcf.gz'),
                           index=True)
        self.assertEqual(list(reader.chroms),['chr11rgn','chr15rgn'])
        for vrt in reader.iter_vrt():
            self.assertIn(vrt.chrom,['chr15rgn','chr11rgn'])

        self.assertEqual(len(list(reader.find_vrt(chrom='chr15rgn'))),4)
        self.assertEqual(len(list(
            reader.find_vrt('chr11rgn',7464,7465))),3)

    def test_from_vcf_missing_values(self):
        buf = io.StringIO()
        header = vcf_header.render(
            samp='\t'.join(['S1','S2']),
            ctg_len={},
            format=[{'name':'AD','number':'R',
                     'type':'Integer','description':'""'},
                    {'name':'DP','number':'1',
                     'type':'Integer','description':'""'}])
        buf.write(header)
        buf.write('chr15\t17017413\t.\tA\tG\t38\t.\t.\tGT\t./.\t0/1\n')
        buf.write('chr15\t17017413\t.\tA\tG\t38\t.\t.\tGT:AD\t./.\t0/1:10,.\n')
        buf.write('chr15\t17017413\t.\tA\tG\t38\t.\t.\tGT:DP\t./.\t1/1:.\n')
        buf.seek(0)
        vs = VCFReader(buf)
        
        v = list(vs.iter_vrt(parse_samples=True))[0]
        self.assertEqual(v.attrib['samples']['S1']['GT'], (None,None))
