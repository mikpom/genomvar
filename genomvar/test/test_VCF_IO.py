import warnings
import unittest
import numpy as np
import io
from genomvar.test import MyTestCase
from genomvar import variant
from genomvar.vcf import VCFReader, BCFReader, VCFWriter, \
    RESERVED_FORMAT, RESERVED_INFO
from genomvar.vcf_utils import header as vcf_header
from pkg_resources import resource_filename as pkg_file

class TestVCFWriterCase(MyTestCase):
    def test_writer_instantiation(self):
        deletion = variant.Del("chr24",5,10)
        writer = VCFWriter(reference=self.chr24.filename)
        row = writer.get_row(deletion)
        self.assertEqual(
            row.REF,
            self.chr24.get(deletion.chrom, deletion.start-1, deletion.end))

        writer.reference.close()

    def test_minimal_VCF_definition_io(self):
        buf = io.StringIO()
        with open(pkg_file('genomvar.test','data/example1.vcf'), 'rt') as fh:
            for line in fh:
                if line.startswith('##fileformat') \
                          or line.startswith('#CHROM') \
                          or not line.startswith('#'):
                    buf.write(line)
        
        buf.seek(0)
        reader = VCFReader(buf)

        outbuf = io.StringIO()
        writer = VCFWriter(format_spec=[RESERVED_FORMAT.GT],
                           samples=reader.samples)
        variants1 = []
        for vrt in reader.iter_vrt(parse_samples=True):
            self.assertTrue(
                isinstance(vrt.attrib['samples']['SAMP1']['GT'], str))
            if vrt.attrib['samples']['SAMP1'].get('GT')=='0/1':
                vrt.attrib['samples']['SAMP1']['GT'] = (0, 1)
            else:
                vrt.attrib['samples']['SAMP1']['GT'] = None
            outbuf.write(str(writer.get_row(vrt)))
            variants1.append(vrt)
        variants1.sort(key=lambda v: v.start)
        
        outbuf.seek(0)
        variants2 = list(VCFReader(outbuf).iter_vrt())
        variants2.sort(key=lambda v: v.start)
        
        for v1, v2 in zip(variants1, variants2):
            self.assertTrue(v1.edit_equal(v2))
        

class TestVCFReaderCase(unittest.TestCase):
    def test_init(self):
        reader = VCFReader(pkg_file('genomvar.test','data/example1.vcf'))
        self.assertEqual(reader.header_len,15)
        dtype = reader._dtype
        self.assertEqual(len(dtype['format']),1)
        self.assertTrue(issubclass(dtype['format']['GT']['dtype'],np.object_),
                        msg='Got type'+str(dtype['format']['GT']['type']))
    
    def test_iter_chrom_rows(self):
        reader = VCFReader(pkg_file('genomvar.test','data/example2.vcf.gz'))
        chroms = set()
        rows = {}
        for chrom,it in reader.iter_rows_by_chrom():
            chroms.add(chrom)
            l = list(it)
            rows[chrom] = len(l)
        self.assertEqual(chroms,{'chr24','chr23'})
        self.assertEqual(rows['chr23'],5)
        self.assertEqual(rows['chr24'],4)

        # same but without reading the chroms
        chroms = set()
        for chrom,it in reader.iter_rows_by_chrom():
            chroms.add(chrom)
        self.assertEqual(chroms,{'chr24','chr23'})

    def test_example3(self):
        reader = VCFReader(pkg_file(
            'genomvar.test','data/example3.vcf'))
        self.assertEqual(list(reader.get_chroms(allow_no_index=True)),
                         ['chr1','chr2','chr10'])
        vrt = list(reader.iter_vrt(parse_info=True,parse_samples=True))
        self.assertGreater(len(vrt),0)

        v = vrt[3]
        self.assertEqual(v.attrib['id'], None)

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
        self.assertEqual(reader.samples, ['SAMP1'])
        vrts = list(reader.iter_vrt(parse_info=True, parse_samples=True))
        vrt1,vrt2 = vrts[:2]

        # Test Ref and Alt
        self.assertEqual([vrt1.start,vrt1.ref,vrt1.alt],
                         [23,'G',''])
        self.assertEqual([vrt2.start,vrt2.ref,vrt2.alt],
                         [24,'','G'])

        # Check row numbers
        self.assertEqual(vrt1.attrib['vcf_notation']['row'],0)
        self.assertEqual(vrt2.attrib['vcf_notation']['row'],0)
        self.assertEqual(vrt1.attrib['allele_num'], 0)
        self.assertEqual(vrt2.attrib['allele_num'], 1)
        # Test INFO
        self.assertEqual(vrt1.attrib['info']['AF'],0.5)
        self.assertEqual(vrt2.attrib['info']['AF'],0.5)

        # Test SAMPLES fields
        self.assertEqual(vrt1.attrib['samples']['SAMP1']['GT'],(0,1,0))
        self.assertEqual(vrt2.attrib['samples']['SAMP1']['GT'],(0,0,1))

    def test_iter_vrt_gzipped(self):
        reader = VCFReader(pkg_file('genomvar.test','data/example2.vcf.gz'),
                           index=True)
        self.assertEqual(list(reader.chroms),['chr23','chr24'])
        for vrt in reader.iter_vrt():
            self.assertIn(vrt.chrom,['chr24','chr23'])

        self.assertEqual(len(list(reader.find_vrt(chrom='chr24'))),4)
        self.assertEqual(len(list(
            reader.find_vrt('chr23',7464,7465))),3)

    def test_from_vcf_missing_values(self):
        buf = io.StringIO()
        format_fields = (RESERVED_FORMAT.AD, RESERVED_FORMAT.DP, RESERVED_FORMAT.GT)
        header = vcf_header.render(
            samples=['S1','S2'],
            ctg_len={},
            format=[{k:getattr(spec, k.upper()) for k in \
                     ('name', 'number', 'type', 'description')} \
                    for spec in format_fields])

        buf.write(header)
        buf.write('chr15\t17017413\t.\tA\tG\t38\t.\t.\tGT\t./.\t0/1\n')
        buf.write('chr15\t17017413\t.\tA\tG\t38\t.\t.\tGT:AD\t./.\t0/1:10,.\n')
        buf.write('chr15\t17017413\t.\tA\tG\t38\t.\t.\tGT:DP\t./.\t1/1:.\n')
        buf.seek(0)
        buf.seek(0)
        vs = VCFReader(buf)
        
        v = list(vs.iter_vrt(parse_samples=True))[0]
        self.assertEqual(v.attrib['samples']['S1']['GT'], (None,None))

    def test_sv_types(self):
        reader = VCFReader(pkg_file('genomvar.test', 'data/example4.vcf.gz'))

        with warnings.catch_warnings(record=True) as wrn:
            # warnings.simplefilter(append=True)
            for cnt,vrt in enumerate(reader.iter_vrt()):
                pass 
            self.assertEqual(cnt, 99)
            self.assertGreater(len(wrn), 1)
            self.assertIn('Structural', str(wrn[-1].message))
            
        

    def test_bcf_format(self):
        reader = BCFReader(
            pkg_file('genomvar.test', 'data/example3.bcf'))

        variants = list(reader.iter_vrt(parse_info=True,parse_samples=True))
        self.assertEqual(len(variants), 8)

        v1 = variants[0]
        self.assertEqual([v1.start, v1.ref,v1.alt],
                         [23, 'G',   ''])
        self.assertEqual([v1.attrib[a] for a in ('id', 'filter')],
                         ['1',['PASS']])
        self.assertEqual(v1.attrib['info']['AF'], 0.5)

        v2 = variants[1]
        self.assertEqual(v2.attrib['samples']['SAMP1']['AD'], 10)

