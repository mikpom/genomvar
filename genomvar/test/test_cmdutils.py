from unittest import TestCase
from pkg_resources import resource_filename as pkg_file
import io
import warnings
import itertools
from genomvar.cmdutils import _cmp_vcf
from genomvar.vcf import VCFRow
import tempfile
from genomvar import UnsortedVariantFileError

class TestStreamCmp(TestCase):
    f1 = pkg_file('genomvar.test','data/example1.vcf.gz')
    f2 = pkg_file('genomvar.test','data/example2.vcf.gz')
    def test_cmp_vcf_files(self):
        def _get_info(info):
            if info=='.':
                return {}
            tokenized = info.split(';')
            kval = map(lambda i: i.split('=',maxsplit=1),tokenized)
            return {k:v for (k,v) in kval}
        out = io.StringIO()
        with warnings.catch_warnings(record=True):
            cnt = _cmp_vcf(self.f1,self.f2,out=out)
        self.assertEqual(cnt[0], 14)
        self.assertEqual(cnt[2], 4)
        self.assertEqual(cnt[1], 12)
        
        out.seek(0)
        noheader = itertools.dropwhile(lambda l: l.startswith('#'),out)
        rows = [VCFRow(*l.strip().split('\t')) for l in noheader]

        row0 = rows[0]
        info = _get_info(row0.INFO)
        self.assertEqual([row0.CHROM,row0.POS,row0.REF,row0.ALT],
                         ['chr23',7462,'G','T'])
        self.assertEqual(info['whichVCF'],'second')
        self.assertEqual(info['ln'],'13')

        #last = rows[-1]
        info = _get_info(rows[-1].INFO)
        self.assertEqual(info['ln'],'30')
        self.assertEqual(info['ln2'],'21')
        
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
        out = io.StringIO()
        with warnings.catch_warnings(record=True):
            with self.assertRaises(UnsortedVariantFileError):
                _cmp_vcf(pkg_file('genomvar.test','data/example1.vcf'),
                         tf.name,out=out)
