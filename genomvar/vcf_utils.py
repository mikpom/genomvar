import os
from jinja2 import Environment,FileSystemLoader

VCF_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",  "INFO",
              "FORMAT",  "SAMPLES"]


class VCFRow(object):
    """Tuple-like object storing variant information in VCF-like form.

    str() returns a string, formatted as a row in VCF file."""
    def __init__(self,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,
                 FORMAT=None,SAMPLES=None,rnum=None):
        self.CHROM = CHROM
        self.POS = int(POS)
        self.ID = ID
        self.REF = REF
        self.ALT = ALT
        self.QUAL = QUAL
        self.FILTER = FILTER
        self.INFO = INFO
        self.FORMAT = FORMAT
        self.SAMPLES = SAMPLES
        self.rnum = rnum

    __slots__ = [*VCF_fields, 'rnum']

    @staticmethod
    def _to_string(v):
        if v is None:
            return '.'
        else:
            return str(v)

    def __repr__(self):
        return '<VCFRow {}:{} {}->{}>'\
            .format(self.CHROM,self.POS,self.REF,self.ALT)
    def __str__(self):
        fields = [self._to_string(getattr(self, a)) for a in VCF_fields[:8]]
        if not self.FORMAT is None:
            fields += [self._to_string(self.FORMAT),str(self.SAMPLES)]
        return '\t'.join(fields)

tmpl_dir = os.path.join(os.path.dirname(__file__),'tmpl')
env = Environment(
    loader=FileSystemLoader(tmpl_dir),
)
header_simple = env.get_template('vcf_head_simple.tmpl')
header = env.get_template('vcf_head.tmpl')
row_tmpl = env.get_template('vcf_row.tmpl')
