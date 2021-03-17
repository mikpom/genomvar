import os
from collections import namedtuple
from jinja2 import Environment,FileSystemLoader
import numpy as np

VCF_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",  "INFO",
              "FORMAT",  "SAMPLES"]
VCF_info_fields = ["NAME",  "NUMBER",  "TYPE",  "DESCRIPTION",  "SOURCE",  "VERSION"]
VCFInfoSpec = namedtuple('VCFInfoSpec', VCF_info_fields)
# Map of VCF types to NumPy types
string2dtype = {'Float':np.float64,'Integer':np.int64,
                'String':np.object_,'Flag':np.bool_,
                'Character':'U1'}
# inverse for writing VCFs
dtype2string = {v:k for k,v in string2dtype.items()}


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

def _make_field_writer_func(tp, num):
    """tp can be a string in VCF notation or internal NumPy types."""
    if isinstance(tp, str):
        try:
            tp = string2dtype[tp]
        except KeyError:
            raise ValueError(
                'Allowed types "Type" in INFO spec are {} while found {}'\
                .format(','.join(string2dtype), tp))
    if tp==np.bool_:
        return lambda k,v: str(k)
    if tp==np.int_:
        tostring = lambda v: str(int(v)) if not np.isnan(v) else '.'
    else:
        tostring = lambda v: str(v) if not v is None else '.'
    if num in ('.', 'G') or (isinstance(num, int) and num>1):
        return lambda k,v: '{}={}'.format(
            k, ','.join(map(tostring, v)) \
            if not v is None else '.')
    elif num in (1, 'A'): # Alleles are split on read-in
        return lambda k,v: '{}={}'.format(
            k, tostring(v) if not v is None else '.')
    else:
        allowed = 'A, G, . or and integer'
        raise ValueError(
            'Allowed types for "Number" INFO spec are {} while found {}'\
            .format(allowed, num))

def issequence(seq):
    """
    Is seq a sequence (ndarray, list or tuple)?

    """
    return isinstance(seq, (np.ndarray, tuple, list))
    
def field_writer_simple(k,v):
    if issequence(v):
        res = '{}={}'.format(k, ','.join(map(str, v)))
    else:
        res = '{}={}'.format(k, str(v))
    return res
