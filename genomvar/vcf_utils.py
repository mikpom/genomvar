import os
from collections import namedtuple
from jinja2 import Environment,FileSystemLoader
import numpy as np

VCF_INFO_OR_FORMAT_DEF_FIELDS = \
    ["NAME",  "NUMBER",  "TYPE",  "DESCRIPTION",  "SOURCE",  "VERSION"]

VCF_INFO_OR_FORMAT_SPEC = namedtuple(
    'VCF_INFO_OR_FORMAT_SPEC',
    VCF_INFO_OR_FORMAT_DEF_FIELDS)

RESERVED_FORMAT_SPECS = \
[VCF_INFO_OR_FORMAT_SPEC(*t, None, None) for t in
 [
    ("AD", "R", "Integer", "Read depth for each allele"),
    ("ADF", "R", "Integer", "Read depth for each allele on the forward strand"),
    ("ADR", "R", "Integer", "Read depth for each allele on the reverse strand"),
    ("DP", "1", "Integer", "Read depth"),
    ("EC", "A", "Integer", "Expected alternate allele counts"),
    ("FT", "1", "String", "Filter indicating if this genotype was “called”"),
    ("GL", "G", "Float", "Genotype likelihoods"),
    ("GP", "G", "Float", "Genotype posterior probabilities"),
    ("GQ", "1", "Integer", "Conditional genotype quality"),
    ("GT", "1", "String", "Genotype"),
    ("HQ", "2", "Integer", "Haplotype quality"),
    ("MQ", "1", "Integer", "RMS mapping quality"),
    ("PL", "G", "Integer", "Phred-scaled genotype likelihoods rounded"\
                          +" to the closest integer")]]
RESERVED_INFO_SPECS = \
    [VCF_INFO_OR_FORMAT_SPEC(*t, None, None) for t in
     [('AA', '1', 'String', 'Ancestral allele'),
      ('AC', 'A', 'Integer', 'Allele count in genotypes, for each ALT allele, in the same order as listed'),
      ('AD', 'R', 'Integer', 'Total read depth for each allele'),
      ('ADF', 'R', 'Integer', 'Read depth for each allele on the forward strand'),
      ('ADR', 'R', 'Integer', 'Read depth for each allele on the reverse strand'),
      ('AF', 'A', 'Float', 'Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)'),
      ('AN', '1', 'Integer', 'Total number of alleles in called genotypes'),
      ('BQ', '1', 'Float', 'RMS base quality'),
      ('CIGAR', 'A', 'String', 'Cigar string describing how to align an alternate allele to the reference allele)'),
      ('DB', '0', 'Flag', 'dbSNP membership'),
      ('DP', '1', 'Integer', 'Combined depth across samples'),
      ('END', '1', 'Integer', 'End position on CHROM (used with symbolic alleles; see below)'),
      ('H2', '0', 'Flag', 'HapMap2 membership'),
      ('H3', '0', 'Flag', 'HapMap3 membership'),
      ('MQ', '1', 'Float', 'RMS mapping quality'),
      ('MQ0', '1', 'Integer', 'Number of MAPQ == 0 reads'),
      ('NS', '1', 'Integer', 'Number of samples with data'),
      ('SB', '4', 'Integer', 'Strand bias'),
      ('SOMATIC', '0', 'Flag', 'Somatic mutation (for cancer genomics)'),
      ('VALIDATED', '0', 'Flag', 'Validated by follow-up experiment'),
      ('1000G', '0', 'Flag', '1000 Genomes mem')]]

VCF_FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",  "INFO",
              "FORMAT",  "SAMPLES"]


# Map of VCF types to NumPy types
string2dtype = {'Float':np.float64,'Integer':np.int64,
                'String':np.object_,'Flag':np.bool_,
                'Character':'U1'}
# inverse for writing VCFs
dtype2string = {v:k for k,v in string2dtype.items()}


class VCFRow(object):
    """Tuple-like object storing variant information in VCF-like form.

    str() returns a string, formatted as a row in VCF file."""
    MANDATORY_FIELDS = VCF_FIELDS[:8]
    
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

    __slots__ = [*VCF_FIELDS, 'rnum']

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
        if not self.FORMAT is None:
            return '\t'.join(self._to_string(getattr(self, a)) for a in \
                      VCF_FIELDS)
        else:
            return '\t'.join(self._to_string(getattr(self, a)) for a in \
                      self.MANDATORY_FIELDS)

tmpl_dir = os.path.join(os.path.dirname(__file__),'tmpl')
env = Environment(
    loader=FileSystemLoader(tmpl_dir),
)
header_simple = env.get_template('vcf_head_simple.tmpl')
header = env.get_template('vcf_head.tmpl')
row_tmpl = env.get_template('vcf_row.tmpl')

def _make_field_writer_func(tp, num, info_or_format):
    """tp can be a string in VCF notation or internal NumPy types."""
    if isinstance(tp, str):
        try:
            tp = string2dtype[tp]
        except KeyError:
            raise ValueError(
                'Allowed types "Type" in INFO spec are {} while found {}'\
                .format(','.join(string2dtype), tp))
    if tp==np.bool_:
        if info_or_format=='info':
            return lambda k,v: str(k)
        elif info_or_format=='format':
            return lambda k,v: '1'
    if tp==np.int_:
        tostring = lambda v: str(int(v)) if not np.isnan(v) else '.'
    else:
        tostring = lambda v: str(v) if not v is None else '.'

    if num in ('.', 'G') or (isinstance(num, int) and num>1) \
                        or (isinstance(num, str) and num.isdigit()):
        if isinstance(num, str) and num.isdigit():
            num = int(num)
        format_value = lambda v: ','.join(map(tostring, v)) \
            if not v is None else '.'
    elif num=='R':
        format_value = lambda v: '.,{}'.format(
            tostring(v) if not v is None else '.')
    elif num in (1, 'A'): # Alleles are split on read-in
        format_value = lambda v: tostring(v) if not v is None else '.'
    else:
        allowed = 'A, G, . or and integer'
        raise ValueError(
            'Allowed types for "Number" INFO spec are {} while found {}'\
            .format(allowed, num))

    # For INFO function will take key and value
    # for SAMPLEs data – only the value
    if info_or_format=='info':
        return lambda k,v : '{}={}'.format(k, format_value(v))
    elif info_or_format=='format':
        return format_value


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

