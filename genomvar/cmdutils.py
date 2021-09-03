import itertools
import warnings
from genomvar.vcf import VCFReader, VCFWriter,\
    dtype2string,_isindexed
from genomvar.vcf_utils import row_tmpl, header as vcf_header
from genomvar.varset import VariantSetFromFile,VariantSetFromFile
from genomvar import variant

def _same_order(chroms1,chroms2):
    common_chroms = set(chroms1).intersection(chroms2)
    _chroms1 = list(filter(lambda e: e in common_chroms,chroms1))
    _chroms2 = list(filter(lambda e: e in common_chroms,chroms2))
    if _chroms1==_chroms2:
        return True
    else:
        return False

def _cmp_vcf(f1,f2,out,match_partial=False,chunk=1000):
    """
    Writes comparison of two VCF files to a specified file handle.
    """
    info = [('VT', 1, 'String', 'Variant type'),
            ('whichVCF', 1, 'String', 'Which input VCF contains the variant; first, second or both'),
            ('ln', 1, 'Integer', 'Line number in input VCF variant originating from'),
            ('ln2', '.', 'Integer', 'If whichVCF is both indicates line numberin the second file')]

    writer = VCFWriter(info_spec=info)
    header = writer.get_header()
    out.write(header)

    if _isindexed(f1):
        vs1 = VariantSetFromFile(f1)
    else:
        warnings.warn('{} not indexed; may impact performance.'.format(f1))
        vs1 = VariantSetFromFile(f1)
    if _isindexed(f2):
        vs2 = VariantSetFromFile(f2)
    else:
        warnings.warn('{} not indexed; may impact performance.'.format(f2))
        vs2 = VariantSetFromFile(f2)

    _which = {0:'first',1:'second',2:'both'}
    nof_vrt = {i:0 for i in _which}
    cb = lambda m: [v.attrib['vcf_notation']['row'] for v in m]
    for which,vrt in vs1._cmp_vrt(vs2,action='all')\
                        .iter_vrt(callback=cb):
        
        nof_vrt[which] += vrt.nof_unit_vrt()
        if which==0:
            lineno = vrt.attrib['vcf_notation']['row']+vs1.vcfreader.header_len+1
        elif which==1:
            lineno = vrt.attrib['vcf_notation']['row']+vs2.vcfreader.header_len+1
        if which==2:
            lineno = vrt.attrib['vcf_notation']['row']+vs1.vcfreader.header_len+1
            lineno2 = [vs2.vcfreader.header_len+n+1 for n in vrt.attrib['cmp']]

        vrt.attrib['info'] = {'whichVCF':_which[which],
                              'ln':lineno,
                              'ln2':lineno2 if which==2\
                                  else None}
        try:
            row = writer.get_row(vrt)
        except ValueError as exc:
            if vrt.is_variant_instance(variant.Haplotype) \
                    or vrt.is_variant_instance(variant.Asterisk):
                continue
            else:
                raise exc

        try:
            out.write(str(row)+'\n')
        except BrokenPipeError:
            exit(1)
    return nof_vrt
