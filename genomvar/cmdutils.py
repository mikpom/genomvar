import itertools
import warnings
from genomvar import UnsortedVariantFileError,\
    DifferentlySortedChromsError,Reference
from genomvar.vcf import VCFReader,header_simple,header_info,\
    row_tmpl,dtype2string,_vcf_row,_isindexed
from genomvar.varset import VariantFileSet,IndexedVariantFileSet

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
    info = [{'name':'VT','number':1,'type':'String','description':
               '"Variant type"'},
            {'name':'whichVCF','number':1,'type':'String','description':
               '"Which input VCF contains the variant; first, second or both"'},
            {'name':'ln','number':1,'type':'Integer','description':
               '"Line number in input VCF variant originating from"'},
            {'name':'ln2','number':'.','type':'Integer','description':
               '"If whichVCF is both indicates line number'\
                +'in the second file"'}]
    header = header_info.render(ctg_len={},info=info)
    out.write(header)

    if _isindexed(f1):
        vs1 = IndexedVariantFileSet(f1)
    else:
        warnings.warn('{} not indexed; may impact performance.'.format(f1))
        vs1 = VariantFileSet(f1)
    if _isindexed(f2):
        vs2 = IndexedVariantFileSet(f2)
    else:
        warnings.warn('{} not indexed; may impact performance.'.format(f2))
        vs2 = VariantFileSet(f2)

    _which = {0:'first',1:'second',2:'both'}
    nof_vrt = {i:0 for i in _which}
    cb = lambda m: [v.attrib['vcf_notation']['row'] for v in m]
    for which,vrt in vs1._cmp_vrt(vs2,action='all')\
                        .iter_vrt(callback=cb):
        
        nof_vrt[which] += vrt.nof_unit_vrt()
        if which==0:
            lineno = vrt.attrib['vcf_notation']['row']+vs1._reader.header_len+1
        elif which==1:
            lineno = vrt.attrib['vcf_notation']['row']+vs2._reader.header_len+1
        if which==2:
            lineno = vrt.attrib['vcf_notation']['row']+vs1._reader.header_len+1
            lineno2 = [vs2._reader.header_len+n+1 for n in vrt.attrib['cmp']]

        vrt.attrib['info'] = {'whichVCF':_which[which],
                              'ln':lineno,
                              'ln2':','.join(map(str,lineno2)) if which==2\
                                  else '.'}
            
        row = _vcf_row(vrt,template=row_tmpl)
        out.write(row+'\n')
    return nof_vrt
