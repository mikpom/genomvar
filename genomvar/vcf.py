"""Module contains classes needed for VCF parsing.

Main class is :class:`VCFReader` which is instantiated from VCF file. 
VCF can be gzipped. Bgzipping and tabix-derived indexing is also 
supported for random coordinate-based access.

:class:`VCFReader` class can iterate over rows which are tuple-like object
containing VCF field strings as attributes without conversion (except
POSition which is converted to int).

Alternatively iteration of variants is supported. It yields
:class:`genomvar.variant.GenomVariant` objects.

    >>> reader = VCFReader('example.vcf.gz')
    >>> vrt = next(reader.iter_vrt(parse_info=True, parse_samples=True))
    >>> print(vrt)
    <GenomVariant: Del chr24:23-24 G/->
    >>> print(vrt.attrib['info'])
    {'NSV': 2, 'AF': 0.5, 'DP4': (11, 22, 33, 44), 'ECNT': 1,\
    'pl': 3, 'mt': 'SUBSTITUTE', 'RECN': 18, 'STR': None}

If ``parse_info`` and ``parse_samples`` parameters are ``True`` Then
INFO and SAMPLEs fields contained in VCF are parsed and split
corresponding to an allele captured by variant object. For performance
reasons these parameters are set to False and these fields are not
parsed.

.bcf format is supported, use :class:`BCFReader`.

"""

import os
import warnings
import heapq
from itertools import dropwhile,groupby,repeat,takewhile,\
    zip_longest,islice,chain,repeat
from collections import namedtuple,OrderedDict,deque
import re
import gzip
import numpy as np
from genomvar import Reference,SINGLETON,StructuralVariantError,\
    variant,ChromSet,VCFSampleMismatch,\
    UnsortedVariantFileError,MAX_END,NoIndexFoundError
from genomvar.utils import rgn_from, grouper
from genomvar.variant import GenomVariant,VariantFactory
import genomvar
from genomvar.vcf_utils import (VCFRow, RESERVED_FORMAT_SPECS,
                                RESERVED_INFO_SPECS,
                                VCF_INFO_OR_FORMAT_SPEC,
                                VCF_INFO_OR_FORMAT_DEF_FIELDS as DATA_DEF_FIELDS,
                                _make_field_writer_func,
                                header as vcf_header,
                                dtype2string, string2dtype,
                                field_writer_simple)

class FIELDSPECS: pass

RESERVED_FORMAT = FIELDSPECS()
for spec in RESERVED_FORMAT_SPECS:
    setattr(RESERVED_FORMAT, spec[0], spec)

RESERVED_INFO = FIELDSPECS()
for spec in RESERVED_INFO_SPECS:
    setattr(RESERVED_INFO, spec[0], spec)
    

# data type of main variants array of VariantSet class
dtype0 = np.dtype([('ind',np.int_),('haplotype',np.int_),('chrom','O'),
                   ('start',np.int_),('end',np.int_),('ref','O'),('alt','O'),
                   ('vartype','O'),('start2',np.int_),('end2',np.int_)])
# data type of VCF fields array of VariantSet class
dtype1 = np.dtype([('start',np.int_),('ref',np.object_),('alt',np.object_),
                   ('id',np.object_),('qual',np.float_),('filter',np.object_),
                   ('row',np.int_)])

def _ensure_sorted(it):
    """Ensures variants are yielded in start-sorted order"""
    try:
        v = next(it)
        start,rstart,cnt,vrt = v
    except StopIteration:
        return

    hp = [(start,rstart,cnt,vrt)]
    cur_end = vrt.end
    heapq.heapify(hp)

    while True:
        try:
            start,rstart,cnt,vrt = next(it)
        except StopIteration:
            break
        if cur_end <= rstart: # Trick is that if in VCF notation
                              # start is to the right of current end
                              # then heap can be pushpopped
            sm = heapq.heappushpop(hp, (start,rstart,cnt,vrt) )
            yield sm[-1]
        else:
            heapq.heappush(hp, (start,rstart,cnt,vrt))
        cur_end = max(vrt.end,cur_end)

    for i in range(len(hp)):
        yield heapq.heappop(hp)[-1]

def _make_converter_func(tp,num,convert=True):
    """Returns a converter function based on TYPE and NUMBER arguments."""
    if convert==False:
        conv = lambda v: v
    else:
        if tp=='U1':
            conv = np.str_
        else:
            conv = tp
    if num=='R':
        def f(v,n):
            try:
                return conv(v[n+1])
            except IndexError:
                return None
            except ValueError as exc:
                if v[n+1]=='.':
                    return None
                else:
                    raise exc
    elif num=='A':
        def f(v,n):
            try:
                return conv(v[n])
            except IndexError:
                return None
            except ValueError as exc:
                if v[n]=='.':
                    return None
                else:
                    raise exc
    elif num==1:
        if issubclass(tp,np.bool_):
            f = lambda v,n: True
        else:
            def f(v,n):
                try:
                    return conv(v)
                except ValueError as exc:
                    if v=='.':
                        return None
                    else:
                        raise exc
                    
    elif num==0:
        if issubclass(tp,np.bool_):
            f = lambda v,n: True
        else:
            f = lambda v,n: None
    elif num=='G' or num=='.':
        f = lambda v,n: tuple(map(tp,v))
    else: # num 2,3 ...
        def f(v,n):
            try:
                return tuple([conv(v[i]) for i in range(num)])
            except IndexError as exc:
                r = []
                for i in range(num):
                    if i>len(v)-1:
                        r.append(None)
                    else:
                        r.append(conv(v[i]))
                return tuple(r)
    return f

no_converter = _make_converter_func(np.object_, 1, convert=False)

def _isindexed(file):
    idx = file+'.tbi'
    if os.path.isfile(idx):
        if os.path.getmtime(idx)>=os.path.getmtime(file):
            return True
        else:
            warnings.warn('Index is outdated for {}'.format(file))
    return False

def _check_VCF_order(it):
    """Raises UnsortedVariantFileError if row iterator is not sorted."""
    prev_row = namedtuple('first','CHROM POS')(object(),0)
    for row in it:
        if row.CHROM==prev_row.CHROM and row.POS<prev_row.POS:
            raise UnsortedVariantFileError
        else:
            prev_row = row
            yield row

gt_cache = {}
def _parse_gt(gt, ind):
    """Function return genotype tuple given genotype string and allele index.
    For performance results are cached based on arguments"""
    try:
        return gt_cache[(gt,ind)]
    except KeyError:
        vals = re.split('[/|]',gt)
        GT = []
        for val in vals:
            if val == '.':
                GT.append(None)
            elif int(val)==ind+1:
                GT.append(1)
            else:
                GT.append(0)
        GT = tuple(GT)
        gt_cache[(gt,ind)] = GT
        return GT

def validate_spec(spec):
    name, 
    
class VCFWriter(object):
    """Class for writing variant to VCF format."""
    def __init__(self, reference=None, info_spec=None, format_spec=None,
                 samples=None):
        if isinstance(reference, str):
            self.reference = Reference(reference)
        else:
            self.reference = reference

        self.samples = samples
        field_specs = {
            'info' : info_spec,
            'format': format_spec
        }

        self.writers = {'info':{}, 'format':{}}
        self.dtype = {'info':{}, 'format':{}}
        for info_or_format, specs in field_specs.items():
            if not specs is None:
                if isinstance(specs, list):
                    specs = specs
                elif isinstance(specs, dict):
                    specs = [{'name':f, **s} for f,s in specs.items()]
                else:
                    raise ValueError(
                        'Spec should be dict or list, got: '+str(specs))
                for spec in specs:
                    if len(spec)<3:
                        raise ValueError(
                            'INFO spec expected to have at '\
                            +'least 3 fields while {} found'\
                            .format(len(spec)))
                    if isinstance(spec, tuple):
                        _spec = VCF_INFO_OR_FORMAT_SPEC(
                            *islice(
                                chain.from_iterable(
                                    [spec, repeat(None)]),
                                6))
                    elif isinstance(spec, dict):
                        _spec = VCF_INFO_OR_FORMAT_SPEC(
                            *[spec.get(f.lower()) for f in DATA_DEF_FIELDS])
                    else:
                        raise ValueError('spec should be a tuple or dict')

                    if info_or_format=='format' and _spec.NAME=='GT':
                        self.writers[info_or_format][_spec.NAME] = \
                            lambda v: '/'.join(map(str, v)) \
                                             if not v is None else './.'
                    else:
                        self.writers[info_or_format][_spec.NAME] = \
                            _make_field_writer_func(
                                _spec.TYPE, _spec.NUMBER, info_or_format)

                    dt = {
                        'name' : _spec.NAME,
                        'number' : _spec.NUMBER,
                        'type' : _spec.TYPE,
                        'dtype' : string2dtype[_spec.TYPE],
                        'description' : _spec.DESCRIPTION,
                        'source' : _spec.SOURCE,
                        'version' : _spec.VERSION
                    }
                    self.dtype[info_or_format][_spec.NAME] = dt

    def get_header(self):
        if hasattr(self, 'dtype'):
            info = [{f.lower():p.get(f.lower(), '.') for f in DATA_DEF_FIELDS} \
                    for p in self.dtype['info'].values()]
            format = [{f.lower():p.get(f.lower(), '.') for f in DATA_DEF_FIELDS} \
                    for p in self.dtype['format'].values()]
        else:
            info = {}
            format = {}
        header = vcf_header.render(
            info=info,
            format=format,
            samples=self.samples,
            ctg_len=self.reference.ctg_len if self.reference else {})
        return header

    @staticmethod
    def _get_filter(v):
        if v is None:
            return
        elif isinstance(v, list):
            return ';'.join(v)
        else:
            return str(v)

    def get_row(self, vrt, **kwds):
        """Formats a variant to a :class:`genomvar.vcf_utils.VCFRow` instance. 

        For indels writer with reference 
        might be needeed to costruct correct REF field.

        Parameters
        ----------
        vrt : Variant instance
            variant to get row for. 

        kwds : VCF fields #FIXME remove kwds here
            optional. If given, these and only these parameters are 
            used to populate corresponding VCF fields: ``id``, 
            ``qual``, ``filter``, ``info``, ``sampdata``. 
            These parameters are taken as is and converted to string
            before returning a VCFRow. 
            If the keyword is not given corresponding field will be 
            populated from ``attrib`` attribute if possible.

            Any other keyword arguments are ignored.

        Returns
        -------
        row : VCFRow
            :class:`genomvar.vcf_utils.VCFRow` object instance

        Notes
        -----
        >>> factory = variant.VariantFactory()
        >>> writer = VCFWriter()
        >>> v1 = factory.from_edit('chr24', 2093, 'TGG', 'CCC')
        >>> row = writer.get_row(v1)
        >>> row
        <VCFRow chr24:2094 TGG->CCC>
        >>> print(row)
        chr24	2094	.	TGG	CCC	.	.	.

        >>> print(writer.get_row(v1, id=123, info='DP=10'))
        chr24	2094	123	TGG	CCC	.	.	DP=10
        """
        if hasattr(vrt, 'attrib'):
            pos, ref, alt = vrt._get_vcf_notation(
                vcf_notation=vrt.attrib.get('vcf_notation'),
                reference=self.reference)
            dt = vrt.attrib

            if 'info' in kwds:
                info = kwds['info']
            elif 'info' in dt and not dt['info'] is None:
                info = self._format_info(dt['info'])
            else:
                info = '.'

            if 'sampdata' in kwds and not kwds['sampdata'] is None:
                format, samples = kwds['sampdata'].split('\t', maxsplit=1)
            elif 'samples' in dt and dt['samples']: # TODO fix samples
                format, samples = self._format_sampdata(dt['samples'])
            else:
                if self.samples is None:
                    format = None
                    samples = None
                else:
                    format = '.'
                    samples = '.'
                
            r = VCFRow(vrt.chrom, pos,
                       kwds.get('id', dt.get('id')),
                       ref.upper(), alt.upper(),
                       kwds.get('qual', dt.get('qual')),
                       self._get_filter(kwds.get('filter', dt.get('filter'))),
                       info, format, samples)
            return r
        else:
            pos, ref, alt = vrt._get_vcf_notation(reference=self.reference)
            return VCFRow(vrt.chrom, pos, kwds.get('id', '.'),
                          ref.upper(),alt.upper(),
                          kwds.get('qual', '.'), kwds.get('filter', '.'),
                          kwds.get('info', '.'))


    def _format_info(self, info):
        if hasattr(self, 'writers'):
            _info = [w(k,info.get(k)) for k,w in \
                     self.writers['info'].items()]
        else:
            _info = [field_writer_simple(k,v) for k,v in info.items()]
        return ';'.join(_info)

    def _format_sampdata1(self, data):
        fields = self.dtype['format']
        ff = [self.writers['format'][f](None if data is None else data.get(f))\
              for f in fields]
        return ';'.join(ff)
    
    def _format_sampdata(self, sampdata):
        format = ':'.join(self.dtype['format'])
        smp = [self._format_sampdata1(sampdata.get(s)) for s in self.samples]
        return format, '\t'.join(smp)
    
class VCFReader(object):
    """Class to read VCF files."""

    def __init__(self, vcf, index=None, reference=None):
        """
        VCFReader is instantiated from VCF file. Optionally,
        index and reference arguments can be given.

        Parameters
        ----------
        index : str or bool
            By default tries ``vcf``.tbi as index file. 
            if index=True and index not found NoIndexFoundError
            is raised.
            A string with path to index can be given. *Default: None*

        reference : Reference or str
            Genome reference. *Default: None*

        Returns
        -------
        reader : VCFReader
            VCFReader object.
        """
        if isinstance(vcf,str):
            self.fl = vcf
            if self.fl.endswith('.gz') or self.fl.endswith('.bgz'):
                self.compressed = True
            else:
                self.compressed = False
            self.idx_file = check_index(self.fl, index)
            if self.compressed:
                self.openfn = gzip.open
            else:
                self.openfn = open
            self.buf = self.openfn(self.fl,'rt')
            self.opened_file = True
        else: # buffer-like object
            self.buf = vcf
            self.opened_file = False
        self._vrt_start = False
        self.reference = reference
        # Init default variant factory
        self._factory = VariantFactory(reference,normindel=False)
        self.vrt_fac = {'nonorm':self._factory}
        self._dtype = {'info':OrderedDict(),'format':OrderedDict()}
        self._samples = []
        for cnt,line in enumerate(self.buf):
            if line.startswith('##INFO'):
                dat = self._parse_dtype(line)
                self._dtype['info'][dat['name']] = dat
            elif line.startswith('##FORMAT'):
                dat = self._parse_dtype(line)
                self._dtype['format'][dat['name']] = dat
            elif line.startswith('#CHROM'):
                # this should be the last header line
                self.header_len = cnt + 1

                vals = line.strip().split('\t')
                if len(vals)>9:
                    self._samples = vals[9:]
                self._vrt_start = True
            elif not line.startswith('#'):
                break

        if self.opened_file:
            self.buf.close()
        self.sample_ind = OrderedDict()
        for ind,sample in enumerate(self._samples):
            self.sample_ind[sample] = ind

    def get_factory(self,normindel=False):
        """Returns a factory based on normindel parameter."""
        try:
            return self.vrt_fac['norm' if normindel else 'nonorm']
        except KeyError:
            vf = VariantFactory(self.reference,normindel=normindel)
            self.vrt_fac['norm' if normindel else 'nonorm'] = vf
            return vf

    @staticmethod
    def _parse_dtype(line):
        """Parses VCF header and returns correct datatypes."""
        rx = '##(?:(?:INFO)|(?:FORMAT))=\<ID=(\S+),Number=([\w.]+),'\
            +'Type=(\w+),Description="(.*)"(:?,Source="(.*)")?(?:,Version=".*")?\>'
        match = re.match(rx, line)
        NAME,NUMBER,TYPE,DESCRIPTION,SOURCE,VERSION = match.groups()

        try:
            num = int(NUMBER)
        except ValueError:
            num = NUMBER
        if isinstance(num,int):
            if num==0:
                size = 1
                tp = np.bool_
            else:
                size = num
                tp = string2dtype[TYPE]
        elif num in ('A','R'):
            size = 1
            tp = string2dtype[TYPE]
        elif num in ('.','G'):
            size = 1
            tp = np.object_
        else:
            raise ValueError('Unknown NUMBER:'+num)
        return {'name':NAME, 'type':TYPE, 'dtype':tp,
                'size':size, 'number':num,
                'description':DESCRIPTION}
        
    def _sample_indices(self,samples):
        ind = {}
        for sample in samples:
            if not sample in self._samples:
                raise VCFSampleMismatch('Sample {} not found'.format(sample))
            else:
                ind[sample] = self._samples.index(sample)
        return ind

    def find_rows(self, chrom=None, start=None, end=None,
                  rgn=None):
        """Yields rows of variant file"""
        cnt = 0
        if not rgn is None:
            chrom,start,end = rgn_from(rgn)
            for row in self._query(chrom,start,end):
                yield row
        elif not chrom is None:
            for row in self._query(chrom,start=start,end=end):
                yield row
        else:
            raise ValueError('rgn or chrom,start,end should be given')

    def iter_rows(self,check_order=False):
        """Yields rows of variant file"""
        
        buf = self.get_buf()
        buf.seek(0)
        return RowIterator(self,buf,check_order)
        
    def iter_rows_by_chrom(self,check_order=False):
        """Yields rows grouped by chromosome"""
        buf = self.get_buf()
        buf.seek(0)
        return RowByChromIterator(self,buf,check_order)

    def _parse_vrt(self, row, factory, parse_info=False, parse_samples=False,
                   parse_null=False):
        """Given the row returns variant objects with alts splitted.
        INFO, SAMPLE data are correctly parsed if needed"""

        alts = row.ALT.split(',')
        # TODO avoid list structure, generator instead
        if parse_info: # Read INFO if necessary
            info = self.dataparser.get_infod(row.INFO, alts=len(alts))
        else:
            info = repeat(None)
        if parse_samples and row.SAMPLES:
            sampdata = self.dataparser.get_sampdatad(row, len(alts))
        else:
            sampdata = repeat(None)


        for ind,(_alt,_info,_sampdata) in enumerate(zip(alts,info,sampdata)):
            try:
                base  = factory.from_edit(row.CHROM,row.POS-1,row.REF,_alt)
            except StructuralVariantError:
                warnings.warn(
                    'Structural variants are not yet supported, skipping row '\
                    +str(row.rnum))
                continue
            infod = _info
            if parse_samples and not _sampdata is None:
                sampd = dict(zip(self._samples, _sampdata))
            else:
                sampd = {}
            # FIXME avoid dictionaries
            attrib = {'info':infod,'samples':sampd,
                      'allele_num':ind,'filter':row.FILTER.split(';'),
                      'id':row.ID if row.ID!='.' else None,
                      'qual':row.QUAL}
            attrib['vcf_notation'] = {'start' : row.POS-1,'ref' : row.REF,
                                      'row' : row.rnum}
            _vrt = GenomVariant(base, attrib=attrib)
            yield _vrt

        if parse_null: # experimental
            if split_info:
                info2 = self._parse_info2(info1,allele_num='null')
            else:
                info2 = info1
            samples2 = {s:{} for s in parse_samples}
            for sm in parse_samples:
                 for key,val in samples1[sm].items():
                     if key=='GT':
                         _gt = [int(i) for i in re.split('[/|]',
                                                         samples1[sm]['GT'])]
                         samples2[sm]['GT'] = tuple([int(i==0) for i in _gt])
                     else:
                         samples2[sm][key] = val

            base  = variant.Null(row.CHROM,row.POS-1,row.REF,'-',
                                row.POS-1+len(row.REF))
            _vrt = GenomVariant(base,attrib={'info':info2,
                                             'samples':samples2,
                                             'filter':row.FILTER.split(';'),
                                             'id':row.ID,'allele_num':'null'})
            yield _vrt

    def get_records(self, parse_info=False, parse_samples=False,
                    normindel=False):
        """
        Returns parsed variant data as a dict of NumPy arrays with structured
        dtype.
        """
        vfac = self.get_factory(normindel=normindel)
        chunk_sz = 1000
        fchunks = {'vrt': {
                       f:[] for f in dtype0.names
                   },
                   'info':{
                       f:[] for f in self._dtype['info']
                   },
                   'sampdata':{
                       s:{f:[] for f in self._dtype['format']} \
                              for s in self._samples},
                   'vcf':{f:[] for f in dtype1.names}
        }
        vnum = 0

        rows = self.iter_rows()
        for row_chunk in grouper(rows,chunk_sz):
            haps = []
            tups = {'vrt':[],'info':[],'vcf':[],
                    'sampdata':{s:[] for s in self._samples}}
            for row in row_chunk:
                if not row: # reached last row
                    break
                alts = row.ALT.split(',')
                if parse_info: # Read INFO if necessary
                    info = self.dataparser.get_info(row.INFO,alts=len(alts))
                else:
                    info = repeat(None)
                if parse_samples and row.SAMPLES:
                    sampdata = self.dataparser.get_sampdata(row,len(alts))
                else:
                    sampdata = repeat(None)
                for _alt,_info,_sampdata in zip(alts,info,sampdata):
                    try:
                        base = vfac._parse_edit(
                            row.CHROM, row.POS-1,row.REF,_alt)
                    except StructuralVariantError:
                        warnings.warn(
                            'Structural variants are not yet supported, '\
                            +'skipping row '+str(row.rnum))
                        continue
                    if base[-1].is_variant_subclass(variant.AmbigIndel):
                        chrom,(start,start2),(end,end2),ref,alt,cls = base
                        tups['vrt'].append( (chrom,start,end,ref,alt,\
                                             cls,start2,end2) )
                        added = 1
                        haps.append(SINGLETON)
                    else:
                        tups['vrt'].append(base)
                        added = 1
                        haps.append(SINGLETON)

                    # Adding VCF notation related fields
                    tups['vcf'].append([row.POS-1,row.REF,row.ALT,row.ID,
                                        row.QUAL if row.QUAL!='.' else None,
                                        row.FILTER,row.rnum]*added)
                    if parse_info:
                        if added==1:
                            tups['info'].append(_info)
                        else:
                            tups.extend([_info]*added)
                    if parse_samples:
                        if added==1:
                            for sn in range(len(self._samples)):
                                tups['sampdata'][self._samples[sn]].append(_sampdata[sn])
                        else:
                            for sn,samp in enumerate(self._samples):
                                tups['sampdata'][samp]\
                                    .extend([_sampdata[sn]]*added)
                    vnum += added

            zipped = {k:list(zip_longest(*tups[k],fillvalue=0)) \
                                 for k in ('vrt','info','vcf')}
            zipped['sampdata'] = {}
            for samp in self._samples:
                zipped['sampdata'][samp] = list(zip(*tups['sampdata'][samp]))
            cur_sz = len(tups['vrt'])
            for ind,field in enumerate(dtype0.names[2:]): # skip ind and hap
                a = np.zeros(dtype=dtype0.fields.get(field)[0],
                             shape=cur_sz)
                try:
                    a[:cur_sz] = zipped['vrt'][ind]
                except IndexError as exc:
                    if field=='start2':
                        break # no ambigous indels
                    else:
                        raise exc
                fchunks['vrt'][field].append(np.resize(a,cur_sz))
            fchunks['vrt']['haplotype'].append(np.asarray(haps,dtype=np.int_))

            # Add VCF notation
            for ind,field in enumerate(dtype1.names): # skip ind and hap
                a = np.zeros(dtype=dtype1.fields.get(field)[0],
                             shape=cur_sz)
                a[:cur_sz] = zipped['vcf'][ind]
                fchunks['vcf'][field].append(np.resize(a,cur_sz))
            if parse_info:
                for ind,field in enumerate(self._dtype['info']):
                    fd = self._dtype['info'][field]
                    dtype = fd['dtype'] if fd['size']==1 \
                        else (fd['dtype'],fd['size'])
                    ar = np.zeros(dtype=dtype,shape=cur_sz)
                    try:
                        ar[:] = zipped['info'][ind]
                    except TypeError as exc:
                        if issubclass(ar.dtype.type,np.int_):
                            dtype = np.float_ if fd['size']==1 \
                                else (np.float_,fd['size'])
                            ar = np.zeros(dtype=dtype,shape=cur_sz)
                            ar[:] = zipped['info'][ind]
                        else:
                            raise exc # reraising if not integer
                    fchunks['info'][field].append(ar)

            if parse_samples:
                for samp in self._samples:
                    for ind,field in enumerate(self._dtype['format']):
                        fd = self._dtype['format'][field]
                        dtype = fd['dtype'] if fd['size']==1 \
                            else (fd['dtype'],fd['size'])
                        ar = np.zeros(dtype=dtype,shape=cur_sz)
                        try:
                            ar[:] = zipped['sampdata'][samp][ind]
                        except TypeError as exc:
                            if issubclass(ar.dtype.type,np.int_):
                                dtype = np.float_ if fd['size']==1 \
                                    else (np.float_,fd['size'])
                                ar = np.zeros(dtype=dtype,shape=cur_sz)
                                ar[:] = zipped['sampdata'][samp][ind]
                            else:
                                raise exc # reraising if not integer
                        fchunks['sampdata'][samp][field].append(ar)
                
        ret = {}
        ret['vrt'] = np.zeros(dtype=dtype0,shape=vnum)
        ret['vcf'] = np.zeros(dtype=dtype1,shape=vnum)

        for field in dtype0.names[1:]: # skipping ind
            try:
                ret['vrt'][field][:] = np.concatenate(fchunks['vrt'][field])
            except ValueError as exc:
                if field == 'start2': # No Ambigous indels
                    break
                elif len(fchunks['vrt'][field])==0: # empty VCF
                    break
                else:
                    raise exc
        for field in dtype1.names: 
            try:
                ret['vcf'][field][:] = np.concatenate(fchunks['vcf'][field])
            except ValueError as exc:
                if len(fchunks['vcf'][field])==0: # empty VCF
                    break
                else:
                    raise exc
        ret['vrt']['ind'] = np.arange(vnum)
        if parse_info:
            dt = self.dataparser.get_dtype('info')
            fields = {}
            for ind,(field,vals) in enumerate(fchunks['info'].items()):
                a = np.concatenate(vals)
                fields[field] = a
                # Need to update dtype due to possible conversion to float
                fdtype = dt[ind]
                dt[ind] = tuple([fdtype[0],a.dtype,*fdtype[2:]])
            ret['info'] = np.zeros(dtype=dt,shape=vnum)
            for field in ret['info'].dtype.names:
                ret['info'][field][:] = fields[field]
        if parse_samples:
            ret['sampdata'] = {}
            for samp in self._samples:
                dt = self.dataparser.get_dtype('format')
                fields = {}
                for ind,(field,vals) in enumerate(
                        fchunks['sampdata'][samp].items()):
                    a = np.concatenate(vals)
                    fields[field] = a
                    # Need to update dtype due to possible conversion to float
                    fdtype = dt[ind]
                    dt[ind] = tuple([fdtype[0],a.dtype,*fdtype[2:]])
                ret['sampdata'][samp] = np.zeros(dtype=dt,shape=vnum)
                for field in ret['sampdata'][samp].dtype.names:
                    ret['sampdata'][samp][field][:] = fields[field]
        return ret


    def iter_vrt(self,check_order=False,parse_info=False,
                 normindel=False,parse_samples=False):
        """
        Yields variant objects.
        
        Parameters
        ----------
        parse_info : bool
            Whether INFO fields should be parsed.  *Default: False*
        parse_samples : bool
            Whether SAMPLEs dat should be parsed.  *Default: False*
        check_order : bool
            If True will raise exception on unsorted VCF rows. *Default: False*
        normindel : bool
            If True insertions and deletions will be left normalized.
            Requires a reference on :class:`VCFReader` instantiation.

        Yields
        -------
        vrt : :class:`genomvar.variant.GenomVariant`
           Variant object 
        """

        if parse_samples==True:
            samps = 'all'
        else:
            samps = self._normalize_samples(parse_samples)
        return VrtIterator(self.iter_rows(check_order),
                           self.get_factory(normindel=normindel),
                           parse_info,samps)

    def iter_vrt_by_chrom(self,parse_info=False,
                          parse_samples=False,normindel=False,
                          check_order=False):
        """
        Generates variants grouped by chromosome
        
        Parameters
        ----------
        parse_info : bool
            Incates whether INFO fields should be parsed.  *Default: False*
        parse_samples : bool
            Incates whether SAMPLE data should be parsed.  *Default: False*
        parse_samples : bool
            If True indels will be normalized. ``VCFReader`` should have been
            instantiated with reference. *Default: False*
        check_order : bool
            If True VCF will be checked for sorting. *Default: False*
        Yields
        -------
        (chrom,it) : tuple of str and iterator
            ``it`` yields :class:`variant.Genomvariant` objects
        """
        if parse_samples==True:
            samps = 'all'
        else:
            samps = self._normalize_samples(parse_samples)
        return VrtByChromIterator(self.iter_rows(check_order),
                                  self.get_factory(normindel=normindel),
                                  parse_info,samps)

    def _normalize_samples(self,parse_samples):
        if isinstance(parse_samples,str):
            if not parse_samples in self._samples:
                raise VCFSampleMismatch('Sample {} not found'\
                                        .format(parse_samples))
            samples = [parse_samples]
        elif isinstance(parse_samples,list):
            for samp in parse_samples:
                if not samp in self._samples:
                    raise VCFSampleMismatch('Sample {} not found'\
                                            .format(samp))

            samples = parse_samples
        elif isinstance(parse_samples,bool):
            if parse_samples:
                samples = self._samples
            else:
                samples = []
        else:
            raise ValueError('parse sample should be str,list-like or bool')
        return samples

    def _variants_from_rows(self,it,factory,parse_info,parse_samples):
        cnt = 0
        for row in it:
            for vrt in self._parse_vrt(row,factory,parse_info,
                                       parse_samples):
                yield (vrt.start,row.POS-1,cnt,vrt)
                
                cnt += 1

    def find_vrt(self, chrom=None, start=0, end=MAX_END,
                 check_order=False, parse_info=False, normindel=False,
                 parse_samples=False):
        """Yields variant objects from a specified region"""
        factory = self.get_factory(normindel=normindel)
        # -1 for the case of insertions which are leftier in VCF notation
        # than they are in `genomvar` notation
        _rows = self.find_rows(
            chrom=chrom, start=start-1 if start else None, end=end)
        if check_order:
            rows = _check_VCF_order(_rows)
        else:
            rows = _rows
        if parse_samples==True:
            samps = 'all'
        else:
            samps = self._normalize_samples(parse_samples)

        _variants = self._variants_from_rows(
            rows, factory, parse_info, samps)
        for vrt in _ensure_sorted(_variants):
            if vrt.start>=end or vrt.end<=start:
                continue
            yield vrt

    def _get_tabix(self):
        if hasattr(self, 'tabix'):
            return self.tabix
        else:
            if self.idx_file:
                try:
                    self.tabix = pysam.TabixFile(
                        filename=self.fl, index=self.idx_file)
                except NameError:
                    import pysam
                    pysam.set_verbosity(0)
                    self.tabix = pysam.TabixFile(
                        filename=self.fl, index=self.idx_file)
                except OSError as exc:
                    msg = exc.args[0]
                    if 'index' in msg and 'not found' in msg:
                        raise NoIndexFoundError
                    else:
                        raise exc
            else:
                raise NoIndexFoundError('No idx file')
        return self.tabix
            
    def _query(self,chrom,start=None,end=None):
        """User tabix index to fetch VCFRow's"""
        if start and not end:
            raise ValueError('"start" was given but "end" was not')
        tabix = self._get_tabix()
        try:
            lines = tabix.fetch(chrom,start,end)
        except ValueError: # pysam raises on wrong chromosome
            return
        for line in lines:
            yield VCFRow(*line.strip().split('\t',maxsplit=9))

    def get_chroms(self,allow_no_index=False):
        """Returns ``ChromSet`` corresponding to VCF. If indexed 
        then index is used for faster access. Alternatively if ``allow_no_index``
        is True the whole file is parsed to get chromosome ordering."""
        try:
            tabix = self._get_tabix()
            self._chroms = ChromSet(tabix.contigs)
            return self._chroms
        except NoIndexFoundError:
            if allow_no_index:
                self._chroms = ChromSet()
                is_comment = lambda l: l.startswith('#')
                with self.openfn(self.fl,'rt') as fh:
                    lines = dropwhile(is_comment,fh)
                    for line in lines:
                        self._chroms.add(line.split('\t',maxsplit=1)[0])
                return self._chroms
            else:
                msg = 'Need index to get chroms'
                raise NotImplementedError(msg)

    def get_buf(self):
        if self.opened_file:
            r = self.openfn(self.fl,'rt')
            return r
        else:
            return self.buf
        
    def close(self):
        if self.opened_file:
            self.buf.close()

    @property
    def samples(self):
        return list(self.sample_ind)
    
    @property
    def chroms(self):
        if hasattr(self,'_chroms'):
            return self._chroms
        else:
            self.get_chroms()
            return self._chroms

    @chroms.setter
    def chroms(self,value):
        raise NotImplementedError

    @chroms.deleter
    def chroms(self):
        delattr(self,'_chroms')

    @property
    def dataparser(self):
        try:
            return self._dataparser
        except AttributeError:
            self._dataparser = _DataParser(self._dtype,self.sample_ind)
            return self._dataparser

    @dataparser.setter
    def dataparser(self,value):
        raise NotImplementedError

    @dataparser.deleter
    def dataparser(self):
        delattr(self,'_dataparser')
        
class BCFReader(VCFReader):
    def __init__(self,bcf,index=False,reference=None):
        # TODO leverage inheritance
        self.fl = bcf
        self.idx_file = None
        if isinstance(index,bool):
            if index:
                idx = self.fl+'.csi'
                if os.path.isfile(idx):
                    self.idx_file = idx
                else:
                    raise OSError('Index not found')
        elif isinstance(index,str):
            if os.path.isfile(index):
                self.idx_file = index
            else:
                raise OSError('{} not found'.format(index))

        self.reference = reference
        # Init default variant factory
        self._factory = VariantFactory(reference,normindel=False)
        self.vrt_fac = {'nonorm':self._factory}
        self._dtype = {'info':OrderedDict(),'format':OrderedDict()}
        buf = self.get_buf()
        for cnt,line in enumerate(str(buf.header).splitlines()):
            if line.startswith('##INFO'):
                dat = self._parse_dtype(line)
                self._dtype['info'][dat['name']] = dat
            elif line.startswith('##FORMAT'):
                dat = self._parse_dtype(line)
                self._dtype['format'][dat['name']] = dat
            elif line.startswith('#CHROM'):
                # this should be the last header line
                self.header_len = cnt + 1

                vals = line.strip().split('\t')
                if len(vals)>9:
                    self._samples = vals[9:]
                else:
                    self._samples = []
                self._vrt_start = True
                break

        self.sample_ind = OrderedDict()
        for ind,sample in enumerate(self._samples):
            self.sample_ind[sample] = ind

    def iter_rows(self, check_order=None):
        return BCFRowIterator(self, self.get_buf(),
                              check_order=check_order)

    def get_buf(self):
        try:
            f = pysam.VariantFile(self.fl)
        except NameError:
            import pysam
            pysam.set_verbosity(0)
            f = pysam.VariantFile(self.fl)
        return f

    @property
    def dataparser(self):
        try:
            return self._dataparser
        except AttributeError:
            self._dataparser = _BCFDataParser(self._dtype,self.sample_ind)
            return self._dataparser

    @dataparser.setter
    def dataparser(self,value):
        raise NotImplementedError

    @dataparser.deleter
    def dataparser(self):
        delattr(self,'_dataparser')

class _DataParser(object):
    """Object for parsing INFO and SAMPLES data."""
    def __init__(self, dtype, sample_ind):
        self.dtype = dtype
        self.converters = {'info':{},'format':{}}
        self.none = {'info':[],'format':[]}
        for info_or_format in ['info','format']:
            for field,props in self.dtype[info_or_format].items():
                tp = props['dtype']
                sz = props['size']
                num = props['number']
                if info_or_format=='format' and field=='GT':
                    self.converters[info_or_format][field] = _parse_gt
                else:
                    self.converters[info_or_format][field] \
                              = _make_converter_func(tp,num)
                self.none[info_or_format].append(None if sz<=1 else [None]*sz)
        # self.converters['format']['GT'] = \
        #     lambda v,a: _parse_gt(v,a)
        self.order = {
            'info':{f:ind for ind,f in enumerate(self.dtype['info'])},
            'format':{f:ind for ind,f in enumerate(self.dtype['format'])}
        }
        self.sample_ind = sample_ind

    def tokenize_info(self,INFO):
        """Function splits info on simple key value pairs"""
        info = []
        for field in INFO.split(';'):
            if not field or field=='.':
                continue
            keyval = field.split('=',maxsplit=1)
            if len(keyval)==1:
                info.append((field, None))
            else:
                key,val=keyval
                try:
                    num = self.dtype['info'][key]['number']
                except KeyError:
                    num=1
                if num in [0,1]:
                    info.append( (key, val) )
                else:
                    info.append( (key, val.split(',')) )
        return info
        
    def tokenize_sampdata(self,FORMAT,SAMPLES):
        def _maybe_split(key,val):
            if key in self.dtype['format'] and \
                       not self.dtype['format'][key]['number'] in [0,1]:
                return val.split(',')
            else:
                return val

        fmt = FORMAT.split(':')
        _samples = SAMPLES.split('\t')
        sampdata = [[(k,_maybe_split(k,v)) for k,v in zip(fmt,d.split(':'))]\
                     for d in _samples]
        return sampdata

    def get_info(self,INFO,alts=1,parse_null=False):
        """Given INFO string and number of alt alleles returns a list
        of lists with data corresponding to alleles, then fields."""
        tokenized = self.tokenize_info(INFO)
        info2 = []
        for an in range(alts):
            info1 = list(self.none['info'])
            for ind,(key,val) in enumerate(tokenized):
                try:
                    conv = self.converters['info'][key]
                except KeyError:
                    conv = no_converter

                try:
                    v = conv(val,an)
                except ValueError as exc:
                    if val=='' or val=='.':
                        v = None
                    else:
                        raise exc
                info1[self.order['info'][key]] = v
            info2.append(info1)

        return info2

    def get_infod(self, INFO, alts=1, parse_null=False):
        """Given INFO string and number of alt alleles returns a list
        of dicts with data corresponding to alleles"""
        tokenized = self.tokenize_info(INFO)
        info2 = []
        for an in range(alts):
            info1 = {}
            for ind,(key,val) in enumerate(tokenized):
                try:
                    conv = self.converters['info'][key]
                except KeyError:
                    conv = no_converter
                try:
                    v = conv(val,an)
                except ValueError as exc:
                    if val=='' or val=='.':
                        v = None
                    else:
                        raise exc
                info1[key] = v
            info2.append(info1)
        return info2

    def get_sampdatad(self, row, alts=1):
        """Given row and number of alt alleles returns a list
        of dicts with data corresponding to alleles"""
        tokenized = self.tokenize_sampdata(row.FORMAT, row.SAMPLES)
        sampdata2 = [[None \
                       for j in range(len(self.sample_ind))] \
                           for i in range(alts)]
        for an in range(alts):
            info1 = {}
            for sample, sn in self.sample_ind.items():
                _sampdata1 = {}
                for key, val in tokenized[sn]:
                    try:
                        conv = self.converters['format'][key]
                    except KeyError:
                        conv = no_converter
                    try:
                        v = conv(val,an)
                    except ValueError as exc:
                        if val=='' or val=='.':
                            v = None
                        else:
                            raise exc
                    _sampdata1[key] = v
                sampdata2[an][sn] = _sampdata1
        return sampdata2

    def get_sampdata(self,row,alts=1):
        """Given SAMPLE data string and number of alt alleles returns a list
        of lists of lists with data corresponding to alleles, then samples,
        then fields."""
        sampdata = [[list(self.none['format']) \
                       for j in range(len(self.sample_ind))] \
                           for i in range(alts)]
        _samples = self.tokenize_sampdata(row.FORMAT, row.SAMPLES)
        
        for an in range(alts):
            for sample,sn in self.sample_ind.items():
                _sampdata = list(self.none['format'])
                for key,val in _samples[sn]:
                    order = self.order['format'][key]
                    v = self.converters['format'][key](val,an)
                    _sampdata[order] = v
                sampdata[an][sn] = _sampdata
        return sampdata

    def get_dtype(self,info_or_format):
        if not info_or_format in ('format','info'):
            raise ValueError
        dtype = []
        for field,props in self.dtype[info_or_format].items():
            if props['size']==1:
                dtype.append( (field,props['dtype']) )
            else:
                dtype.append( (field,props['dtype'],props['size']) )
        return dtype

class _BCFDataParser(_DataParser):
    def __init__(self,dtype,sample_ind):
        self.dtype = dtype
        self.converters = {'info':{},'format':{}}
        self.none = {'info':[],'format':[]}
        for info_or_format in ['info','format']:
            for field,props in self.dtype[info_or_format].items():
                tp = props['dtype']
                sz = props['size']
                num = props['number']
                self.converters[info_or_format][field] = _make_converter_func(tp,num,convert=False)
                self.none[info_or_format].append(None if sz<=1 else [None]*sz)
        self.converters['format']['GT'] = \
            lambda v,a: tuple([e==a for e in v])
        self.order = {
            'info':{f:ind for ind,f in enumerate(self.dtype['info'])},
            'format':{f:ind for ind,f in enumerate(self.dtype['format'])}
        }
        self.sample_ind = sample_ind

    def tokenize_info(self, INFO):
        return INFO.items()

    def tokenize_sampdata(self, FORMAT, SAMPLES):
        return [v.items() for v in SAMPLES.values()]
                    
class RowIterator:
    def iterate(self):
        cnt = 0
        for line in dropwhile(lambda l: l.startswith('#'),self.fh):
            row = VCFRow(*line.strip().split('\t',maxsplit=9),
                         rnum=cnt)
            yield row
            cnt += 1
        self.close()

    def __init__(self,reader,fh,check_order):
        self.fh = fh
        self.check_order = check_order
        self.reader = reader
        if self.check_order:
            self.iterator = _check_VCF_order(self.iterate())
        else:
            self.iterator = self.iterate()

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)

    def close(self):
        self.fh.close()

class RowByChromIterator(RowIterator):
    def iterate(self):
        rows = super().iterate()
        if self.check_order:
            rows = _check_VCF_order(rows)
        for chrom,it in groupby(rows,key=lambda r: r.CHROM):
            yield chrom,it
        
    def __init__(self,*args,**kwds):
        super().__init__(*args,**kwds)
        
class BCFRowIterator(RowIterator):
    def __init__(self, *args, **kwds):
        super().__init__(*args, *kwds)

    def iterate(self):
        cnt = 0
        for rec in self.fh:
            row = VCFRow(
                rec.contig,rec.pos,rec.id,
                rec.ref, ','.join(rec.alts), rec.qual, ','.join(list(rec.filter)),
                rec.info, rec.format, rec.samples,
                rnum=cnt)
            yield row
            cnt += 1
        self.close()

class VrtIterator():
    def iterate(self):
        variants = self.reader._variants_from_rows(
                self.rows,self.factory,self.parse_info,self.samps)
        for vrt in variants:
            yield vrt
        self.close()
        
    def __init__(self,row_iterator,factory,parse_info,samps,
                 ensure_sorted=True):
        try:
            self.fh = row_iterator.fh
        except AttributeError:
            self.fh = None
        self.rows = row_iterator
        self.reader = row_iterator.reader
        self.factory = factory
        self.parse_info = parse_info
        self.samps = samps
        if ensure_sorted:
            self.iterator = _ensure_sorted(self.iterate())
        else:
            self.iterator = self.iterate()
        
    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)

    def close(self):
        if self.fh is None:
            return
        self.fh.close() 
    
class VrtByChromIterator(VrtIterator):
    def iterate(self):
        vrt = super().iterate()
        for chrom,it in groupby(vrt,key=lambda r: r[-1].chrom):
            yield chrom,_ensure_sorted(it)

    def __init__(self, *args, **kwds):
        kwds['ensure_sorted'] = False
        super().__init__(*args, **kwds)

        
def _get_reader(file, index=False, reference=None):
    if isinstance(file, str) and file.endswith('.bcf'):
        return BCFReader(file, reference=reference,
                         index=index)
    else:
        return VCFReader(file, reference=reference,
                         index=index)

def check_index(file, index):
    if index is True or index is None:
        idx = file+'.tbi'
        if os.path.isfile(idx):
            return idx
        elif index is True:
            raise NoIndexFoundError('no index found for '+str(file))
        else:
            return None
    elif isinstance(index, str):
        if os.path.isfile(index):
            return index
        else:
            raise ValueError('{} not found'.format(index))
    
