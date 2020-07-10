"""
Module ``varset`` defines variant set classes which are containers of variants
supporting set-like operations, searching variants, export to NumPy and some
other goodies.

There are several classes depending on whether all the variants should be loaded
into memory or not and whether the object should be mutable:

    - :class:`~genomvar.varset.VariantSet` is immutable, in-memory, random
      access supporting class

    - :class:`~genomvar.varset.MutableVariantSet` mutable, in-memory, random
      access

    - :class:`~genomvar.varset.VariantFileSet` immutable, VCF file storage, no
      random access

    - :class:`~genomvar.varset.IndexedVariantFileSet` immutable, VCF file
      storage, random access
"""
import warnings
from copy import copy
import itertools
import collections
import numpy as np
from numpy.lib import recfunctions as rfn
from bx.intervals.intersection import Interval as bxiv
from bx.intervals.intersection import IntervalTree
from genomvar import variant, variantops as ops, MAX_END
from genomvar.variant import VariantBase,AmbigIndel,\
    GenomVariant,Haplotype,VariantFactory
from genomvar.utils import rgn_from,zip_variants,\
    chunkit,no_ovlp
from genomvar.vcf import VCFReader,header_simple,row_tmpl,_vcf_row,\
    dtype0,dtype1
from genomvar import OverlappingHaplotypeVars,\
    Reference,DifferentlySortedChromsError,\
    singleton,ChromSet

class VariantSetBase(object):
    """Base class for variant set subclasses.

    Attributes
    ----------
    reference: :class:`genomvar.Reference`
        Object representing reference
    chroms: :class:`genomvar.ChromSet`
        list-like representation of chromosomes in the set
    ctg_len: dict
        dict mapping chromosome names to corresponding lengths.
        
        This keeps track of the rightmost position of variants in the set
        if reference was not provided. If reference was given ``ctg_len`` 
        is taken from length of chromosomes in reference.
    """

    def __init__(self,variants=None,data=None,reference=None):
        if not reference is None:
            if isinstance(reference,Reference):
                self.reference=reference
            elif isinstance(reference,str):
                self.reference = Reference(reference,cache_dst=100000)
            else:
                raise TypeError('reference not understood')
            self.ctg_len = self.reference.ctg_len
            self._chroms = ChromSet(self.reference.get_chroms())
        else:
            self.reference = None
            self.ctg_len = {}
            self._chroms = ChromSet()
        
        self._data = data if data else {}
        self._factory = VariantFactory(reference=self.reference,
                                       normindel=False)
        self._fac = {'nonorm':self._factory}
        
    def get_chroms(self):
        return self.chroms
    
    def _find_iv(self,chrom,start=0,end=MAX_END):
        if chrom not in self._data:
            return iter([])
        for iv in self._data[chrom].find(start,end):
            if iv.value in self._obs:
                continue
            yield iv

    def find_vrt(self,chrom=None,start=0,end=MAX_END,rgn=None,
                 expand=False,check_order=False):
        """Finds variants in specified region.

        If not ``chrom`` nor ``rgn`` parameters are given it 
        traverses all variants.

        Parameters
        ----------

        chrom: str, optional
            Chromosome to search variants on. If not given, all chroms
            are traversed and ``start`` and ``end`` parameters ignored.
        start: int, optional
            Start of search region (0-based). *Default: 0*
        end: int, optional
            End of search region (position ``end`` is excluded). *Default: 
            Chromosome length*
        rgn: str, optional
            If given chrom, start and end parameters are
            ignored.  A string specifying search region, e.g.
            ``chr15:365,008-365,958``.
        expand: bool, optional
            Indicates whether to traverse haplotypes. *Default: False*

        Returns
        -------
        List of :class:`genomvar.variant.VariantBase` subclass objects
        """
        # here only need to parse arguments correctly and delegate
        # actual job to _find_vrt
        if not rgn is None:
            chrom,start,end = rgn_from(rgn)
            for vrt in self._find_vrt(chrom,start,end,
                                       expand=expand):
                yield vrt
        elif not chrom is None:
            for vrt in self._find_vrt(chrom,start=start,end=end,
                                       expand=expand):
                yield vrt
        else:
            for b in self.iter_vrt(expand=expand):
                yield b

    def iter_vrt_by_chrom(self,check_order=False):
        """Iterates variants grouped by chromosome."""
        for chrom in self.chroms:
            yield (chrom,self.find_vrt(chrom=chrom))

    def iter_vrt(self,expand=False):
        """
        Iterate over all variants in the set.
        
        Parameters
        ----------
        expand : bool, optional
            If ``True`` variants in encountered haplotypes are yielded.
            *Default: False*
        Yields
        -------
        :class:`genomvar.variant.GenomVariant`
        """
        for chrom in sorted(self.chroms):
            for vrt in self._find_vrt(chrom,start=0,end=MAX_END,
                                       expand=expand):
                yield vrt

    def ovlp(self,vrt0,match_ambig=False):
        """
        Returns variants altering the same positions as ``vrt0``.
        
        Parameters
        ----------

        vrt0 : VariantBase-like
            Variant to search overlap with.
        match_ambig : bool
            If False ambigous indels treated as regular indels. *Default: False*

        Returns
        -------
        list of variants
        """
        chrom = vrt0.chrom
        ovlp = []
        if not match_ambig and vrt0.is_instance(variant.AmbigIndel):
            start = vrt0.act_start
            end = vrt0.act_end
        else:
            start = vrt0.start
            end = vrt0.end

        # Insertions need special treatment due to arbitrariness of
        # position definition
        if vrt0.is_instance(variant.Ins):
            for _vrt in self._find_vrt(chrom,start,end):
                if _vrt.is_instance(variant.AmbigIndel) and not match_ambig:
                    _vrt_end = _vrt.act_end
                    _vrt_start = _vrt.act_start
                else:
                    _vrt_end = _vrt.end
                    _vrt_start = _vrt.start
                # insertions in the same spot do overlap
                # if _vrt.is_instance(variant.Ins) and _vrt_start==start:
                #     ovlp.append(_vrt)
                #     continue
                # next is no overlap case
                if max(_vrt_start,start)>=min(_vrt_end,end):
                    continue
                # next is for case like
                # CT|TTCTTG
                #   T
                #    G
                # but if it is insertion than can add to overlap
                # due to identical position
                elif end<=_vrt_start+1 and not _vrt.is_instance(variant.Ins):
                    continue
                else:
                    ovlp.append(_vrt)
        else:
            ovlp += self._ovlp1(chrom,start,end,
                                match_ambig=match_ambig)

        return ovlp

    def _ovlp1(self,chrom,start,end,match_ambig=False):
        """Returns variants in the specified region. Locations of 
        indels are taken care of."""
        ovlp = []
        
        for _vrt in self._find_vrt(chrom,start,end):
            if _vrt.is_instance(variant.AmbigIndel) and not match_ambig:
                _vrt_end = _vrt.act_end
                _vrt_start = _vrt.act_start
            else:
                _vrt_end = _vrt.end
                _vrt_start = _vrt.start
            # if _vrt is insertion and the search region is further
            if _vrt.is_instance(variant.Ins) and _vrt_end<=start+1:
                continue
            # check maybe no overlp after ambiguity removal
            if max(_vrt_start,start)>=min(_vrt_end,end):
                continue
            else:
                ovlp.append(_vrt)
        return ovlp

    def match(self,vrt,match_partial=True,match_ambig=False):
        """
        Find variants mathing vrt.
        
        This method aims to answer common question whether variant ``vrt``
        is present in the set.

        Parameters
        ----------
        vrt : VariantBase-like
            Variant to search match for. 
        match_partial: bool, optional
            If ``False`` common positions won't be searched in non-identical
            MNPs (relevant if ``vrt`` is MNP or SNP or haplotype, containing
            these). *Default: True*
        match_ambig: bool, optional
            If ``False`` ambigous indels will be treated as regular indels
            (relevant if ``vrt`` is Ambigindel or haplotype, containing ones).
            *Default: False*
        Returns
        -------
        list or dict with matching variants
            dict is returned only of ``vrt`` is a Haplotype
        """
        ovlp = list(self._find_vrt(vrt.chrom,vrt.start,vrt.end))
        return ops.matchv(vrt,ovlp,match_ambig=match_ambig,
                          match_partial=match_partial)

    def diff_vrt(self,other,match_partial=True,match_ambig=False):
        """
        Generate variants different between ``self`` and ``other``.
        
        This creates a lazily evaluated object of class :class:`CmpSet`. This
        object can be further inspected with ``.iter_vrt()`` or ``.region()``

        Parameters
        ----------
        other : variant set
            Variant set to compare against ``self``

        match_partial: bool, optional
            If ``False`` common positions won't be searched in non-identical
            MNPs. *Default: True*

        match_ambig: bool, optional
            If ``False`` ambigous indels will be treated as regular indels.
            *Default: False*
        Returns
        -------
        :class:`CmpSet`
            Transient comparison set
        """
        return self._cmp_vrt(other,action='diff',match_ambig=match_ambig,
                          match_partial=match_partial)

    def comm_vrt(self,other,match_partial=True,match_ambig=False):
        """
        Generate variants common between ``self`` and ``other``.
        
        This creates a lazily evaluated object of class :class:`CmpSet`. This
        object can be further inspected with ``.iter_vrt()`` or ``.region()``

        Parameters
        ----------
        other : variant set
            Variant set to compare against ``self``
        match_partial: bool, optional
            If ``False`` common positions won't be searched in non-identical
            MNPs. *Default: True*
        match_ambig: bool, optional
            If ``False`` ambigous indels will be treated as regular indels.
            *Default: False*
        Returns
        -------
        :class:`CmpSet`
            Transient comparison set
        """
        return self._cmp_vrt(other,action='comm',match_ambig=match_ambig,
                          match_partial=match_partial)


    def _cmp_vrt(self,other,action,match_partial=True,match_ambig=False):
        return CmpSet(self,other,action,match_partial=match_partial,
                      match_ambig=match_ambig)
    
    def to_vcf(self,fh):
        """Writes a minimal VCF with variants from the set. 

        INFO and SAMPLE data from source data is not preserved.

        Parameters
        ----------

        fh : handle or str
           Where to write variants.

        Returns
        -------
        None
        """
        header = header_simple.render(ctg_len=self.ctg_len)
        fh.write(header)
        for vrt in self.find_vrt(expand=True):
            if vrt.is_instance(variant.Haplotype):
                continue
            row = _vcf_row(vrt,template=row_tmpl,
                               reference=self.reference)
            fh.write(row+'\n')

    def get_factory(self,normindel=False):
        """Returns a VariantFactory object. It will inherit the same reference
        as was used for instantiation of ``self``"""
        try:
            return self._fac['norm' if normindel else 'nonorm']
        except KeyError:
            vf = VariantFactory(self.reference,normindel=normindel)
            self._fac['norm' if normindel else 'nonorm'] = vf
            return vf

class VariantSet(VariantSetBase):
    """
    Immutable class for storing variants in-memory.

    Supports random access by variant location.
    Can be instantiated from a VCF file (see
    :meth:`VariantSet.from_vcf`) or from iterable
    containing variants (see :meth:`VariantSet.from_variants`).

    Attributes
    ----------
    reference: Reference
        Genome reference.
    ctg_len: dict
        dict mapping Chromosome name to corresponding lengths.
        
        This keeps track of the rightmost position of variants in the set
        if reference was not provided. If reference was given ``ctg_len`` 
        is taken from length of chromosomes in the reference.
    """
    def __init__(self,variants,vcf_notation=None,info=None,
                 sampdata=None,reference=None):
        super().__init__(reference=reference)
        self._variants = variants
        self._info = info
        self._sampdata = sampdata
        self._vcf_notation = vcf_notation
        for ind,(chrom,start,end) in enumerate(zip(self._variants['chrom'],
                            self._variants['start'],self._variants['end'])):
            self._data.setdefault(chrom,IntervalTree())\
                     .add_interval(bxiv(start,end,ind))
            if not reference:
                self.ctg_len[chrom] = max(self.ctg_len.get(chrom,0),end)
                self.chroms.add(chrom)
            else:
                try:
                    if end > self.ctg_len[chrom]:
                        msg = 'Variant at {}:{}-{} is out of reference range'\
                            .format(chrom,start,end)
                        raise ValueError(msg)
                except KeyError as exc:
                    if not chrom in self.chroms:
                        msg = 'Chromosome {} not in reference'.format(chrom)
                        raise ValueError(msg)

    @staticmethod
    def _get_vrt(row):
        """Constructs a variant from row"""
        
        ind,hap,chrom,start,end,ref,alt,cls,start2,end2 = row
        try:
            return cls(chrom=chrom,start=start,end=end,ref=ref,alt=alt)
        except TypeError:
            return cls(chrom=chrom,start=(start,start2),
                       end=(end,end2),ref=ref,alt=alt)

    @staticmethod
    def _get_tups(vrt_iter):
        """Constructs tuples from variants."""
        cnt = 0
        for vrt in vrt_iter:
            if vrt.is_instance(variant.Haplotype):
                ind = cnt
                yield (ind,singleton,*vrt.tolist())
                cnt += 1
                for child in vrt.variants:
                    a = child.tolist()
                    yield (cnt,ind,*a)
                    cnt += 1
            else:
                a = vrt.tolist()
                b = (cnt,singleton,*a)
                yield b
                cnt += 1

    @staticmethod
    def _rows_from_variants(variants):
        """Returns rows given variants"""
        cnt = 0
        for vrt in variants:
            if vrt.vtp.is_subclass(variant.Haplotype):
                cnt += len(list(vrt._variants)) + 1
            else:
                cnt += 1
        _variants = np.zeros(dtype=dtype0,shape=cnt)
        zipped = list(itertools.zip_longest(*VariantSet._get_tups(variants),
                                            fillvalue=0))
        for ind,field in enumerate(dtype0.fields):
            try:
                _variants[field] = zipped[ind]
            except IndexError:
                continue
        return _variants

    def _find_vrt(self,chrom,start,end,expand=False):
        """Uses IntervalTree to find variants and yields them. 
        Implementation is class-specific."""
        haplotypes = set()
        start = start if start else 0
        end = end if end else MAX_END
        try:
            intervals = self._data[chrom].find(start,end)
        except KeyError:
            return
        for iv in intervals:
            row = self._variants[iv.value].tolist()
            vartype=row[7]
            if vartype==variant.Haplotype:
                haplotypes.add(row[-1])
                if not self._info is None:
                    attrib = {'info':self._info[iv.value]}
                else:
                    attrib = {}
                if not self._sampdata is None:
                    attrib.update({'samples':{s:self._sampdata[s][iv.value]\
                                              for s in self._sampdata}})
                # TODO find haplotypes
                base = self._factory.from_data(row[1:7])
                yield GenomVariant(base=base,attrib=attrib)
            else:
                if row[1]!=singleton and not expand: # part of a haplotype
                    continue
                if not self._info is None:
                    attrib = {'info':self._info[iv.value]}
                else:
                    attrib = {}
                if not self._vcf_notation is None:
                    attrib['vcf_notation'] = dict(zip(dtype1.names,
                                    self._vcf_notation[iv.value].tolist()))
                if not self._sampdata is None:
                    attrib.update({'samples':{s:self._sampdata[s][iv.value]\
                                              for s in self._sampdata}})
                base = VariantSet._get_vrt(self._variants[iv.value].tolist())
                yield GenomVariant(base=base,attrib=attrib)

    def copy(self):
        """Returns a copy of variant set"""
        cpy = self.__new__(self.__class__)
        cpy.__init__(variants=self._variants.copy(),
                     info=self._info.copy(),
                     sampdata={s:self._sampdata[s].copy() \
                               for s in self._sampdata} )
        return cpy
    def nof_unit_vrt(self):
        """Returns number of simple genomic variations in a variant set. 

        SNPs and indels count as 1. MNPs as number of affected nucleotides.
        """
        N = 0
        for start,end,cls in zip(self._variants['start'],self._variants['end'],
                                 self._variants['vartype']):
            if cls.is_subclass(variant.Haplotype):
                continue
            else:
                if cls.is_subclass(variant.Indel):
                    N += 1
                elif cls.is_subclass(variant.MNP):
                    N += (end-start)
        return N

    @classmethod
    def from_vcf(cls,vcf,reference=None,parse_info=False, parse_samples=False,
                 normindel=False, parse_null=False):
        """
        Parse VCF variant file and return VariantSet object
        
        Parameters
        ----------
        vcf : str
            VCF file to read data from.
        parse_info: bool, optional
            If ``True`` INFO fields will be parsed.
            *Default: False*.
        parse_samples: bool, str or list of str, optional
            If ``True`` all sample data will be parsed.  If string, sample 
            with this name will be parsed.  If list of strings only these
            samples will be parsed. *Default: False*
        normindel: bool, optional
            If ``True`` indels will be normalized. ``reference`` should also
            be provided.  *Default: False*
        reference: :class:`genomvar.Reference`, optional
            Reference genome. *Default: None*
        parse_null: bool, optional
            If ``True`` null variants will also be parsed (experimental).
            *Default: False*
        Returns
        -------
        VariantSet
        """
        reader = VCFReader(vcf,reference=reference)
        records = reader.get_records(parse_info=parse_info,
                                     parse_samples=parse_samples,
                                     normindel=normindel)
        vset = cls.__new__(cls)
        vset.__init__(records['vrt'],records['vcf'],records.get('info'),
                      records.get('sampdata'), reference=reference)
        return vset

    @classmethod
    def from_variants(cls,variants,reference=None):
        """
        Instantiate VariantSet from iterable of variants.
        
        Parameters
        ----------
        variants : iterable
            Variants to add.
        reference : Reference, optional
            Genome reference *Default: None*
        Returns
        -------
        VariantSet
        """
        _variants = cls._rows_from_variants(variants)
        vset = cls.__new__(cls)
        vset.__init__(_variants,reference=reference)
        return vset

    def to_records(self,nested=True):
        """
        Export ``self`` to NumPy ndarray.
        
        Parameters
        ----------

        nested : bool, optional
            Matters only if ``self`` was instatiated with INFO or SAMPDATA.
            When ``True`` dtype of the return will be nested. If ``False``
            no recursive fields will be present and "INFO_" or "SAMP_%SAMP%" will
            be prepended to corresponding fields *Default: True*

        Returns
        -------
        NumPy ndarray with structured dtype
        """
        dtype = [(f,self._variants.dtype.fields[f][0]) for f in \
                      ('chrom','start','end','ref','alt') ]
        dtype.append( ('vartype','U10') )
        dtype.append( ('phase_group',np.int_) )
        
        if not self._info is None:
            if nested:
                dtype.append( ('INFO',self._info.dtype) )
            else:
                for field,(_dtype,offset) in self._info.dtype.fields.items():
                    dtype.append( ('INFO_'+field,_dtype) )
        if not self._sampdata is None:
            if nested:
                samp_dtype = [(s,self._sampdata[s].dtype) for s in self._sampdata]
                dtype.append( ('SAMPLES',samp_dtype) )
            else:
                for samp in self._sampdata:
                    for field,(_dtype,offset) in \
                                self._sampdata[samp].dtype.fields.items():
                        dtype.append( ('SAMPLES_{}_{}'.format(samp,field),
                                       _dtype) )

        ret = np.zeros(shape=self._variants.shape,dtype=dtype)

        for field in ('chrom','start','end','ref','alt','vartype'):
            ret[field] = self._variants[field]
        ret['phase_group'] = self._variants['haplotype']
        ret['vartype'] = [v.__name__ for v in self._variants['vartype']]

        mask = (ret['vartype']==variant.AmbigIns)\
            |(ret['vartype']==variant.AmbigDel)
        ret[mask]['start'] = self._variants[mask]['start2']
        ret[mask]['end'] = self._variants[mask]['end2']

        # now actually populating the data
        if not self._info is None:
            if nested:
                ret['INFO'] = self._info
            else:
                for field in self._info.dtype.fields:
                    ret['INFO_'+field] = self._info[field]

        if not self._sampdata is None:
            if nested:
                for samp in self._sampdata:
                    ret['SAMPLES'][samp] = self._sampdata[samp]
            else:
                for samp in self._sampdata:
                    for field in self._sampdata[samp].dtype.fields:
                        ret['SAMPLES_{}_{}'.format(samp,field)] = \
                                    self._sampdata[samp][field]

        
        nohp_mask = self._variants['vartype']!=variant.Haplotype
        ret = ret[nohp_mask]

        ret['ref'][ret['ref']==''] = '-'
        ret['alt'][ret['alt']==''] = '-'
        return ret


    @property
    def chroms(self):
        """Chromosome set"""
        try:
            return self._chroms
        except AttributeError:
            chroms,ind = np.unique(self._variants['chrom'],return_index=True)
            self._chroms = ChromSet(self._variants['chrom'][np.sort(ind)])
            return self._chroms
    
    def sort_chroms(self,key=None):
        """
        Sorts chromosomes in the set
        
        This can influence order in which methods return variants and 
        variant set comparison.

        Parameters
        ----------
        key : callable, optional
            key to use for chromosome sorting. If not given will sort
            lexicographically. *Default: None*

        Returns
        -------
        None
        """
        if not key:
            ind = np.argsort(self._variants['chrom'])
        else:
            vals = [key(c) for c in self._variants['chrom']]
            ind = np.argsort(vals)
        
        self._variants = np.take(self._variants,ind)
        if self._info:
            self._info = np.take(self._info,ind)
        if self._sampdata:
            for s,val in self._sampdata.items():
                self._sampdata[s] = np.take(self._sampdata[s],ind)
        chroms,ind = np.unique(self._variants['chrom'],return_index=True)
        self._chroms = ChromSet(self._variants['chrom'][np.sort(ind)])
                    
    def _cmp(self,other,action,match_ambig=False,match_partial=True):
        """Returns a new variation set depending on operation
        specified.
        """
        def iterate_aligned(vs1,vs2):
            ind1 = []
            ind2 = []
            for ind,chrom,start,end in zip(*[vs1._variants[k] \
                            for k in ('ind','chrom','start','end')]):
                if chrom not in vs2.chroms:
                    ind1.append(ind)
                    ind2.append(None)
                    continue
                ovlp = vs2._data[chrom].find(start,end)
                if not ovlp:
                    ind1.append(ind)
                    ind2.append(None)
                else:
                    for iv in ovlp:
                        ind1.append(ind)
                        ind2.append(iv.value)
            _ar1 = np.take(vs1._variants,ind1)
            _ar2 = np.take(vs2._variants,
                           list(filter(lambda v: not v is None,ind2)))
            vars1 = zip(*[_ar1[f] for f in vs1._variants.dtype.fields])
            vars2 = zip(*[_ar2[f] for f in vs2._variants.dtype.fields])
            none = (None,)*len(vs2._variants.dtype.fields)
            for ind1,ind2 in zip(ind1,ind2):
                try:
                    if ind2 is None:
                        yield next(vars1),none
                    else:
                        yield next(vars1),next(vars2)
                except StopIteration:
                    return
            # iterate_aligned ends here

        if not isinstance(other,VariantSet):
            msg = '{} is not a VariantSet'.format(other)
            raise TypeError(msg)

        aligned = iterate_aligned(self,other)
        ind2take = []
        vrt2add = []
        for ind1,subit in itertools.groupby(aligned,key=lambda p: p[0][0]):
            var1,var2 = next(subit)
            if var2[0] is None: # array of Nones
                if action=='diff':
                    ind2take.append(var1[0])
                elif action=='comm':
                    continue
            else:
                ind1,h,chrom1,start1,end1,ref1,alt1,vt1,*rest = var1
                if vt1.is_subclass(variant.MNP):
                    locus = [var2,*[t[1] for t in subit]]
                    # For performance treat specially case SNP-SNP match
                    if vt1.is_subclass(variant.SNP) and len(locus)==1 \
                             and locus[0][7].is_subclass(variant.SNP):
                        ind2,h,chrom2,start2,end2,ref2,alt2,vt2,*rest = locus[0]
                        if start1==start2 and alt1==alt2:
                            if action=='comm':
                                ind2take.append(ind1)
                        else:
                            if action=='diff':
                                ind2take.append(ind1)
                    else: # now through variant objects
                        #print('!! here')
                        vrt1 = VariantSet._get_vrt(var1)
                        locus_vrt = [VariantSet._get_vrt(v) for v in locus]
                        mt = ops.matchv(vrt1,locus_vrt,match_partial=match_partial,
                                        match_ambig=match_ambig)
                        comp = ops.cmpv(vrt1,mt,action)
                        vrt2add.extend([(var1[0],v) for v in comp])
                elif vt1.is_subclass(variant.Indel):
                    if not vt1.is_subclass(variant.AmbigIndel):
                        for _,(ind2,h,chrom2,start2,end2,ref2,alt2,vt2,*rest)\
                               in itertools.chain([(var1,var2)],subit):
                            if vt2==vt1 and start1==start2 \
                                       and ref1==ref2 and alt1==alt2:
                                if action=='comm':
                                    ind2take.append(ind1)
                                break
                        else:
                            if action=='diff':
                                ind2take.append(ind1)
                    else:
                        v1 = vt1(chrom=chrom1,start=(start1,rest[0]),
                                 end=(end1,rest[1]),
                                 ref=ref1,alt=alt1)
                        for _,row in itertools.chain([(var1,var2)],subit):
                            vrt = VariantSet._get_vrt(row)
                            if (match_ambig and v1.ambig_equal(vrt)) \
                               or (not match_ambig and v1.edit_equal(vrt)):
                                if action=='comm':
                                    ind2take.append(ind1)
                                break
                        else:
                            if action=='diff':
                                ind2take.append(ind1)
                            
        new = np.take(self._variants,ind2take)
        
        if vrt2add:
            append = VariantSet._rows_from_variants([t[1] for t in vrt2add])
            try:
                new = np.concatenate([new,append])
            except TypeError as exc:
                if new.shape[0]==0:
                    new = append
                else:
                    raise exc
        new['ind'] = np.arange(new.shape[0])
        if not self._info is None:
            info = np.take(self._info,ind2take)
            if vrt2add:
                extra = np.take(self._info,[t[0] for t in vrt2add])
                info = np.concatenate([info,extra])
        else:
            info = None
        if not self._sampdata is None:
            sampdata = {}
            for samp in self._sampdata:
                _sampdata = np.take(self._sampdata[samp],ind2take)
                if vrt2add:
                    extra = np.take(self._sampdata[samp],
                                    [t[0] for t in vrt2add])
                    _sampdata = np.concatenate([_sampdata,extra])
                sampdata[samp] = _sampdata
        else:
            sampdata = None
        rt = self.__new__(self.__class__)
        rt.__init__(new,info=info,sampdata=sampdata)
        return rt
    
    def diff(self,other,match_partial=True,
             match_ambig=False):
        """
        Returns a new variant set object which has variants
        present in ``self`` but absent in the ``other``.

        Parameters
        ----------

        match_partial: bool
            If ``False`` common positions won't be searched in non-identical
            MNPs. *Default: True*
        match_ambig: bool
            If ``False`` ambigous indels will be treated as regular indels.
            *Default: False*

        Returns
        -------
        diff : VariantSet
            different variants

        Notes
        -----
        Here is an example of ``diff`` operation on two variant sets s1 and s2.
        Positions identical to REF are empty::

          REF         GATTGGTAC
          ---------------------
          s1          C  CCC  T
                         CCC  G
          ---------------------
          s2                  T
                        AC    T
          =====================
          s1.diff(s2) C   CC
                          CC  G
          ---------------------
          s2.diff(s1)   A

        """
        return self._cmp(other,'diff',match_partial=match_partial,
                         match_ambig=match_ambig)

    def comm(self,other,match_partial=True,match_ambig=False):
        """
        Returns a new variant set object which has variants
        present both ``self`` and ``other``.

        Parameters
        ----------
        match_partial: bool
            If ``False`` common positions won't be searched in non-identical
            MNPs. *Default: True*
        match_ambig: bool
            If ``False`` ambigous indels will be treated as regular indels.
            *Default: False*

        Returns
        -------
        comm : VariantSet
            Common variants.

        Notes
        -----
        Here is an example of ``diff`` operation on two variant sets s1 and s2.
        Positions identical to REF are empty::

          REF         GATTGGTAC
          ---------------------
          s1          C  CCC  T
                         CCC  G
          ---------------------
          s2                  T
                        AC    T
          =====================
          s1.comm(s2)    C    T
                         C
        """
        return self._cmp(other,'comm',
                         match_partial=match_partial,
                         match_ambig=match_ambig)

class MutableVariantSet(VariantSetBase):
    """
    Mutable class storing variants supporting addition of variants with
    specified genotype.
    
    Supports random access by variant location.
    Can be instantiated from a VCF file (see
    :meth:`genomvar.varset.VariantSet.from_vcf`) or from iterable containing
    variants (see :meth:`VariantSet.from_variants`).

    Attributes
    ----------
    reference: Reference
        Genome reference.
    chroms: ChromSet
        list-like representation of chromosome in the set.
    ctg_len: dict
        dict mapping Chromosome name to corresponding lengths.
        
        This keeps track of the rightmost position of variants in the set
        if reference was not provided. If reference was given ``ctg_len`` 
        is taken from length of chromosomes in reference.
    """
    def __init__(self,variants=None,data=None,reference=None,max_ploidy=2):
        super().__init__(variants=variants,data=data,reference=reference)
        self._variants = variants if variants else {}
        self._obs = {}
        self.max_ploidy = max_ploidy

    def copy(self):
        """Returns a copy of variant set"""
        cpy = self.__new__(self.__class__)
        data = {}
        for chrom,itree in self._data.items():
            data[chrom] = IntervalTree()
            for iv in itree.find(0,self.ctg_len[chrom]):
                data[chrom].add_interval(bxiv(iv.start,iv.end,iv.value))
        cpy.__init__(self._variants.copy(),data,self.reference)
        return cpy

    @classmethod
    def from_vcf(cls,vcf,reference=None,normindel=False,
                 parse_info=False,sample=None,
                 max_ploidy=2,parse_null=False):
        """
        Parse VCF variant file and return VariantSet object
        
        Parameters
        ----------
        vcf : str
            VCF file to read data from.
        parse_info: bool, optional
            If ``True`` INFO fields will be parsed.
            *Default: False*.
        parse_samples: bool, str or list of str, optional
            If ``True`` all sample data will be parsed.  If string, sample 
            with this name will be parsed.  If list of string only these
            samples will be parsed. *Default: False*
        normindel: bool, optional
            If ``True`` indels will be normalized. ``reference`` should also
            be provided.  *Default: False*
        reference: genomvar.Reference, optional
            Reference genome. *Default: None*
        max_ploidy: int, optional
            Indicates how many variants can occupy a single locus 
            (experimental). Default: 2
        parse_null: bool, optional
            If ``True`` null variants will also be parsed (experimental).
            *Default: False*
        Returns
        -------
        VariantSet
        """
        reader = VCFReader(vcf)
        vset = cls.__new__(cls)
        vset.__init__(reference=reference,max_ploidy=max_ploidy)

        if reference is None and normindel:
                raise ValueError('Can\'t norm indels without reference')

        for vrt in reader.iter_vrt(parse_samples=sample if sample else False,
                                   parse_info=parse_info):
            vset.add_vrt(vrt, GT=vrt.attrib['samples'][sample]['GT'] \
                         if sample else (1,0),
                         allow_adjust_genotype=True,attrib=vrt.attrib)

        return vset

    def _add_iv(self,chrom,start,end,id):
        """Insert new interval with ID into ``_data``"""
        iv = bxiv(start,end,value=id)
        self._data.setdefault(chrom,IntervalTree()).add_interval(iv)
        self.ctg_len[chrom] = max(end,self.ctg_len.get(chrom,0))
        return iv

    def _process_gt(self,GT,ovlp,allow_adjust_genotype=False):
        """Function takes a genotype and a list of overlapping
        variants and returns a resulting genotype depending on genotype
        parameter.
        """

        # other_gt is an array of
        # genotypes from overlap
        other_gt = np.zeros(shape=(len(ovlp),
                    max(len(GT),*[len(o.GT) for o in ovlp])))
        for ind,o in enumerate(ovlp):
            other_gt[ind,0:len(o.GT)] = o.GT

        s = np.sum(other_gt,axis=0)
        # If the positions where GT has 1 others have 1 then
        # genotypes will have to be adjust.
        if s[np.nonzero(GT)[0]].sum()>0: # have to change GT
            # zr holds the positions where genotypes from overlap
            # have zeros
            if not allow_adjust_genotype:
                msg = 'GT {:s} not compatible with variants {:s}'\
                    .format(str(GT),str(ovlp))
                raise OverlappingHaplotypeVars(msg)
            zr = np.nonzero(s==0)[0]
            if len(zr)>0: # There is an empty spot to squeeze in new genotype
                GT = np.zeros(other_gt.shape[1],dtype=np.int_)
                GT[zr[0]] = 1
                GT = tuple(GT)
            else: # No space to squeze in GT
                if max([len(v.GT) for v in ovlp])+1>self.max_ploidy:
                    msg = "Can't increase ploidy more than maximum of "\
                        +"{} when inserting around {}"\
                        .format(self.max_ploidy,ovlp)
                    raise OverlappingHaplotypeVars(msg)
                for cur_vrt in ovlp:
                    cur_vrt.GT = (0,)+cur_vrt.GT
                GT = (1,)+(0,)*other_gt.shape[1]
        # Genotypes are compatible as is
        elif s[np.nonzero(GT)[0]].sum()==0:
            pass
        return GT

    def _find_vrt(self,chrom,start=0,end=MAX_END,expand=False):
        """Auxillary function invoked by ``find_vrt``"""
        for iv in self._find_iv(chrom,start,end):
            vrt = self._variants[iv.value]
            yield vrt
            if vrt.is_instance(variant.Haplotype) and expand:
                for _vrt in vrt.find_vrt(start,end):
                    yield _vrt

    def delete_vrt(self,vrt):
        """
        Deletes a variant from variant set.
        """
        if isinstance(vrt,GenomVariant):
            self._obs[vrt.key] = True

    def add_hap_variants(self,variants,GT,attrib=None,
                         allow_adjust_genotype=False):
        """Add variants as a single haplotype defined by genotype GT
        """
        hap = Haplotype.from_variants(variants)
        if attrib:
            for v,attr in zip(hap.variants,attrib):
                if attr:
                    v.attrib = attr
        vrt = self.add_vrt(hap,GT=GT,
                           allow_adjust_genotype=allow_adjust_genotype)
        return vrt
        
    def nof_unit_vrt(self,collapse_homozigosity=True):
        """
        Returns number of simple genomic variations in a variant set. 
        
        SNPs and indels count as 1. MNPs as number of affected nucleotides.

        Parameters
        ----------
        collapse_homozigosity: bool
            If ``True`` same unit variants on different chromosomes 
            are counted as 1. Otherwise counted separately. *Default: True*
        Returns
        -------
        int
        """
        N = 0
        variants = self.find_vrt(expand=True)
        if collapse_homozigosity:
            it1,it2,it3 = itertools.tee(variants,3)
            mnps = filter(lambda v: v.is_instance(variant.MNP),it1)
            insertions = filter(lambda v: v.is_instance(variant.Ins),it2)
            deletions = filter(lambda v: v.is_instance(variant.Del),it3)
            for grp in no_ovlp(mnps):
                if len(grp)==1:
                    N += grp[0].nof_unit_vrt()
                else:
                    N += ops.nof_snp_vrt(grp)

            for _,it in itertools.groupby(insertions,lambda v: (v.start,v.alt)):
                N += 1
            for _,it in itertools.groupby(deletions,lambda v: (v.start,v.ref)):
                N += 1
            for grp in no_ovlp(mnps):
                if len(grp)==1:
                    N += grp[0].nof_unit_vrt()
                else:
                    N += ops.nof_snp_vrt(grp)
        else:
            for vrt in filter(lambda v: not v.is_instance(variant.Haplotype),
                              variants):
                N += vrt.nof_unit_vrt()
        return N

    def add_vrt(self,vrt,GT,allow_adjust_genotype=False,
                attrib=None):
        """
        Insert a variant into variant set.
        
        Parameters
        ----------
        vrt : VariantBase-like
            variant to insert
        GT : tuple
            Genotype of the variant, e.g. (0,1)
        allow_adjust_genotype : bool
            If the variant will overlap by genotype with variants already
            present in the set, setting to ``True`` will allow new 
            variant to fit.
        attrib : dict, optional
            Dictionary with attributes to attached to return
        Returns
        -------
        variant.GenomVariant
        """
        ovlp = self.ovlp(vrt,match_ambig=False)

        if ovlp:
            GT2 = self._process_gt(GT,ovlp,
                                   allow_adjust_genotype=allow_adjust_genotype)
        else:
            GT2 = GT

        # Initializing GenomVariant object with the VariantBase
        ret = GenomVariant(base=vrt if isinstance(vrt,VariantBase) \
                           else vrt.base,GT=GT2,attrib=attrib)

        # Add variant to interval tree and variant dict
        self._add_iv(ret.chrom,ret.start,ret.end,ret.key)
        self._chroms.add(ret.chrom)
        self._variants[ret.key] = ret
        
        return ret

    def sort_chroms(self,key=None):
        self._chroms = ChromSet(sorted(self._chroms,key=key))

    @property
    def chroms(self):
        """Chromosome set"""
        return self._chroms

    def _cmp(self,other,action,match_partial=True,match_ambig=False):
        """Returns a new variation set depending on operation
        specified.
        """
        rt = self.__new__(self.__class__)
        rt.__init__(reference=self.reference)

        for chrom in self.chroms:
            for vrt_self in self.find_vrt(chrom):
                # TODO handle null vrt
                if vrt_self.is_instance(variant.Null):
                    continue
                match = other.match(vrt_self,match_partial=match_partial,
                                    match_ambig=match_ambig)
                vrt2add = ops.cmpv(vrt_self,match,action=action)
                if not vrt2add:
                    continue
                if vrt_self.is_instance(variant.Haplotype):
                    rt.add_hap_variants(vrt2add,GT=vrt_self.GT,
                                        allow_adjust_genotype=True,
                                        allow_adjust_ploidy=True,
                                        attrib = [vrt_self.attrib]*len(vrt2add))
                else:
                    for _vrt in vrt2add:
                        nop = rt.add_vrt(_vrt, GT=vrt_self.GT,
                                         allow_adjust_genotype=True,
                                         attrib=vrt_self.attrib)
        return rt

    def diff(self,other,match_partial=True,
             match_ambig=False):
        """
        Returns a new variant set object which has variants
        present in ``self`` but absent in the ``other``.

        Parameters
        ----------

        match_partial: bool
            If ``False`` common positions won't be searched in non-identical
            MNPs. *Default: True*
        match_ambig: bool
            If ``False`` ambigous indels will be treated as regular indels.
            *Default: False*

        Returns
        -------
        diff : MutableVariantSet
            Different variants.

        Notes
        -----
        Here is an example of ``diff`` operation on two variant sets s1 and s2.
        Positions identical to REF are empty::

          REF         GATTGGTAC
          ---------------------
          s1          C  CCC  T
                         CCC  G
          ---------------------
          s2                  T
                        AC    T
          =====================
          s1.diff(s2) C   CC
                          CC  G
          ---------------------
          s2.diff(s1)   A

        """
        return self._cmp(other,'diff',match_partial=match_partial,
                         match_ambig=match_ambig)

    def comm(self,other,match_partial=True,match_ambig=False):
        """
        Returns a new variant set object which has variants
        present both ``self`` and ``other``.

        Parameters
        ----------
        match_partial: bool
            If ``False`` common positions won't be searched in non-identical
            MNPs. *Default: True*
        match_ambig: bool
            If ``False`` ambigous indels will be treated as regular indels.
            *Default: False*

        Returns
        -------
        comm : MutableVariantSet
            Common variants

        Notes
        -----
        Here is an example of ``diff`` operation on two variant sets s1 and s2.
        Positions identical to REF are empty::

          REF         GATTGGTAC
          ---------------------
          s1          C  CCC  T
                         CCC  G
          ---------------------
          s2                  T
                        AC    T
          =====================
          s1.comm(s2)    C    T
                         C
        """
        return self._cmp(other,'comm',
                         match_partial=match_partial,
                         match_ambig=match_ambig)

class CmpSet(object):
    """
    Class ``CmpSet`` is a lazy implementation of variant set comparison which is
    returned by methods ``diff_vrt`` and ``comm_vrt`` of any variant set class
    instance. Useful for comparisons including large VCF files while keeping
    memory profile low.
    """
    def __init__(self,left,right,op,match_partial=True,match_ambig=False):
        self.left = left
        self.right = right
        self.action = op
        self.match_partial = match_partial
        if (isinstance(self.left,VariantFileSet) \
                 or isinstance(self.right,VariantFileSet)) and match_ambig:
            raise ValueError('Cant match ambigous indels on VariantFileSet'\
                             +' subclasses')
        self.match_ambig = match_ambig

    def _ordered_chroms(self,check_consistency=True):
        """Returns a list indicating how chromosomes should be iterated
        to compare the sets.

        Returns
        -------
        ret : (int, str)
           if ``int`` is 0 left should be advanced if 1 right,
           if 2 then both. ``str`` is the chromosome.
        """
        def _split_at(it,pred):
            buf = []
            for item in it:
                if pred(item):
                    yield buf
                    buf = []
                else:
                    buf.append(item)
            yield buf

        lchrom = self.left.get_chroms()
        rchrom = self.right.get_chroms()
        
        comm_chrom = list(filter(lambda c: c in rchrom, lchrom))
        if check_consistency:
            comm_chrom2 = list(filter(lambda c: c in lchrom, rchrom))
            if comm_chrom != comm_chrom2:
                raise DifferentlySortedChromsError
        

        chroms = []
        
        pred = lambda e: e in comm_chrom
        for left,right,common in zip(_split_at(lchrom,pred),
                                     _split_at(rchrom,pred),
                                     comm_chrom+[None]):
            chroms.extend([(0,e) for e in left]\
                          +[(1,e) for e in right]\
                          +([(2,common)] if not common is None else []))

        return chroms
    
    def region(self,chrom=None,start=0,end=MAX_END,rgn=None):
        """
        Generates variants corresponding to comparison within a region.
        
        Region can be specified using chrom,start,end or rgn with notation like
        ``chr15:1000-2000``

        Parameters
        ----------
        chrom : str
            Chromosome
        start : int, optional
            Start of region. *Default: 0*
        end: int
            End of region. *Default: max chrom length*
        rgn : str
            Region to yield variants from

        Yields
        ------
        :class:`genomvar.variant.GenomVariant`
        """
        if not rgn is None:
            if chrom:
                raise ValueError('rgn or chrom,start,end should be specified')
            chrom,start,end = rgn_from(rgn)
        kwds = {'chrom':chrom,'start':start,'end':end,
                  'check_order':True,}
        for vs in (self.left,self.right):
            if not isinstance(vs,IndexedVariantFileSet) \
                              and not isinstance(vs,MutableVariantSet) \
                              and not isinstance(vs,VariantSet):
                raise ValueError('{} does not support region access'\
                                 .format(repr(vs)))
        return self._iter_chrom_vrt(self.left.find_vrt(**kwds),
                                    self.right.find_vrt(**kwds))

    def iter_vrt(self,callback=None):
        """
        Generates variants corresponding to comparison.

        Parameters
        ----------
        callback: callable
            A function to be called on a match in right variants set 
            for every variant of left variant set. Result will be stored 
            in ``vrt.attrib['cmp']`` of yielded variant

        Yields
        ------
        variants : :class:`genomvar.variant.GenomVariant` or \
            tuple (int, :class:`genomvar.variant.GenomVariant`)
            if ``action`` is ``diff`` or ``comm`` function yields instances of 
            :class:`genomvar.variant.GenomVariant`. If ``action`` is ``all``
            function yields tuples where ``int`` is 0 if variant is present 
            in first set; 1 if in second; 2 if in both.
        """
        lchrom_it = self.left.iter_vrt_by_chrom(check_order=True)
        rchrom_it = self.right.iter_vrt_by_chrom(check_order=True)
        _vs = {0:[self.left], 1:[self.right], 2:[self.left,self.right]}
        if self.action=='diff':
            for which,chrom in self._ordered_chroms():
                if which==0:
                    lchrom,it = next(lchrom_it)
                    for b in it:
                        yield b
                elif which==1:
                    next(rchrom_it)
                elif which==2:
                    lchrom,lit = next(lchrom_it)
                    rchrom,rit = next(rchrom_it)
                    for b in self._iter_chrom_vrt(lit,rit,callback=callback):
                        yield b
        elif self.action=='comm':
            for which,_ in self._ordered_chroms():
                if which==0:
                    next(lchrom_it)
                elif which==1:
                    next(rchrom_it)
                elif which==2:
                    lchrom,lit = next(lchrom_it)
                    rchrom,rit = next(rchrom_it)
                    for b in self._iter_chrom_vrt(lit,rit,callback=callback):
                        yield b
        elif self.action=='all':
            for which,_ in self._ordered_chroms():
                if which==0:
                    lchrom,it = next(lchrom_it)
                    for b in it:
                        yield (0,b)
                elif which==1:
                    rchrom,it = next(rchrom_it)
                    if self.action=='all':
                        for b in it:
                            yield (1,b)
                elif which==2:
                    lchrom,lit = next(lchrom_it)
                    rchrom,rit = next(rchrom_it)
                    for b in self._iter_chrom_vrt(lit,rit,callback=callback):
                        yield b
            
                
    def _iter_chrom_vrt(self,lit,rit,callback=None):
        """Subroutine of :meth:`CmpSet.iter_vrt` acting on sorted variants from
        the same chromosome."""
        def _compare(vrt,ovlp,action):
            return ops.cmpv(vrt,mt,action,callback=callback)
        
        cur_end = {}
        cur_right = []
        cur_pos_right = None
        for ind,vrt,ovlp in zip_variants(lit,rit):
            if self.action=='diff':
                if ind!=0:
                    continue
                if not ovlp:
                    yield vrt
                else:
                    mt = ops.matchv(vrt,ovlp)
                    for vrt in _compare(vrt,mt,self.action):
                        yield vrt
            elif self.action=='comm':
                if ind!=0 or not ovlp:
                    continue
                else:
                    mt = ops.matchv(vrt,ovlp)
                    for vrt in _compare(vrt,ovlp,self.action):
                        yield vrt
            elif self.action=='all':
                if not ovlp:
                    yield(ind,vrt)
                else:
                    mt = ops.matchv(vrt,ovlp)
                    for op in _compare(vrt,mt,'diff'):
                        yield (ind,op)
                    if ind!=0:
                        continue
                    for op in _compare(vrt,mt,'comm'):
                        yield (2,op)
                
class VariantFileSet(VariantSetBase):
    """
    Variant set representing variants contained in underlying VCF file specified
    on initialization.

    Attributes
    ----------
    file : str
        Variant file storing variants.
    """
    # Setting methods not supposed to be implemented to None
    # because official Python docs say so
    comm = None
    diff = None
    _find_vrt = None

    def __init__(self,file,reference=None,parse_info=False,
                 parse_samples=False):
        self._reader=VCFReader(file)
        self.parse_info = parse_info
        self.parse_samples = parse_samples
        super().__init__(reference=reference)

    def get_factory(self):
        """Returns a VariantFactory object. It will inherit the same reference
        as was used for instantiation of ``self``"""
        return super().get_factory(normindel=False)

    def iter_vrt_by_chrom(self,check_order=False):
        """Iterates variants grouped by chromosome."""
        return self._reader.iter_vrt_by_chrom(parse_info=self.parse_info,
                                              parse_samples=self.parse_samples,
                                              check_order=check_order)
    def iter_vrt(self,expand=False):
        """
        Iterate over all variants in underlying VCF.
        
        Parameters
        ----------
        expand : bool, optional
            If ``True`` variants in encountered haplotypes are yielded.
            *Default: False*
        Yields
        -------
        variant : :class:`genomvar.variant.GenomVariant`
        """
        return self._reader.iter_vrt(parse_info=self.parse_info,
                                     parse_samples=self.parse_samples)

    def nof_unit_vrt(self):
        # """
        # Returns number of simple genomic variations in a variant set. 
        
        # Returns
        # -------
        # N : int

        # Notes
        # -----
        # See documentation on ref:`variant_scoring` for details.
        # """
        N = 0
        variants = self.find_vrt(expand=True)
        for vrt in filter(lambda v: not v.is_instance(variant.Haplotype),
                          variants):
            N += vrt.nof_unit_vrt()
        return N

    def get_chroms(self):
        return self._reader.get_chroms(unindexed=True)

class IndexedVariantFileSet(VariantFileSet,MutableVariantSet):
    """
    Variant set representing variants contained in underlying VCF file specified
    on initialization. By default it searches for Tabix index and uses it 
    for random access of variants.
    """
    def __init__(self,file,index=None,parse_info=False,parse_samples=False,
                 reference=None):
        super().__init__(file,reference=reference,parse_info=parse_info,
                         parse_samples=parse_samples)
        self._reader=VCFReader(file,index=True)
        self._chroms = set(self._reader.chroms)
        
    def _find_vrt(self,chrom,start,end,expand=False):
        """Auxillary function for ``find_vrt``"""
        for vrt in self._reader.find_vrt(
                chrom,start,end,parse_info=self.parse_info,
                parse_samples=self.parse_samples):
            yield vrt
