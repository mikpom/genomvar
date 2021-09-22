"""
Module ``varset`` defines variant set classes which are containers of variants
supporting set-like operations, searching variants, export to NumPy...

There are two classes depending on whether all the variants are loaded
into memory or not:

    - :class:`~genomvar.varset.VariantSet` in-memory

    - :class:`~genomvar.varset.VariantSetFromFile` can read VCF files
      on demand, needs index for random access
"""
import warnings
from copy import copy
import itertools
from functools import lru_cache
import collections
import numpy as np
from rbi_tree.tree import ITree
from genomvar import variant, variantops as ops, MAX_END
from genomvar.variant import VariantBase,AmbigIndel,\
    GenomVariant,Haplotype,VariantFactory
from genomvar.utils import rgn_from,zip_variants,\
    chunkit,no_ovlp
from genomvar.vcf import (_get_reader, VCFReader, VCFWriter,
                          dtype0,dtype1, )
from genomvar.vcf_utils import (header as vcf_header,
                                row_tmpl, 
                                field_writer_simple,
                                dtype2string, string2dtype,
                                VCF_INFO_OR_FORMAT_DEF_FIELDS)
from genomvar import OverlappingHaplotypeVars,\
    NoIndexFoundError, \
    Reference, DifferentlySortedChromsError,\
    DuplicateVariants, SINGLETON, ChromSet

class _VariantSetBase(object):
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


    def __init__(self,variants=None,data=None,reference=None,
                 samples=None):
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

        self._samples = samples
        self._data = data if data else {}
        self._factory = VariantFactory(reference=self.reference,
                                       normindel=False)
        self._fac = {'nonorm':self._factory}
        
    def get_chroms(self):
        return self.chroms
    
    def _find_iv(self,chrom,start=0,end=MAX_END):
        if chrom not in self._data:
            return iter([])
        for ivl in self._data[chrom].find(start,end):
            yield ivl

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
            max chrom length*
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

    def _find_chrom_vrt(self,chrom,expand=False):
        for vrt in self._find_vrt(chrom,expand=expand):
            yield vrt

    def iter_vrt_by_chrom(self,**kwds):
        """Iterates variants grouped by chromosome."""
        for chrom in self.chroms:
            yield (chrom,self._find_chrom_vrt(chrom))

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
            for vrt in self._find_chrom_vrt(chrom,expand=expand):
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
        if not match_ambig and vrt0.is_variant_instance(variant.AmbigIndel):
            start = vrt0.act_start
            end = vrt0.act_end
        else:
            start = vrt0.start
            end = vrt0.end

        # Insertions need special treatment due to arbitrariness of
        # position definition
        if vrt0.is_variant_instance(variant.Ins):
            for _vrt in self._find_vrt(chrom,start,end):
                if _vrt.is_variant_instance(variant.AmbigIndel) and not match_ambig:
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
                elif end<=_vrt_start+1 and not _vrt.is_variant_instance(variant.Ins):
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
            if _vrt.is_variant_instance(variant.AmbigIndel) and not match_ambig:
                _vrt_end = _vrt.act_end
                _vrt_start = _vrt.act_start
            else:
                _vrt_end = _vrt.end
                _vrt_start = _vrt.start
            # if _vrt is insertion and the search region is further
            if _vrt.is_variant_instance(variant.Ins) and _vrt_end<=start+1:
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
        
        The method checks whether variant ``vrt``
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
    
    def to_vcf(self, out, reference=None, info_spec=None, format_spec=None,
               samples=None):
        """Writes a minimal VCF with variants from the set. 

        INFO and SAMPLE data from source data is not preserved.

        Parameters
        ----------

        out : handle or str
           If string then it's path to file, otherwise a handle 
           to write variants.

        reference : Reference or str, optional
           If string then it's path to reference FASTA. Otherwise
           object of :class:`genomvar.Reference` is expected. Not
           necessary if self was instantiated with a reference.
           If given, this argument takes precedence over the ref 
           used at instantiation.

        info_spec : list of tuples containing specification of 
            an INFO field in Order per VCF spec, i.e. `NAME`,  
            `NUMBER`,  `TYPE`,  `DESCRIPTION`,  `SOURCE`,  `VERSION`.
            Only the first three fiels (name, number and type) 
            are required. If variant set was instantiated from VCF
            file spec obtained from it will be used. Providing 
            empty sequence will omit INFO fields in output.
        
        format_spec : spec of FORMAT fields, similar to info_spec

        Returns
        -------
        None

        Examples
        --------

        >>> with open(fname, 'wt') as fh:
        ...     vs.to_vcf(fh, info_spec=[('DP4', 4, 'Integer'),
        ...                              ('NSV', 1, 'Integer')])
        """
        if isinstance(out, str):
            fh = open(out,'wt')
            opened = True
        else:
            fh = out
            opened = False

        specs = {'info_spec' : info_spec,
                 'format_spec' : format_spec}
        FIELDS = VCF_INFO_OR_FORMAT_DEF_FIELDS
        
        for spec in specs:
            if specs[spec] is None:
                # in this case check maybe VariantSet was
                # instantiated with dtype
                key = 'info' if spec=='info_spec' else 'format'
                if not self.dtype is None and key in self.dtype:
                # if key in self.dtype:
                    specs[spec] = [tuple([v.get(f.lower()) for f in FIELDS]) \
                        for v in self.dtype[key].values()]

        if not samples is None:
            if self._samples:
                for sample in samples:
                    if sample not in self._samples:
                        raise ValueError(
                            'Sample {} not found in {}'.format(
                                sample, self._samples))
        writer = VCFWriter(
            reference if not reference is None else self.reference,
            samples=samples if not samples is None else self._samples,
            **specs)

        fh.write(writer.get_header())
        self._write_variants(fh, writer)

        if opened:
            fh.close()

    def _write_variants(self, fh, writer):
        for vrt in self.iter_vrt(expand=True):
            try:
                row = writer.get_row(vrt)
            except ValueError as exc: # TODO more specific error
                if vrt.is_variant_instance(variant.Haplotype) \
                        or vrt.is_variant_instance(variant.Asterisk):
                    continue
                else:
                    raise exc
            fh.write(str(row)+'\n')

    def get_factory(self,normindel=False):
        """Returns a VariantFactory object. It will inherit the same reference
        as was used for instantiation of ``self``"""
        try:
            return self._fac['norm' if normindel else 'nonorm']
        except KeyError:
            vf = VariantFactory(self.reference,normindel=normindel)
            self._fac['norm' if normindel else 'nonorm'] = vf
            return vf

class VariantSet(_VariantSetBase):
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
    def __init__(self,variants, attrib=None, vcf_notation=None, info=None,
                 sampdata=None, reference=None, dtype=None, samples=None):
        super().__init__(reference=reference, samples=samples)
        self._variants = variants
        # if not attrib is None:
        self._attrib = attrib
        # if not dtype is None:
        self.dtype = dtype
        self._info = info
        self._sampdata = sampdata
        self._vcf_notation = vcf_notation
        fields = ('ind','haplotype','chrom','start','end')
        for ind,hap,chrom,start,end in zip(
                    *[self._variants[f] for f in fields]):
            if hap!=SINGLETON:
                continue
            t = self._data.setdefault(chrom, ITree())
            t.insert(start, end, value=ind)
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

    def _genom_vrt_from_row(self,row):
        base = self._vrt_from_row(row)
        attrib = self._get_attrib(row[0])
        return variant.GenomVariant(base=base,attrib=attrib)

    # def _vrt_from_row(self,row):

    def _get_attrib(self, rnum):
        if not self._attrib is None:
            return self._attrib[rnum]

        EXPOSED_ATTRIB = ['id', 'qual', 'filter']
        if not self._info is None:
            attrib = {'info':self._info[rnum]}
        else:
            attrib = {}
        if not self._vcf_notation is None:
            d = dict(zip(dtype1.names,
                         self._vcf_notation[rnum].tolist()))
            d2 = {k:d.pop(k) for k in EXPOSED_ATTRIB}
            attrib['vcf_notation'] = d
            attrib.update(d2)
        if not self._sampdata is None:
            attrib.update({'samples':{s:self._sampdata[s][rnum]\
                                      for s in self._sampdata}})
        return attrib

        
    def _vrt_from_row(self,row):
        ind,hap,chrom,start,end,ref,alt,cls,start2,end2 = row
        try:
            return cls(chrom=chrom,start=start,end=end,ref=ref,alt=alt)
        except TypeError as exc:
            if cls.is_variant_subclass(variant.AmbigIndel):
                return cls(chrom=chrom,start=(start,start2),
                           end=(end,end2),ref=ref,alt=alt)
            elif cls.is_variant_subclass(variant.Haplotype):
                vrt = list(self._get_hap_variants(ind))
                return variant.Haplotype(chrom,vrt)
            else:
                raise exc
                
    # @lru_cache(maxsize=1000)
    # def _get_vrt(self,rnum):
    #     """Constructs a variant from row number.
    #     If haplotypes its subvariants will be GenomVariant instances."""
    #     return self._vrt_from_row(self._variants[rnum])

    @lru_cache(maxsize=1000)
    def _get_vrt(self,rnum):
        """Constructs a variant from row number.
        If haplotypes its subvariants will be GenomVariant instances."""
        ind,hap,chrom,start,end,ref,alt,cls,start2,end2 = self._variants[rnum]
        try:
            return cls(chrom=chrom,start=start,end=end,ref=ref,alt=alt)
        except TypeError as exc:
            if cls.is_variant_subclass(variant.AmbigIndel):
                return cls(chrom=chrom,start=(start,start2),
                           end=(end,end2),ref=ref,alt=alt)
            elif cls.is_variant_subclass(variant.Haplotype):
                vrt = list(self._get_hap_variants(ind))
                return variant.Haplotype(chrom,vrt)
            else:
                raise exc

    @lru_cache(maxsize=1000)
    def _get_genom_vrt(self,rnum):
        """Like _get_vrt but also adds attributes."""
        base = self._get_vrt(rnum)
        attrib = self._get_attrib(rnum)
        return GenomVariant(base=base,attrib=attrib)

    @staticmethod
    def _get_tups(vrt_iter):
        """Constructs tuples from variants."""
        cnt = 0
        for vrt in vrt_iter:
            if vrt.is_variant_instance(variant.Haplotype):
                ind = cnt
                yield (ind, SINGLETON, *vrt.tolist(),
                       getattr(vrt, 'attrib', None))
                cnt += 1
                for child in vrt.variants:
                    a = child.tolist()
                    yield (cnt, ind, *a, getattr(vrt, 'attrib', None))
                    cnt += 1
            else:
                a = vrt.tolist()
                b = (cnt, SINGLETON, *a, getattr(vrt, 'attrib', None))
                yield b
                cnt += 1

    @staticmethod
    def _rows_from_variants(variants):
        """Returns rows given variants"""
        zipped = list(zip(
            *VariantSet._get_tups(variants))) # 0 for ambig
        
        _variants = np.zeros(dtype=dtype0,shape=len(zipped[0]))
        for ind,field in enumerate(dtype0.fields):
            try:
                _variants[field] = zipped[ind]
            except IndexError:
                continue
        attrib = zipped[-1]
        return _variants, attrib

    def _get_hap_variants(self,hapid):
        ind = hapid+1
        while ind<self._variants.shape[0]:
            row = self._variants[ind].tolist()
            if row[1]==hapid:
                yield self._get_genom_vrt(row[0])
            else:
                break
            ind += 1
        return

    def sample(self, size):
        """
        Returns a random sample of variants.
        
        Haplotypes (if present) maybe returned equally likely as non-haplotype
        variants and their subvariants.

        Parameters
        ----------
        size : int
            size of the sample

        Returns
        -------
        marray : list
            List of variants.
        """
        return [self._get_genom_vrt(r) for r in \
                np.random.choice(self._variants['ind'], size)]
        

    def iter_vrt(self,expand=False):
        rows = self._rows_iter_vrt(expand=expand)
        for rown in rows:
            yield self._get_genom_vrt(rown)

    def _rows_iter_vrt(self, expand=False):
        """Returns indexes of rows iterating all variants"""
        ivl = itertools.chain.from_iterable(
                [self._data[c].iter_ivl() for c in self.chroms])
        ind = [i[2] for i in ivl]
        vrt = self._variants[ind]
        idx = vrt[vrt['haplotype']==SINGLETON]['ind']
        return idx
            

    def _find_vrt(self,chrom,start=0,end=MAX_END,expand=False):
        """Uses ITree to find variants and yields them.
        Implementation is variant set class-specific."""
        try:
            intervals = self._data[chrom].find(start,end)
        except KeyError:
            return
        for s,e,rown in intervals:
            vrt = self._get_genom_vrt(rown)
            yield vrt
            if vrt.is_variant_instance(variant.Haplotype) and expand:
                for _vrt in vrt.find_vrt(start,end):
                    yield _vrt

    def _find_chrom_vrt(self,chrom,expand=False):
        try:
            ivl = self._data[chrom].iter_ivl()
        except KeyError:
            return iter([])
        inds = [i[2] for i in ivl]
        vrt = self._variants[inds]
        for row in zip(*[vrt[f] for f in self._variants.dtype.fields]):
            if row[1]!=SINGLETON:
                continue
            if row[-3]==variant.Haplotype:
                hap = self._get_genom_vrt(row[0])
                yield hap
                if expand:
                    for vrt in hap.variants:
                        yield vrt
            else:
                vrt = self._genom_vrt_from_row(row)
                yield vrt
                
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
        Mixed variants (if present) are count also as 1. 
        """
        N = 0
        for start,end,cls in zip(self._variants['start'],self._variants['end'],
                                 self._variants['vartype']):
            if cls.is_variant_subclass(variant.Haplotype):
                continue
            else:
                if cls.is_variant_subclass(variant.Indel):
                    N += 1
                elif cls.is_variant_subclass(variant.MNP):
                    N += (end-start)
                elif cls.is_variant_subclass(variant.Mixed):
                    N += 1
        return N

    def drop_duplicates(self,return_dropped=False):
        """
        Remove non-unique variants.
        
        Returns a variant set where no pair of variants
        is edit-equal. Haplotypes are left as is.

        Parameters
        ----------
        return_dropped : bool
            Whether to return a list of dropped variants.  *Default: False*

        Returns
        -------
        vs : VariantSet
            Variant set with unique variants.
        dropped : list of variants, optional
            if ``return_dropped`` is True dropped variants are also returned.
        """
        ret = {}
        mask = (self._variants['haplotype']==SINGLETON)\
            &(self._variants['vartype']!=variant.Haplotype)
        dedup = self._variants[mask]
        
        vrt = {}
        fields = {}
        # Insertions
        fields['ins'] = ['chrom','start','alt']
        reg_ins = dedup[dedup['vartype']==variant.Ins].copy()
        ambig_ins = dedup[dedup['vartype']==variant.AmbigIns].copy()
        ambig_ins['start'] = ambig_ins['start2']
        vrt['ins'] = np.concatenate([reg_ins,ambig_ins])

        # Deletions
        fields['dels'] = ['chrom','start','end']
        reg_dels = dedup[dedup['vartype']==variant.Del].copy()
        ambig_dels = dedup[dedup['vartype']==variant.AmbigDel].copy()
        ambig_dels['start'] = ambig_dels['start2']
        ambig_dels['end'] = ambig_dels['end2']
        vrt['dels'] = np.concatenate([reg_dels,ambig_dels])

        # Rest
        fields['rest'] = ['vartype','chrom','start','end','alt']
        rest = dedup[(~np.isin(dedup['vartype'],[variant.AmbigIns,variant.Ins,
                                           variant.Del,variant.AmbigDel]))]\
                                           .copy()
        rest['vartype'] = [c.__name__ for c in rest['vartype']]
        vrt['rest'] = rest

        uniq = {}
        for tp in ('ins','dels','rest'):
            u,idx = np.unique(vrt[tp][fields[tp]],return_index=True)
            uniq[tp] = np.take(vrt[tp],idx)['ind']
            
        ind = np.concatenate((self._variants[~mask]['ind'],
                              *uniq.values()))

        vrt = self._variants[ind]
        vrt['ind'] = np.arange(vrt.shape[0])

        vcf_not = self._vcf_notation[ind] if not self._vcf_notation \
                     is None else None
        info = self._info[ind] if not self._info is None else None
        sampdata = {s:self._sampdata[s][ind] for s in self._sampdata} \
            if not self._sampdata is None else None
        ret = VariantSet(vrt,vcf_not,info,sampdata)
        if return_dropped:
            dropped_ind = self._variants['ind']\
                [~np.isin(self._variants['ind'],ind)]
            dropped = [self._get_genom_vrt(i) for i in dropped_ind]
            return ret,dropped
        else:
            return ret

    @classmethod
    def from_vcf(cls,vcf,reference=None,parse_info=False, parse_samples=False,
                 normindel=False,duplicates='ignore',parse_null=False):
        """
        Parse VCF variant file and return VariantSet object
        
        Parameters
        ----------
        vcf : str
            VCF file to read data from.
        reference: :class:`genomvar.Reference`, optional
            Reference genome. *Default: None*
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
        duplicates: {'ignore','keepfirst','raise'}
            How to treat duplicate variants. If ``keepfirst`` only the first
            is taken. If ``raise`` error is raised. *Default: ``ignore`` and 
            do not check for duplicates*.
        parse_null: bool, optional
            If ``True`` null variants will also be parsed (experimental).
            *Default: False*
        Returns
        -------
        VariantSet
        """
        reader = _get_reader(vcf, reference=reference)
        records = reader.get_records(parse_info=parse_info,
                                     parse_samples=parse_samples,
                                     normindel=normindel)
        if duplicates!='ignore':
            records = cls._check_duplicates(records, duplicates=duplicates)
        # else:
        #     vrt = records['vrt']
        #     vcf_notation = records['vcf']
        #     info = records.get('info')
        #     sampdata = records.get('sampdata')
            
        vset = cls.__new__(cls)
        vset.__init__(records['vrt'],
                      vcf_notation=records['vcf'],
                      info=records.get('info'),
                      sampdata=records.get('sampdata'),
                      reference=reference,
                      dtype=reader.dataparser.dtype,
                      samples=reader.samples)
        return vset

    @classmethod
    def from_variants(cls, variants, reference=None):
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
        _variants, attrib = cls._rows_from_variants(variants)
        # print('!!', _variants, attrib)
        vset = cls.__new__(cls)
        vset.__init__(_variants, attrib=attrib, reference=reference)
        return vset

    def to_records(self, nested=True):
        """
        Export ``self`` to NumPy ndarray.
        
        Parameters
        ----------

        nested : bool, optional
            Matters only if ``self`` was instatiated with INFO or SAMPDATA.
            When ``True`` dtype of the return will be nested. If ``False``
            no recursive fields will be present and "info_" or "SAMP_%SAMP%" will
            be prepended to corresponding fields *Default: True*

        Returns
        -------
        NumPy ndarray with structured dtype
        """
        dtype = [(f,self._variants.dtype.fields[f][0]) for f in \
                      ('chrom','start','end','ref','alt') ]
        dtype.append( ('vartype','U10') )
        dtype.append( ('phase_group',np.int_) )

        if not self._attrib is None:
            dtype.append( ('attrib', np.object_) )
        else:
            if not self._info is None:
                if nested:
                    dtype.append( ('info',self._info.dtype) )
                else:
                    for field,(_dtype,offset) in self._info.dtype.fields.items():
                        dtype.append( ('info_'+field,_dtype) )
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
        if not self._attrib is None:
            ret['attrib'] = self._attrib
        else:
            if not self._info is None:
                if nested:
                    ret['info'] = self._info
                else:
                    for field in self._info.dtype.fields:
                        ret['info_'+field] = self._info[field]

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

    @staticmethod
    def _iterate_aligned(vs1,vs2):
        def _iterate(ar1,ar2):
            vars1 = zip(*[ar1[f] for f in vs1._variants.dtype.fields])
            vars2 = zip(*[ar2[f] for f in vs2._variants.dtype.fields])
            while True:
                try:
                    yield next(vars1),next(vars2)
                except StopIteration:
                    return
        ind1 = []
        ind2 = []
        fields = ('ind','haplotype','chrom','start','end')
        nomatch = []
        vrt = vs1._variants[np.argsort(vs1._variants['chrom'])]
        for ind,hap,chrom,start,end in zip(
                *[vrt[k] for k in fields]):
            if hap!=SINGLETON:
                continue
            if chrom not in vs2.chroms:
                nomatch.append(ind)
                continue
            ovlp = vs2._data[chrom].find(start,end)
            if not ovlp:
                nomatch.append(ind)
            else:
                ind1.extend([ind]*len(ovlp))
                ind2.extend([v[2] for v in ovlp])

        if ind2:
            ar1 = np.take(vs1._variants,ind1)
            ar2 = np.take(vs2._variants,ind2)
            it = _iterate(ar1, ar2)
        else:
            it = iter([])
        return nomatch,it

    def _cmp(self,other,action,match_ambig=False,match_partial=True):
        """Returns a new variant set depending on operation
        specified.
        """
        def _quickly_match(v1,v2):
            """Return True if match, False if not and None if unclear
            without object instantiation"""
            ind1,h1,chrom1,start1,end1,ref1,alt1,vt1,*rest = v1
            ind2,h2,chrom2,start2,end2,ref2,alt2,vt2,*rest = v2
            if vt1.is_variant_subclass(variant.Haplotype) \
                          or vt2.is_variant_subclass(variant.Haplotype):
                return

            if vt1.is_variant_subclass(variant.SNP):
                if vt2==vt1:
                    if start1==start2 and alt1==alt2:
                        return True
                    else:
                        return False
                elif vt2==variant.MNP:
                    return
                else: # Indel or Mixed
                    return False
            elif vt1==variant.Del:
                if vt2==vt1:
                    if start1==start2 and end1==end2:
                        return True
                    else:
                        return False
                elif vt2.is_variant_subclass(variant.AmbigDel):
                    return
                else:
                    return False
            elif vt1==variant.Ins:
                if vt2==vt1:
                    if start1==start2 and alt1==alt2:
                        return True
                    else:
                        return False
                elif vt2.is_variant_subclass(variant.AmbigIns):
                    return
                else:
                    return False
            else:
                return
            # _quicky_match ends here

        if not isinstance(other,VariantSet):
            msg = '{} should be instance of VariantSet'.format(other)
            raise TypeError(msg)

        nomatch, aligned = self._iterate_aligned(self,other)
        if action=='diff':
            ind2take = nomatch
        else:
            ind2take = []
        vrt2add = []
        for ind1,subit in itertools.groupby(aligned,key=lambda p: p[0][0]):
            vars1, locus = list(zip(*subit))
            var1 = vars1[0] # vars1 repeats the same variant len(locus) times
            if len(locus)==1:
                maybe_equal = _quickly_match(var1,locus[0])
            else:
                maybe_equal = None

            if not maybe_equal is None:
                if maybe_equal:
                    if action=='comm':
                        ind2take.append(ind1)
                elif action=='diff':
                    ind2take.append(ind1)
            else: # last resort, instantiating the objects
                vrt1 = self._get_vrt(var1[0])
                locus_vrt = [other._get_vrt(r[0]) for r in locus]
                mt = ops.matchv(
                    vrt1, locus_vrt, match_partial=match_partial,
                    match_ambig=match_ambig)
                comp = ops.cmpv(vrt1,mt,action)
                # TODO preserve haplotypes
                vrt2add.extend([(var1[0],v) for v in comp])

        new = np.take(self._variants,ind2take)
        if vrt2add:
            append, append_attrib = VariantSet._rows_from_variants([t[1] for t in vrt2add])
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


    def _write_variants(self, fh, writer):
        idx = self._rows_iter_vrt(expand=True)
        infostringsgen = self._format_infostrings(writer, idx)
        sampdatagen = self._format_sampdata(writer, idx)
        variants = [self._get_genom_vrt(i) for i in idx] # generator
        for vrt, infostring, sampdata in \
                           zip(variants, infostringsgen, sampdatagen):
 
            try:
                row = writer.get_row(vrt, info=infostring, sampdata=sampdata)
            except ValueError as exc: # TODO more specific error
                if vrt.is_variant_instance(variant.Haplotype) \
                        or vrt.is_variant_instance(variant.Asterisk):
                    continue
                else:
                    raise exc
            fh.write(str(row)+'\n')

    def _format_infostrings(self, writer, idx):
        def _use_writer(writer, field, values):
            for v in values:
                yield writer(field, v)

        if not self._info is None: # case of NumPy array stored data
            fwriters = writer.writers['info']
            field_names = list(writer.writers['info'])
            # writers = [writer.writers[f] for f in field_names]
            fields = [_use_writer(fwriters[f], f,  self._info[idx][f]) \
                      for f in field_names]
            
            records = zip(*fields)
            
            formatted = (';'.join(r) for r in records)
        elif not self._attrib is None: # case for dictionary stored data
            field_names = list(writer.writers['info']) if hasattr(writer, 'writers') \
                else []
            if len(field_names)==0:
                formatted = ['.' for i in idx]
            else:
                
                records = (a.get('info', {f:'' for f in field_names}) \
                       for a in [self._attrib[i] for i in idx])
                formatted = (';'.join(
                    [writer.writers['info'][f](f, rec[f]) for f in field_names]) \
                             for rec in records)
        else:
            formatted = ('.' for i in idx)

        return formatted

    def _format_sampdata(self, writer, idx):
        def _use_writer(writer, values):
            for v in values:
                yield writer(v)

        # FIXME make storage scheme more explicit
        if not self._sampdata is None: # case of NumPy array stored data
            ffields = list(writer.writers) # FIXME
            # writers = [writer.writers[f] for f in ffields]
            fields = [_use_writer(writer.writers[f], self._info[idx][f]) \
                      for f in ffields]
            records = zip(*fields)
            formatted = (';'.join(r) for r in records)
        elif hasattr(self, '_attrib') and not self._attrib is None: # case for dictionary stored data
            ffields = list(writer.writers['format']) if hasattr(writer, 'writers') \
                else []
            
            if len(ffields)==0:
                formatted = itertools.repeat(None)
            else:
                fn = ':'.join(ffields)
                dflt = {f:None for f in ffields}
                fwriters = writer.writers['format']
                by_sample_strings = []
                for sample in writer.samples:
                    sample_strings = []
                    for _i in idx:
                        dct = self._attrib[_i]['samples'].get(sample, dflt)
                        string = ':'.join(
                            [fwriters[f](dct.get(f)) for f in ffields])
                        sample_strings.append(string)
                    by_sample_strings.append(sample_strings)
                formatted = [fn+'\t'+'\t'.join(a) \
                             for a in zip(*by_sample_strings)]
        else:
            formatted = itertools.repeat(None)

        return formatted
    
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
        self.match_ambig = match_ambig

    def _ordered_chroms(self,check_consistency=True):
        """Returns a list indicating how chromosomes should be iterated
        to compare the sets.

        Returns
        -------
        (which, chrom) : (int, str)
           if ``which`` is 0 left should be advanced, if 1 then right,
           if 2 then both. ``chrom`` is the chromosome.
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
            if not isinstance(vs,VariantSetFromFile) \
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
        try:
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
                        for b in self._iter_chrom_vrt(
                                lit,rit,callback=callback):
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
                        for b in self._iter_chrom_vrt(
                                lit,rit,callback=callback):
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
                        for b in self._iter_chrom_vrt(
                                lit,rit,callback=callback):
                            yield b
        finally:
            lchrom_it.close()
            rchrom_it.close()
                
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
                
class VariantSetFromFile(_VariantSetBase):
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

    def __init__(self, file, reference=None, index=None, parse_info=False,
                 parse_samples=False):
        """
        VariantSetFromFile is instantiated from files containing
        variants (VCF). 

        Optionally, index and reference arguments can be given.

        Parameters
        ----------
        index : str or bool
            By default tries ``vcf``.tbi as index file. 
            if index=True and index not found NoIndexFoundError
            is raised.
            A string with path to index can be given. *Default: None*

        reference : Reference or str
            Genome reference. *Default: None*

        parse_info : bool
            Indicates if info fields should be parsed. *Default: False*

        parse_samples : bool
            Indicates if sample fields should be parsed. *Default: False*

        Returns
        -------
        reader : VCFReader
            VCFReader object.
        """
        self.vcfreader=_get_reader(file, index)
        self.parse_info = parse_info
        self.parse_samples = parse_samples
        super().__init__(reference=reference)
        if index:
            self._chroms = set(self.vcfreader.chroms)

    def get_factory(self):
        """Returns a VariantFactory object. It will inherit the same reference
        as was used for instantiation of ``self``"""
        return super().get_factory(normindel=False)

    def iter_vrt_by_chrom(self,check_order=False):
        """Iterates variants grouped by chromosome."""
        return self.vcfreader.iter_vrt_by_chrom(
                parse_info=self.parse_info,
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
        return self.vcfreader.iter_vrt(
            parse_info=self.parse_info,
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
        for vrt in filter(lambda v: not v.is_variant_instance(variant.Haplotype),
                          variants):
            N += vrt.nof_unit_vrt()
        return N

    def get_chroms(self):
        return self.vcfreader.get_chroms(allow_no_index=True)

    def _find_vrt(self, chrom, start, end, expand=False):
        """Auxillary function for ``find_vrt``"""
        try:
            for vrt in self.vcfreader.find_vrt(
                    chrom, start, end, parse_info=self.parse_info,
                    parse_samples=self.parse_samples):
                yield vrt
        except NoIndexFoundError:
            raise ValueError('Index is required to random access variants')

    @property
    def chroms(self):
        """Chromosome set"""
        return self._chroms

# class VariantSetFromFile(VariantSetFromFile):
#     """
#     Variant set representing variants contained in underlying VCF file specified
#     on initialization. By default it searches for Tabix index and uses it 
#     for random access of variants.
#     """
#     def __init__(self,file,index=None,parse_info=False,parse_samples=False,
#                  reference=None):
#         super().__init__(file,reference=reference,parse_info=parse_info,
#                          parse_samples=parse_samples)
#         self._reader=VCFReader(file,index=True)
#         self._chroms = set(self._reader.chroms)
        
