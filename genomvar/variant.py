"""
Module :mod:`genomvar.variant` contains classes representing genomic
alterations.

The hierarchy of the classes used in the package is the following::

                      VariantBase
                      /      |   \\
     AmbigIndel <-Indel      V   |
       |     |      / \     MNP  V
       | ----+--- Del  Ins   |  Haplotype
       | |   |          |    V
       V V   |          V   SNP
    AmbigDel -----> AmbigIns

All variants have ``start`` and ``end`` attributes defining a range they
act on and can be searched overlap for.

To test whether a variant is instance of some type
:meth:`~VariantBase.is_instance` method can be used.  Variant equality
can be tested using :meth:`~VariantBase.edit_equal`.

Objects can be instantiated directly, e.g.::

    >>> vrt = variant.MNP('chr1',154678,'GT')
    >>> print(vrt)
    <MNP chr1:154678-154680 NN/GT>

This will create an MNP which substitutes positions 154678 and 154679 on
chromosome 1 for ``GT``.

Alternatively variants can be created using :class:`VariantFactory`
objects.  This class can work with VCF-like notation. For example, ::

    >>> fac = VariantFactory()
    >>> vrt = fac.from_edit('chr15',575,'TA','T')
    >>> print(vrt)
    <Del chr15:576-577 A/->

Position is **0-based** so this creates a deletion of A in 577th
position of chromosome 15.

Alternatively, limited subset of HGVS notation is supported (numbering
in HGVS strings is **1-based** following the spec)::

    >>> vrt = fac.from_hgvs('chr1:g.15C>A')
    >>> print(vrt)
    <SNP chr1:14 C/A>

Variant sets defined in :mod:`genomvar.varset` use class :class:`GenomVariant`.
Objects of this class contain genomic alteration (attribute ``base``)
and optionally, genotype (attribute ``GT``) and other attributes
commonly found in VCF files (attribute ``attrib``).  Attribute
``base`` is an object of some ``VariantBase`` subclass (SNPs, Deletions
etc.).
"""
from collections import OrderedDict
import itertools
import re
from bx.intervals.intersection import Interval as bxiv
from bx.intervals.intersection import IntervalTree
from genomvar.utils import _strip_ref_alt
from genomvar import Reference,MAX_END

hgvs_regex = {'SNP':'([0-9]+)([AGTC])>([AGTC])',
              'Del':'([0-9]+)_([0-9]+)del',
              'Ins':'[0-9]+_([0-9]+)ins([AGTC]+)',
              'DelIns':{0:'([0-9]+)_([0-9]+)delins([ATGC]+)',
                        1:'([0-9]+)delins([ATGC]+)'}}


class VariantBase(object):
    """
    Base class for the other genomic variant classes.
    """
    def __init__(self,chrom,start,end,ref,alt):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.ref = ref
        self.alt = alt

    __slots__ = 'chrom','start','end','ref','alt','_key'

    def __str__(self):
        return '<{:s} {:s}:{:d}-{:d} {:s}/{:s}>'\
            .format(type(self).__name__,self.chrom,self.start,self.end,
                    self.ref if self.ref else '-',
                    self.alt if self.alt else '-')
    def __repr__(self):
        return '{}(chrom="{}",start={},end={},ref="{}",alt="{}")'\
            .format(type(self).__name__,
                    self.chrom,self.start,self.end,
                    self.ref,self.alt)

    def edit_equal(self,other):
        """
        Returns True if ``self`` represents the same alteration as the ``other``
        """
        return self.vtp==other.vtp \
            and self.chrom==other.chrom\
            and self.start==other.start \
            and self.end==other.end \
            and self.alt==other.alt

    def __getattribute__(self,name):
        return object.__getattribute__(self,name)

    @classmethod
    def is_subclass(cls,other):
        """
        Checks is variant class is subclass of ``other``.  See module
        :mod:`genomvar.variant` documentaion for class hierarchy.
        """
        if issubclass(cls,other):
            return True
        else:
            return False

    def is_instance(self,cls):
        """Returns True if object's variant class is subclass of
        ``cls``.  See module :mod:`genomvar.variant` documentation for
        class hierarchy.
        """
        if issubclass(type(self),cls):
            return True
        else:
            return False

    def tolist(self):
        return [self.chrom,self.start,self.end,
                self.ref,self.alt,type(self)]
        
    @property
    def vtp(self):
        """Variant class"""
        return type(self)
        
    @property
    def key(self):
        """Tuple representation of variant."""
        try:
            return self._key
        except AttributeError:
            self._key = (self.chrom,self.start,self.end,
                         getattr(self,'ref',None),
                         getattr(self,'alt',None))
            return self._key

    @key.setter
    def key(self,value):
        raise NotImplementedError

    @key.deleter
    def key(self,value):
        delattr(self,'_key')


class MNP(VariantBase):
    """
    Multiple-nucleotide polymorphism.  Substitute N nucleotides of the
    reference for N other nucleotides.

    For instantiation it requires chromosome,position and 
    alternative sequence. ``end`` will inferred from ``start`` and ``alt``.
    ``ref`` is also optional.

    >>> from genomvar import variant
    >>> variant.MNP('chr1',154678,'GT')
    MNP("chr1",154678,"GT")
    
    """
    def __init__(self,chrom,start,alt,end=None,ref=None):
        if not end:
            end = start+len(alt)
        super().__init__(chrom,start,end,ref,alt)

    def __str__(self):
        return '<{:s} {:s}:{:d}-{:d} {:s}/{:s}>'\
            .format(type(self).__name__,self.chrom,self.start,self.end,
                    self.ref if self.ref else 'N'*(self.end-self.start),
                    self.alt if self.alt else '-')

    def __repr__(self):
        return 'MNP("{}",{},"{}")'\
            .format(self.chrom,self.start,self.alt)

    def nof_unit_vrt(self):
        return self.end-self.start

class SNP(MNP):
    """
    Single-nucleotide polymorphism.

    For instantiation it requires chromosome, position and 
    alternative sequence.

    >>> from genomvar import variant
    >>> variant.SNP('chr1',154678,'T')
    SNP("chr1",154678,"T")
    
    """
    def __init__(self,chrom,start,alt,end=None,ref=None):
        super().__init__(chrom=chrom,start=start,end=start+1,
                         ref=ref,alt=alt)
    def __repr__(self):
        return 'SNP("{}",{},"{}")'\
            .format(self.chrom,self.start,self.alt)

    def __str__(self):
        return '<{:s} {:s}:{:d} {:s}/{:s}>'\
            .format(type(self).__name__,self.chrom,self.start,
                    self.ref if self.ref else '-',
                    self.alt if self.alt else '-')

    def nof_unit_vrt(self):
        return 1

class Indel(VariantBase):
    """Abstract class to accomodate insertion/deletion subclasses"""
    def nof_unit_vrt(self):
        return 1


class Ins(Indel):
    """
    Insertion of nucleotides. For instantiation ``chrom``, 
    ``start``  and 
    inserted sequence (``alt``) are required.

    Start and end denote the nucleotide after the inserted sequence,
    i.e. ``start`` is 0-based number of a nucleotide after insertion,
    ``end`` is ``start+1`` by definition.
    
    >>> from genomvar.variant import Ins
    # Insertion of TA before position chr2:100543.
    >>> print(Ins('chr2',100543,'TA'))
    <Ins chr2:100543 -/TA>
    """
    def __init__(self,chrom,start,alt,end=None,ref=None):
        super().__init__(chrom,start,end=start+1,ref=ref,alt=alt)

    def __repr__(self):
        return 'Ins("{}",{},"{}")'\
            .format(self.chrom,self.start,self.alt)

    def __str__(self):
        return '<{:s} {:s}:{:d} -/{:s}>'\
            .format(type(self).__name__,self.chrom,self.start,
                    self.alt if self.alt else '-')

    def edit_equal(self,other,**kwds):
        return other.is_instance(Ins) \
            and self.chrom==other.chrom\
            and self.start==(other.act_start if other.is_instance(AmbigIns) \
                       else other.start)\
            and self.alt==other.alt

class Del(Indel):
    """
    Deletion of nucleotides. For instantiation ``chrom``,
    ``start`` (0-based),
    and ``end`` (position ``end`` is excluded) are required.
    
    >>> from genomvar.variant import Del
    # Deletion of 3 nucleotides starting at chr3:7843488 (0-based)
    >>> print(Del('chr3',7843488,7843488+3))
    <Del chr3:7843488-7843491 NNN/->
    """
    def __init__(self,chrom,start,end,ref=None,alt=None):
        super().__init__(chrom,start,end,ref,alt)

    def __repr__(self):
        return 'Del("{}",{},{})'\
            .format(self.chrom,self.start,self.end)

    def __str__(self):
        return '<{:s} {:s}:{:d}-{:d} {:s}/->'\
            .format(type(self).__name__,self.chrom,self.start,self.end,
                    self.ref if self.ref else 'N'*(self.end-self.start))

    def edit_equal(self,other,**kwds):
        v =  other.is_instance(Del) \
            and self.chrom==other.chrom\
            and self.start==(other.act_start if other.is_instance(AmbigDel) \
                       else other.start)\
            and self.end==(other.act_end if other.is_instance(AmbigDel) \
                       else other.end)
        return v

class AmbigIndel(Indel):
    """Class representing indel which position is ambigous.  Ambiguity means the
    indel could be applied in any position of some region resulting in the same
    alternative sequence.
    """
    def __str__(self):
        return '<{:s} {:s}:{:d}-{:d}({:d}-{:d}) {:s}/{:s}>'\
            .format(type(self).__name__,
                    self.chrom,self.start,self.end,
                    self.act_start,self.act_end,
                    self.ref if self.ref else '-',
                    self.alt if self.alt else '-')

    def shift_equal(self,other):
        a = (other.act_start-self.act_start) % len(self.seq)
        self_shift = self.seq[a:]+self.seq[:a]
        if self_shift==other.seq and self.chrom==other.chrom \
                 and self.vtp==other.vtp:
            return True
        else:
            return False

    def ambig_equal(self,other):
        return self.edit_equal(other) or self.shift_equal(other)

    def tolist(self):
        return [self.chrom,self.start,self.end,
                self.ref,self.alt,type(self),
                self.act_start,self.act_end]

class AmbigIns(AmbigIndel,Ins):
    """
    Class representing indel which position is ambigous.  Ambiguity
    means the indel could be applied in any position of some region
    resulting in the same alternative sequence.

    Let the reference file ``test.fasta`` contain a toy sequence::

      >seq1
      TTTAATA

    Consider a variant extending  3 ``T`` s in the beginning by one 
    more T. It can be done in several places so the corresponding
    insertion can be given as and ``AmbigIns`` object::

      >>> from genomvar import Reference
      >>> from genomvar.variant import VariantFactory
      >>> fac = VariantFactory(Reference('test.fasta'),normindel=True)
      >>> print( fac.from_edit('seq1',0,'T','TT') )
      <AmbigIns seq1:0-4(1-2) -/T>
    
    Positions 1 and 2 are actual start and end meaning that T
    is inserted before nucleotide located 1-2.
    Positions 0-4 indicate that start and end can be extended to 
    these values resulting in the same alteration.
    """
    def __init__(self,chrom,start,end,alt,ref=None):
        start,self.act_start = start
        end,self.act_end = end
        VariantBase.__init__(self,chrom=chrom,start=start,end=end,
                             ref=ref,alt=alt)
        
    def edit_equal(self,other):
        return other.is_instance(Ins) \
            and self.chrom==other.chrom\
            and self.act_start==(other.act_start \
                      if other.is_instance(AmbigIns) else other.start)\
            and self.alt==other.alt

    @property
    def seq(self):
        return self.alt

    def __repr__(self):
        return 'AmbigIns("{}",start=({},{}),end=({},{}),alt="{}")'\
            .format(self.chrom,self.start,self.act_start,self.end,
                    self.act_end,self.alt)

class AmbigDel(AmbigIndel,Del):
    """
    Class representing del which position is ambigous.  Ambiguity means
    the same number of positions could be deleted in some range
    resulting in the same aternative sequence.

    Let the reference file ``test.fasta`` contain a toy sequence::

        >seq1
        TCTTTTTGACTGG

    >>> fac = VariantFactory(Reference('test.fasta'),normindel=True) >>>
    print( fac.from_edit('seq1',1,'CTTTTTGAC','C') ) <AmbigDel
    seq1:1-11(2-10) TTTTTGAC/->

    Deletion of TTTTTGAC starts at 2 and ends on 9th nucleotide
    (including 9th resulting in range 2-10).  1-11 denote that start and
    end can be extended to these values resulting in the same
    alteration.
    """
    def __init__(self,chrom,start,end,ref=None,alt=None):
        start,self.act_start = start
        end,self.act_end = end
        VariantBase.__init__(self,chrom=chrom,start=start,end=end,
                             ref=ref,alt=alt)

    def __repr__(self):
        return 'AmbigDel("{}",start=({},{}),end=({},{}))'\
            .format(self.chrom,self.start,self.act_start,
                    self.end,self.act_end)

    def edit_equal(self,other):
        return other.is_instance(Del) \
            and self.chrom==other.chrom\
            and self.act_start==(other.act_start \
                      if other.is_instance(AmbigDel) else other.start)\
            and self.alt==other.alt

    @property
    def seq(self):
        return self.ref

class Mixed(VariantBase):
    """Combination of Indel and MNP. Usage of this class is discouraged
    and exists for compatibility."""
    def __init__(self,chrom,start,end,alt,ref=None):
        super().__init__(chrom=chrom,start=start,end=end,ref=ref,alt=alt)

    def __repr__(self):
        return '{}("{}",{},{},alt="{}")'\
            .format(type(self).__name__,
                    self.chrom,self.start,self.end,
                    self.alt)

    def nof_unit_vrt(self):
        return min(len(self.alt),self.end-self.start) + 1

class Haplotype(VariantBase):
    """
    An object representing genome variants on the
    same chromosome (or contig).

    Can be instantiated from a list of GenomVariant objects using
    :meth:`Haplotype.from_variants` class method.
    """
    def __init__(self,chrom,start,end,variants=None):
        super().__init__(chrom,start,ref='-',alt='-',end=end)
        self._variants = OrderedDict([(v.key,v) for v in variants]) \
                 if variants else OrderedDict()
        self.data = IntervalTree()
        for vrt in self._variants.values():
            self.data.add_interval(bxiv(vrt.start,vrt.end,vrt.key))

    def __str__(self):
        return '<Haplotype {:s}:{:d}-{:d} of {} variants>'\
            .format(self.chrom,self.start,
                    self.end,len(self._variants))
    def __repr__(self):
        return self.__str__()

    def nof_unit_vrt(self):
        return sum([v.nof_unit_vrt() for v in self.variants])

    def find_vrt(self,start=0,end=MAX_END):
        """
        Finds variants in specified region of a haplotype.
        
        Parameters
        ----------
        start : int, this description is getting a little bit 
            Start. *Default: 0*

            Here is longer description of this parameter
            which can involve elaborate details
        end : int
            End. *Default: end of haplotype*
        yields
        -------
        vrt : variants
        """
        for iv in self.data.find(start,end):
            yield self._variants[iv.value]

    @classmethod
    def from_variants(cls,variants):
        """
        Create haplotype from an iterable of variants

        Parameters
        ----------
        variants : list of variants
            Variants to instantiate haplotype from.
        Returns
        -------
        hap : Haplotype
            Haplotype variant object.
        """

        if len(variants)==0:
            raise ValueError('variants should be non-empty list')

        # infer chromosome
        chroms = set([o.chrom for o in variants])
        if len(chroms)>1:
            raise ValueError('variants on different chromosomes')
        elif len(chroms)==1:
            chrom = list(chroms)[0]
        else:
            raise ValueError('Chromosome not set')
        start = min([v.start for v in variants])
        end = max([v.end for v in variants])

        hap = cls.__new__(cls)
        hap.__init__(chrom,start,end,variants=variants)

        return hap

    @property
    def variants(self):
        for vrt in self._variants.values():
            yield vrt

class Null(VariantBase):
    def nof_unit_vrt(self):
        return 0

class Asterisk(VariantBase):
    def nof_unit_vrt(self):
        return 0

class GenomVariant(object):
    """
    This is a variant class holding gemomic alteration 
    (a VariantBase subclass object)
    and extra attributes. Underlying variant is accessible as ``base`` attribute.
    On top it has genotype as ``GT`` attribute and ``attrib`` containing
    extra attributes parsed from VCF files.  All attributes of
    underlying ``base`` (start, end etc.) are accessible from
    GenomVariant object.
    """
    def __init__(self,base,GT=None,attrib=None):
        if not isinstance(base,VariantBase):
            raise TypeError('Should be VariantBase or str')
        self.base = base
        self.GT = GT
        self.attrib = attrib if attrib else {}

    def __getattr__(self,name):
        if name=='base':
            raise AttributeError
        return object.__getattribute__(self.base,name)

    def __str__(self):
        return '{} GT:{}>'.\
            format(self.base.__str__()[:-1],
                   getattr(self,'GT',None))

    def __repr__(self):
        return 'GenomVariant({}, GT={})'.format(repr(self.base),str(self.GT))

    def edit_equal(self,other):
        """
        Check if GenomVariant holds the same genome alteration
        as ``other``.
        
        Parameters
        ----------
        variant : GenomVariant or any VariantBase
            variant to compare
        Returns
        -------
        equality : bool
            True if the same genomic alteration.
        """
        
        return self.base.edit_equal(other)

class VariantFactory(object):
    """
    Factory class used to create Variant objects.  Can be instantiated
    without arguments or an instance of :class:`genomvar.Reference` can
    be given. Additionnaly if ``normindel=True`` (reference should be
    given) indels will be checked for ambiguity, and if ambigous instance
    of :class:`genomvar.variant.AmbigIndel` will be returned.

    >>> from genomvar import Reference
    >>> from genomvar.variant import VariantFactory
    >>> reference = Reference('test.fasta')
    >>> fac = VariantFactory(reference,normindel=True)
    """
    def __init__(self,reference=None,normindel=False):
        if isinstance(reference,Reference):
            self.reference = reference
        elif isinstance(reference,str):
            self.reference = Reference(reference)
        elif reference is None:
            self.reference = None
        else:
            raise ValueError('reference not understood')
        if normindel:
            if self.reference is None:
                raise ValueError('reference needed for indel normalization')
        self.normindel = normindel

    def from_hgvs(self,st):
        """Returns a variant from a HGVS notation string. Only some of
        possible HGVS notations are supported as listed below::

          >>> fac = VariantFactory()
          >>> print(fac.from_hgvs('chr1:g.15C>A'))
          <SNP chr1:14 C/A>
          >>> print(fac.from_hgvs('chrW:g.19_21del'))
          <Del chrW:18-21 NNN/->
          >>> print(fac.from_hgvs('chrZ:g.10_11insCCT'))
          <Ins chrZ:10 -/CCT>
          >>> print(fac.from_hgvs('chr23:g.10delinsGA'))
          <Mixed chr23:9-10 -/GA>
          >>> print(fac.from_hgvs('chr24:g.145_147delinsTGG'))
          <MNP chr24:144-147 NNN/TGG>
        """
        chrom,rest = st.split(':',maxsplit=1)
        tp,posedit = rest.split('.',maxsplit=1)
        if '>' in posedit:
            match = re.match(hgvs_regex['SNP'],posedit)
            pos,ref,alt = match.groups()
            start = int(pos)-1
            return SNP(chrom=chrom,start=start,alt=alt,ref=ref)
        elif posedit.endswith('del'):
            match = re.match(hgvs_regex['Del'],posedit)
            pos1,pos2 = match.groups()
            start = int(pos1)-1
            end = int(pos2)
            if not self.normindel:
                return Del(chrom=chrom,start=start,end=end)
            else:
                start2,end2,ref,alt,cls = \
                      self._maybe_norm_del(chrom=chrom,start=start,end=end,
                                ref=self.reference.get(chrom,start,end))
                return cls(chrom,start2,end2,ref,alt)
        elif 'delins' in posedit:
            m = re.match(hgvs_regex['DelIns'][0],posedit)
            if m:
                pos1,pos2,ins = m.groups()
                start = int(pos1)-1
                end = int(pos2)
            else:
                m2 = re.match(hgvs_regex['DelIns'][1],posedit)
                if m2:
                    pos,ins = m2.groups()
                    start = int(pos)-1
                    end = start+1
                else:
                    raise ValueError('Unable to parse '+st)
            if end-start==len(ins): # MNP
                return MNP(chrom=chrom,start=start,end=end,alt=ins)
            else: #Mixed
                return Mixed(chrom=chrom,start=start,end=end,alt=ins)
        elif 'ins' in posedit:
            match = re.match(hgvs_regex['Ins'],posedit)
            pos,alt = match.groups()
            start = int(pos)-1
            if not self.normindel:
                return Ins(chrom,start=start,alt=alt)
            else:
                start2,end2,ref,alt,cls = \
                      self._maybe_norm_ins(chrom=chrom,start=start,
                                             alt=alt)
                return cls(chrom,start2,end2,ref,alt)
        else:
            raise ValueError('Unable to parse '+st)

    def from_edit(self,chrom,start,ref,alt):
        """
        Takes chrom, start, ref and alt (similar to a row in VCF file,
        but start is 0-based) and returns a variant object.

          >>> vrt = fac.from_edit('chr15',575,'TA','T')
          >>> print(vrt)
          <Del chr15:576-577 A/->
        
        Method attempts to strip ref and alt sequences from left and right
        (in this order), for example,
        
          >>> vrt = fac.from_edit('chr15',start=65,ref='CTG',alt='CTC')
          >>> print(vrt)
          <SNP chr15:67-68 G/C>
        """
        chrom,start,end,ref,alt,cls = \
               self._parse_edit(chrom,start,ref,alt)
        return cls(chrom=chrom,start=start,end=end,ref=ref,alt=alt)
    
    def _parse_edit(self,chrom,start,ref,alt):
        """
        Given a chrom, start, ref and alt record it returns a tuple
        for variant instantiation.
        """
        # store the start for later use
        start0 = start + min(len(ref),len(alt))

        if alt=='*':
            return (chrom,start,start+len(ref),ref,alt,Asterisk)

        if len(ref)==1 and len(alt)==1 and ref!=alt:
            tp=SNP
            end = start+len(ref)
        else:
            ref,alt,shift = _strip_ref_alt(ref,alt)
            start += shift
            if not ref: # Ins
                if not self.normindel:
                    return (chrom,start,start+1,'',alt,Ins)
                else:
                    start2,end2,ref,alt,cls = \
                          self._maybe_norm_ins(chrom=chrom,start=start,
                                            alt=alt)
                    return (chrom,start2,end2,ref,alt,cls)
            elif not alt: #Del
                if not self.normindel:
                    return (chrom,start,start+len(ref),ref,'',Del)
                else:
                    start2,end2,ref,alt,cls = \
                          self._maybe_norm_del(chrom=chrom,start=start,
                                               end=start+len(ref),ref=ref)
                    return (chrom,start2,end2,ref,alt,cls)
                
            elif len(ref)==len(alt): #MNP
                end = start+len(ref)
                tp = MNP if len(ref)>1 else SNP
            else: #Mixed
                return (chrom,start,start+len(ref),ref,alt,Mixed)

        return (chrom,start,end,ref,alt,tp)

    def _maybe_norm_ins(self,chrom,start,alt,max_len=MAX_END):
        """Check insertion for normalization"""
        left_norm = False
        right_norm = False
        # _start will be the first position to check for
        # equality with the unit
        _start = start
        it = itertools.cycle(reversed(alt))
        while _start >= 1:
            _start -= 1
            rl = str(self.reference.get(chrom,_start,_start+1))
            ul = next(it)
            if rl != ul:
                _start += 1
                break
            left_norm = True

        # Checking end
        it = itertools.cycle(alt)
        _end = start
        while _end < max_len:
            rl = str(self.reference.get(chrom,_end,_end+1))
            ul = next(it)
            if rl != ul:
                _end += 1
                break
            _end += 1
            right_norm = True

        if left_norm or right_norm:
            return ((_start,start),(_end,start+1),'',alt,AmbigIns)
        else:
            return (_start,_start+1,'',alt,Ins)

    def _maybe_norm_del(self,chrom,start,end,ref=None,max_len=MAX_END):
        """Check deletion if needs to be normalized"""

        if not ref:
            ref = self.reference(chrom,start,end)
        left_norm = False
        right_norm = False
        # _start will be the first position to check for
        # equality with the unit
        _start = start
        it = itertools.cycle(reversed(ref))
        while _start >= 0:
            _start -= 1
            rl = str(self.reference.get(chrom,_start,_start+1))
            ul = next(it)
            if rl != ul:
                _start += 1
                break
            left_norm = True

        # Checking end
        it = itertools.cycle(ref)
        _end = start + len(ref)
        while _end < max_len:
            rl = str(self.reference.get(chrom,_end,_end+1))
            ul = next(it)
            if rl != ul:
                break
            _end += 1
            right_norm = True

        if left_norm or right_norm:
            return ((_start,start),(_end,end),ref,'',AmbigDel)
        else:
            return (_start,_end,ref,'',Del)
