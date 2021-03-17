"""
# genomvar #
Package allows to work with genomic variations in Python.
"""
import collections
import pyfaidx
import numpy as np

class OverlappingHaplotypeVars(Exception):
    def __init__(self,message,ovlp=None):
        self.message = message
        self.ovlp = ovlp
class VCFSampleMismatch(Exception):pass
class UnsortedVariantFileError(Exception): pass
class DifferentlySortedChromsError(Exception): pass
class DuplicateVariants(Exception): pass
class NoVCFNotationError(Exception): pass
class StructuralVariantError(NotImplementedError): pass
class NoIndexFoundError(Exception): pass

class ChromSet(collections.abc.Set):
    '''Chromosome set remembers the order of chromosome addition.'''
    def __init__(self, iterable=None):
        self.elements =  collections.OrderedDict()
        if not iterable is None:
            for value in iterable:
                self.elements[value] = 0

    def __iter__(self):
        return iter(self.elements)

    def __contains__(self, value):
        return value in self.elements

    def __len__(self):
        return len(self.elements)

    def add(self,value):
        self.elements[value] = 0

    def __repr__(self):
        return 'ChromSet({})'.format(repr(list(self.elements)))

class Reference(object):
    """A wrapper around pyfaidx FastaFile class.
    It supports cached fetching of reference sequence.
    Speed up is achieved by minimizing disk IO
    when consequtive queries are close to one another.
    """
    def __init__(self, filename, cache_dst=1000000):
        self.filename = filename
        self.REF = pyfaidx.Fasta(self.filename)
        self.ctg_len = {k:len(v) for k,v in self.REF.records.items()}
        self.cache_dst = cache_dst
        self.cache = {}

    def get(self,chrom,start,end):
        """Reference object uses cache to improve performance on consequtive
        locations by reducing IO."""
        if start<0:
            raise IndexError('start should be >=0')
        if end - start >= self.cache_dst: # Large sequence is requested
            return str(self.REF[chrsom][start:end])
        c = self.cache.get(chrom,None)
        if c: # if chrom in cache, try it
            start2 = start - c[0] # c[0] is a start of seq in cache
            end2 = end - c[0]
            if start2>=0 and end2<2*self.cache_dst-1: # within cached
                ret = c[1][start2:end2]
                return ret
            else: # new location on chrom, new cache for chrom
                cstart = max(start-self.cache_dst,0)
                cend = cstart+2*self.cache_dst
                refseq = str(self.REF[chrom][cstart:cend])
                self.cache[chrom] = (cstart,refseq)
                return refseq[start-cstart:end-cstart]
        else:
            cstart = max(start-self.cache_dst,0)
            cend = cstart+2*self.cache_dst
            refseq = str(self.REF[chrom][cstart:cend])
            self.cache[chrom] = (cstart,refseq)
            return refseq[start-cstart:end-cstart]

    def close(self):
        self.REF.close()

    def get_chroms(self):
        return list(self.REF.records)

SINGLETON = -1

MAX_END = np.iinfo(np.int32).max # Looks like it is max for bx-python
