import re
import itertools
from collections import deque
import heapq

def chunkit(it,size):
    cur_chunk = []
    for el in it:
        cur_chunk.append(el)
        if len(cur_chunk)>=size:
            yield cur_chunk
            cur_chunk = []
    yield cur_chunk
    return
    
def grouper(it,n,fill=None):
    args = [it]*n
    return itertools.zip_longest(*args,fillvalue=fill)

class zip_variants():
    class Exhasted(Exception):
        pass

    def __init__(self,left,right,direction='both'):
        self.deque = [deque(),deque()]
        self.cache = deque()
        _left = zip(itertools.repeat(0),left)
        _right = zip(itertools.repeat(1),right)
        self.merged = heapq.merge(_left,_right,key=lambda i:i[1].start)

    def __iter__(self):
        return self

    def _trim_deques(self,val):
        for deque in self.deque:
            if not deque:
                continue
            while True:
                if deque and deque[0].end<val.start:
                    deque.popleft()
                else:
                    break

    def _next_(self):
        try:
            deque,val = self.cache.popleft()
        except IndexError:
            try:
                deque,val = next(self.merged)
                self.deque[deque].append(val)
            except StopIteration:
                raise self.Exhasted
        self._trim_deques(val)
        return deque,val

    def __next__(self):
        try:
            deque,vrt = self._next_()
        except self.Exhasted:
            raise StopIteration

        # ~deque+2 is merely a fast way to convert 0 to 1 and 1 to 0
        if self.deque[~deque+2] and \
                      self.deque[~deque+2][-1].start>vrt.end:
            v = (deque,vrt,list(self.deque[~deque+2]))
            return v
        else:
            while True:
                try:
                    deque2,vrt2 = next(self.merged)
                    self.deque[deque2].append(vrt2)
                except StopIteration:
                    return (deque,vrt,list(self.deque[~deque+2]))
                    break

                self.cache.append((deque2,vrt2))
                if vrt2.start>vrt.end:
                    v = (deque,vrt,list(self.deque[~deque+2]))
                    return v
                    break
                else:
                    continue

def _strip_ref_alt(ref,alt):
    # Check identical in the beginning
    i = 0
    mn = min(len(ref),len(alt))
    while i < mn:
        if ref[i] == alt[i]:
            i += 1
            continue
        else:
            break
    shift=i

    # Strip positions identical at the end
    j=1
    while j <= mn-i:
        if ref[-j] == alt[-j]:
            j += 1
            continue
        else:
            break
    j -= 1
    ref = ref[i:len(ref)-j]
    alt = alt[i:len(alt)-j]
    return ref,alt,shift

def rgn_from(s):
    m = re.match('(.+):([0-9,]+)-([0-9,]+)',s)
    if m:
        chrom,start,end = m.groups()
        return chrom,int(start.replace(',','')),\
            int(end.replace(',',''))
    else:
        raise ValueError('Could not parse region: {}'.format(s))

def no_ovlp(it):
    for chrom,chrom_it in itertools.groupby(it,lambda v: v.chrom):
        cur_chunk = []
        for vrt in chrom_it:
            if len(cur_chunk)>0:
                if vrt.start < max([i.end for i in cur_chunk]):
                    cur_chunk.append(vrt)
                else:
                    yield cur_chunk
                    cur_chunk = [vrt]
            else:
                cur_chunk = [vrt]
        if len(cur_chunk)>0:
            yield cur_chunk
        
