import numpy as np
from genomvar import variant

def _get_vartype(ref1,alt1):
    """ref1 and alt1 assumed to be stripped of identical nucleotides.
    Function returns Ins,Del,SNP or MNP"""
    if len(ref1)==0:
        return variant.Ins
    elif len(alt1)==0:
        return variant.Del
    else:
        if len(ref1)==len(alt1):
            if len(ref1)==1:
                return variant.SNP
            elif len(ref1)>1:
                return variant.MNP
            else:
                raise ValueError('Empty ref and alt given')
        else:
            return variant.Haplotype

def diff_slice(s1,s2):
    """
    >>> s1='ABCD'; s2='AXYD'; diff_slice(s1,s2)
    (1,3)
    >>> s1='TT'; s2='GT'; diff_slice(s1,s2)
    (0,1)
    """
    if len(s1) != len(s2):
        raise ValueError('Unequal lengths')
    start = 0
    end = len(s1)
    for pos in range(len(s1)):
        if s1[pos] == s2[pos]:
            start += 1
        else:
            break
    for pos in range(len(s1)):
        if s1[-(pos+1)] == s2[-(pos+1)]:
            end -= 1
        else:
            break
    if start > end:
        raise ValueError('Weird')
    return (start,end)

def ligate_mnp(mnps):
    """
    ligates a list of MNPs into one inserting '-'.
    mnps assumed to be non-overlapping.
    """
    _mnps = sorted(mnps,key=lambda o: o.start)
    prev = _mnps[0]
    alt = list(prev.alt)
    ref = list(prev.ref)

    for mnp in _mnps[1:]:
        alt += ['-']*(mnp.start-prev.end)
        alt += list(mnp.alt)
        ref += ['-']*(mnp.start-prev.end)
        ref += list(mnp.ref)
        prev = mnp

    chrom = _mnps[0].chrom
    start = _mnps[0].start
    end = _mnps[-1].end

    return  VariantBase(chrom=chrom,start=start,end=end,ref=''.join(ref),
                  alt=''.join(alt),vartype='mnp')

def novlp_grp(variants):
    """Returns non-overlapping groups of variants.
    Variants assumed to be sorted.

    For example,assuming the following variants::

      R         TTATATT
                -------
      v1,v3     AGT TG
      v2,v4      G    A

    >>> novlp_grp([v1,v2,v3,v4])
    [[MNP""<1:1-4 TTA->AGT GT=None>,SNP""<1:2-3 T->G GT=None>],
     [MNP""<1:5-7 AT->TG GT=None>],
     [SNP""<1:7-8 T->A GT=None>]]
    """
    if len(variants)==1:
        yield variants
    else:
        grp = []
        cend = -1
        start1 = variants[0].start
        for vrt in variants: #sorted(variants,key=lambda o: o.start):
            start2 = vrt.start
            if start2<start1:
                raise ValueError('Unsorted: {} seen after {}'\
                                 .format(start2,start1))
            if len(grp)>0:
                if start2 < cend:
                    if vrt.end > cend:
                        cend = vrt.end
                    grp.append(vrt)
                else:
                    yield grp
                    grp = [vrt]
                    cend = vrt.end
            else:
                grp = [vrt]
                cend = vrt.end
            start1 = start2
        yield grp

def novlp_cmbtns(variants):
    """Returns non-overlapping combinations of variants.
    Variants assumed to be sorted.

    For example,assuming the following variants::

      Ref       TTATATT
                -------
      vrt0      AGTCTGA
      vrt1,vrt3 AGT TG
      vrt2,vrt4  G    A

    >>> novlp_cmbtns([vrt1,vrt2,vrt3,vrt4])
    [[MNP""<1:1-4 TTA->AGT GT=None>,
         MNP""<1:5-7 AT->TG GT=None>,SNP""<1:7-8 T->A GT=None>],
     [SNP""<1:2-3 T->G GT=None>,
         MNP""<1:5-7 AT->TG GT=None>,SNP""<1:7-8 T->A GT=None>]]
    """

    if len(variants) == 1:
        return [variants]
    rt = []
    first = variants[0]
    first_end = first.end
    start1 = -1
    for vrt in variants:
        start2 = vrt.start
        if start2 < start1:
            raise ValueError('Should be sorted')
        # breaking when next cluster starts
        if start2>=first_end:
            break

        # now processing those recursively
        right = list(filter(lambda o: o.start>=first_end,variants))
        if len(right)==0:
            rt.append([vrt])
        elif len(right) == 1:
            rt.append([vrt,right[0]])
        else:
            cmb = novlp_cmbtns(right)
            for cmbtn in cmb:
                rt.append([vrt]+cmbtn)
        start1 = start2
    return rt

def _mnp_isect(mnp0,mnps):
    """
    This return com positions between mnp0 and mnps.
    mnps assumed to be non-overlapping
    """

    isect = ['-']*len(mnp0.alt)

    for mnp in mnps:
        _start = max(mnp0.start,mnp.start)
        _end = min(mnp0.end,mnp.end)
        #print('_start','_end',_start,_end)
        isect[_start-mnp0.start:_end-mnp0.start] = \
            mnp.alt[_start-mnp.start:_end-mnp.start]
    #print('isect',isect)

    rt = ['-']*len(mnp0.alt)
    for i in range(len(isect)):
        if i>=len(rt):
            break
        if isect[i] == mnp0.alt[i]:
            rt[i] = isect[i]
    return rt

def _cmp_mnps(vrt,variants,action):
    """
    # R         TTATATT
    #           -------
    # v0        AGTCTGA
    # v1,v3     AGT TG
    # v2,v4      G    A

    >>> _cmp_mnps(v0,[v1,v2,v3,v4],'comm')
    [SNP""<1:1-4 TTA->AGT GT=None>,SNP""<1:5-8 ATT->TGA GT=None>]

    >>> _cmp_mnps(v0,[v1,v2,v3,v4],'diff')
    [SNP""<1:4-5 T->C GT=None>]
    """
    alt = ['-']*len(vrt.alt)
    _key = lambda cmbtn: len(list(filter(lambda e: e!='-',
                                         _mnp_isect(vrt,cmbtn))))
    for grp in novlp_grp(variants):
        ovlp = _mnp_isect(vrt,grp)
        #print('ovlp',ovlp)
        for ind,nucl in enumerate(ovlp):
            if nucl!='-':
                alt[ind] = nucl

    #print('alt',alt)

    # Find contigous blocks depending on action
    blocks = []
    cur_start = 0
    pos = 0
    if action == 'comm':
        func = lambda l: l!='-'
    else:
        func = lambda l: l=='-'
    while True:
        if func(alt[pos]):
            cur_start = pos
            while pos<len(alt) and func(alt[pos]):
                pos += 1
            blocks.append( (cur_start,pos) )
            pos+=1
        else:
            pos+= 1
        if pos >= len(alt):
            break

    rt = []
    for block in blocks:
        cls = variant.SNP if (block[1]-block[0])==1 else variant.MNP
        _vrt = cls(chrom=vrt.chrom,
                   start=vrt.start+block[0],end=vrt.start+block[1],
                   ref=''.join(vrt.ref[block[0]:block[1]]) if vrt.ref \
                   else None,
                   alt=''.join(vrt.alt[block[0]:block[1]]))
        rt.append(_vrt)

    return rt

def apply_vrt(s,vrt,start=None):
    if not start:
        start = 0
    _s = dict(enumerate(s))
    for ind in range(vrt.start-start,vrt.end-start):
        _s.pop(ind)
    _s[vrt.start-start] = vrt.alt
    r = [_s[i] for i in sorted(_s)]
    return ''.join(r)

def gt_ovlp(gt1,gt2):
    return tuple([all(t) for t in zip(gt1,gt2)])

def are_complementary(*gt):
    shape = len(gt[0])
    if sum([any(t) for t in zip(*gt)]) == shape:
        return True
    else:
        return False

def nof_snp_vrt(mnps):
    rng = [min([o.start for o in mnps]),
           max([o.end for o in mnps])]

    dt = np.zeros(shape=(len(mnps),rng[1]-rng[0]),dtype='S1')

    for ind,vrt in enumerate(mnps):
        dt[ind,vrt.start-rng[0]:vrt.end-rng[0]] \
            = list(vrt.alt)

    cnt = 0
    for i in range(dt.shape[1]):
        unq = np.unique(dt[:,i])
        cnt += len(unq)
    return cnt
    

def compl_blk(*variants):
    rng = [min([o.start for o in variants]),max([o.end for o in variants])]

    ref = np.zeros(shape=(rng[1]-rng[0],),dtype='U1')
    dt = np.zeros(shape=(len(variants[0].GT),rng[1]-rng[0]),dtype='U1')


    for ind,vrt in enumerate(variants):
        dt[np.array(vrt.GT)==1,vrt.start-rng[0]:vrt.end-rng[0]] \
            = list(vrt.alt)
        ref[vrt.start-rng[0]:vrt.end-rng[0]] = list(vrt.ref)

    blocks = []
    left = 0
    right = 0
    while right < dt.shape[1]:
        unq = np.unique(dt[:,right])
        if not (unq.shape==(1,) and unq[0]!='' and are_complementary([o.GT for o in variants])):
            if left != right:
                blocks.append( (left,right) )
                right += 1
                left = right
            else:
                right += 1
                left += 1
        else:
            right += 1
    if left != right:
        blocks.append( (left,right) )

    ret_variants = []
    for blk in blocks:
        alt = ''.join(dt[0,blk[0]:blk[1]])
        ref = ''.join(ref[blk[0]:blk[1]])
        vrt = VariantBase(chrom=variants[0].chrom,start=blk[0]+rng[0],end=blk[1]+rng[0],
                     ref=ref,alt=alt,vartype='mnp')
        ret_variants.append(vrt)

    return ret_variants

def _matchv2(v1,v2,match_partial=True,match_ambig=False):
    """Take v1 (non-compound) and v2 (maybe Compound) and returns a 
    list of [(sub-v2,score),...]"""
    
    if v2.is_instance(variant.Haplotype):
        ovlp = v2.find_vrt(v1.start,v1.end)
    else:
        ovlp = [v2]
    m = []
    if v1.is_instance(variant.MNP):
        # vrt holds matching MNPs
        if match_partial:
            vrt = list(filter(lambda o: o.is_instance(variant.MNP),ovlp))
            s = _mnp_isect(v1,vrt)
            if len(vrt)>0:
                if v1.is_instance(variant.SNP):
                    assert len(vrt)==1
                score = len(s)-s.count('-')
                if score>0:
                    m.append( (vrt,score) )
        else:
            vrt = list(filter(
                lambda o: o.is_instance(variant.MNP) and \
                o.start==v1.start and o.end==v1.end \
                and o.alt==v1.alt,ovlp))
            if vrt:
                m.append( (vrt[0],v1.end-v1.start) )
    elif v1.is_instance(variant.AmbigIndel):
        if match_ambig:
            vrt = list(filter(lambda o: v1.ambig_equal(o),ovlp))
        else:
            vrt = list(filter(lambda o: v1.edit_equal(o),ovlp))
        if len(vrt)>0:
            m.append( (vrt,1) )
    elif v1.is_instance(variant.Indel):
        vrt = list(filter(
            lambda o: v1.edit_equal(o),ovlp))
        if len(vrt)>0:
            m.append( (vrt,1) )
    return m

def matchv(target,locus,match_partial=True,match_ambig=False):
    """
    Returns a dictionary with the best VariantBase IDs from locu
    matching variant target. target can be a haplotype.

    Dict is structured like this:
    {vrt_id (from vrt):[vrt_self]}

    In the return dictionary all vrt result from traversing the
    haplotypes if any encountered.
    """
    
    if len(locus)==0:
        return {}
    if not target.is_instance(variant.Haplotype):
        if target.is_instance(variant.MNP):
            tot_dt = []
            for grp in novlp_grp(locus):
                grp_dt = []
                grp_score = -1
                for cmb in novlp_cmbtns(grp):
                    cmb_score = 0
                    _dt = []
                    for locusv in cmb:
                        m = _matchv2(target,locusv,match_ambig=match_ambig,
                                     match_partial=match_partial)
                        # m is {target_key:([vrt_key from locus],score)}
                        for locus_vrts,vrt_score in m:
                            cmb_score += vrt_score
                            _dt.extend(locus_vrts)
                    if cmb_score > grp_score:
                        grp_dt = _dt
                        grp_score = cmb_score
                tot_dt.extend(grp_dt)
            return tot_dt
        elif target.is_instance(variant.Indel):
            # this could be done using _matchv2 but dealing
            # with it explicitely for performance
            if target.is_instance(variant.AmbigIndel) and match_ambig:
                return list(filter(lambda v: target.ambig_equal(v),locus))
            else:
                return list(filter(lambda v: target.edit_equal(v),locus))
        else:
            return list(filter(lambda v:v.edit_equal(target),locus))
    else:
        tot_dt = {}
        for grp in novlp_grp(locus):
            grp_dt = {}
            grp_score = -1
            for cmb in novlp_cmbtns(grp):
                cmb_score = 0
                _dt = {}
                for locusv in cmb:
                    m = {}
                    for vrt in target.find_vrt(locusv.start,locusv.end):
                        sm = _matchv2(vrt,locusv,match_ambig=match_ambig,
                                     match_partial=match_partial)
                        if sm:
                            m[vrt.key] = sm
                    # m is {target_key:([vrt_key from locus],score)}
                    for target_key in m:
                        for locus_vrts,vrt_score in m[target_key]:
                            cmb_score += vrt_score
                            _dt.setdefault(target_key,[]).extend(locus_vrts)
                if cmb_score > grp_score:
                    grp_dt = _dt
                    grp_score = cmb_score

            for target_key,vals in grp_dt.items():
                if target_key in tot_dt:
                    tot_dt[target_key] += vals
                else:
                    tot_dt[target_key] = vals
        return tot_dt

def _cmpv2(target,match,action,callback=None):
    """Takes non-Compound target and a list of matching and returns
    a comparison result"""
    if not match:
        raise ValueError('no match')
    if target.is_instance(variant.MNP):
        vrt2add = []
        if isinstance(target,variant.GenomVariant):
            for vrt in _cmp_mnps(target,match,action=action):
                _vrt = variant.GenomVariant(vrt,
                                attrib=getattr(target,'attrib',None))
                if callback:
                    _vrt.attrib['cmp'] = callback(match)
                vrt2add.append(_vrt)
        else:
            vrt2add.extend(_cmp_mnps(target,match,action=action))
        return vrt2add
    else: # target.is_instance(variant.Indel):
        if action=='diff':
            return [] # Generally should not be here because cmpv
                      # returned already
        elif action=='comm':
            if callback:
                cv = callback(match)
                try:
                    target.attrib['cmp'] = cv
                except AttributeError as exc:
                    if isinstance(target,variant.VariantBase):
                        _vrt = variant.GenomVariant(target)
                        _vrt.attrib['cmp'] = cv
                        return _vrt
                    else:
                        raise exc
                return [target]
            else:
                return [target]
    
def cmpv(target,match,action,callback=None):
    """
    Get result of comparing target and match. 

    Match is expected to be formatted like output of method
    :meth:`varset.VariantSet.match`, i.e. a dictionary if 
    target is a haplotype and a list otherwise.

    Parameters
    ----------
    target : VariantBase subclass or GenomVariant
        A variant to consider.
    match : list or dict
        Matching variants to find common or subtract.
    action : {'diff','comm'}
        Action to perform.
    Returns
    -------
    cmp_result : list of :class:`variant.GenomVariant`.
    """
    if not target.is_instance(variant.Haplotype):
        if not match:
            if action=='diff':
                return [target]
            else:
                return []
        else:
            return _cmpv2(target,match,action=action,callback=callback)
    else:
        vrt2add = []
        for target_vrt in target.variants:
            if not target_vrt.key in match:
                if action=='diff':
                    vrt2add.append(target_vrt)
                elif action=='comm':
                    continue
            else:
                vrt2add.extend(_cmpv2(target_vrt,match[target_vrt.key],
                                      action=action,callback=callback))
        return vrt2add


