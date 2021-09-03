Variant sets
************
.. automodule:: genomvar.varset

Important operation on variant sets is their comparison, answering a
common question of finding variants present in one set and absent in
the other. In-memory variant sets support
:meth:`~genomvar.varset.VariantSet.diff` and
:meth:`~genomvar.varset.VariantSet.comm` methods returning in-memory
variant sets of corrresponding type::

  >>> fac = VariantFactory()
  >>> mnp = fac.from_edit('1',100000,'AG','TT')
  >>> vs1 = varset.VariantSet.from_variants([mnp])

Now the same two variants but split::

  >>> snp1 = fac.from_edit('1',100000,'A','T')
  >>> snp2 = fac.from_edit('1',100001,'G','T')
  >>> vs2 = varset.VariantSet.from_variants([snp1,snp2])

Evaluate the difference using ``comm`` method::
  
  >>> list(vs1.comm(vs2,match_partial=True).iter_vrt())
  [GenomVariant(MNP("1",100000,"TT"), GT=None)]
  >>> list(vs1.comm(vs2,match_partial=False).iter_vrt())
  []

Note that when ``match_partial`` is ``False`` there are no common
variants.

Alternatively methods ``diff_vrt`` and ``comm_vrt`` avaible for all
variant set classes and yield variants. These methods return
intermediate comparison object of class
:class:`~genomvar.varset.CmpSet`, which can be further inspected. Lets assume
some file ``test.vcf`` contains only the two lines corresponding to
variants in ``vs2`` above::

  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
  1	100001	.	A	T	100	PASS	DP=10
  1	100002	.	G	T	100	PASS	DP=10

Lets create a toy variant set wrapper for this file and compare it
to ``vs1``::

  >>> vs3 = varset.VariantFileSet('test.vcf')
  >>> for vrt in vs3.comm_vrt(vs1).iter_vrt():
  ...     print(vrt)
  <SNP 1:100000 A/T GT:None>
  <SNP 1:100001 G/T GT:None>

This does not require any indexes because ``test.vcf`` is assumed to
be sorted by coordinate.
