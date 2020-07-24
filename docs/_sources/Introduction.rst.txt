Introduction / Quick start
##########################

Package ``genomvar`` works with genomic variants and implements
set-like operations on them. It supports import from VCF files and
export to NumPy.

For documentation see `here <https://mikpom.github.io/genomvar/>`_.

Installation
============

Requirements:

1. Python >=3.6
2. rbi-tree
3. jinja2
4. pysam

To install::

  pip install genomvar

Sample usage
============

Case 1
------

Common task in genome variant analysis is: there are two VCF files (for
example obtained from variant caller #1 and caller #2)
and the differences should be investigated.

First we read the VCF files
into genomvar :class:`genomvar.varset.VariantSet` objects which 
hold the variants with underlying data contained in INFO fields:

.. code-block:: python

  >>> from genomvar.varset import VariantSet
  >>> vs1 = VariantSet.from_vcf('caller1.out.vcf.gz',parse_info=True)
  >>> vs2 = VariantSet.from_vcf('caller2.out.vcf.gz',parse_info=True)

To find variants detected by caller #1 but not caller #2 ``diff``
method is used. Then differences are exported to ``numpy`` for futher
analysis:

.. code-block:: python

  >>> diff = vs1.diff(vs2)
  >>> recs = diff.to_records() # recs is a numpy structured dtype array
  >>> recs[['chrom','start','end','ref','alt','vartype']]
  [('chr1',  1046755,  1046756, 'T', 'G', 'SNP')
   ('chr1',  1057987,  1057988, 'T', 'C', 'SNP')
    ...,
   ('chr19', 56434340, 56434341, 'A', 'G', 'SNP')
   ('chrY', 56839067, 56839068, 'A', 'G', 'SNP')]
  >>> recs['INFO']['DP'].mean() # recs['INFO']['DP'] is a numpy ndarray
  232.18819746028257

Case 2
------

There is a smaller variant file obtained from the data and a bigger one
usually obtained from a database. Variants in the former should be "annotated"
with some data associated with variants in the latter.

This case is different from the previous in that DB file might not
comfortably fit into memory. Class
:class:`~genomvar.varset.IndexedVariantFileSet` can be used for this
purpose:

.. code-block:: python

    >>> vs = varset.VariantSet.from_vcf('vcf_of_interest.vcf')
    >>> dbSNP = varset.IndexedVariantFileSet('DBSNP.vcf.gz')
    >>> annots = []
    >>> for vrt in vs.iter_vrt():
    >>>     m = dbSNP.match(vrt)
    >>>     annots.append(m[0].attrib['id'] if m else None)
    >>> annots
    [None, None, 'rs540057607', 'rs367710686', 'rs940651103', ...]


Here :meth:`~genomvar.varset.VariantSet.match` method is used. It
searches for variants with the same genomic alteration as argument
variant and returns a list of those.  Then VCF ``ID`` field can be
accessed from those matching variants in ``attrib['id']`` (dbSNP rs
numbers in this particular case).
