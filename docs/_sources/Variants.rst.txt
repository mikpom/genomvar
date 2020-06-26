Variants
********

.. automodule:: genomvar.variant
   :noindex:
      
Variant classes
===============

Basic classes are SNP representing a single nucleotide
polimorphism and its generalization MNP – multiple nucleotide
polimorphism.

.. autoclass:: genomvar.variant.SNP
   :noindex:

.. autoclass:: genomvar.variant.MNP
   :noindex:
      
There are separate classes of insetion and deletion inheriting from
abstract class ``Indel``. 

.. autoclass:: genomvar.variant.Ins
   :noindex:
      
.. autoclass:: genomvar.variant.Del
   :noindex:
      
There is a special flavor of indels for cases when it can be
applied in several places resulting in the same alternate sequence,
termed to as  ``ambigous`` indels. They are :class:`AmbigDel`,
:class:`AmbigIns` which on top of regular deletion or insertion
attributes 
contain information about a region they can be applied to. For 
instantion :class:`VariantFactory` with a reference is needed.

.. autoclass:: genomvar.variant.AmbigIns
   :noindex:
      
.. autoclass:: genomvar.variant.AmbigDel   
   :noindex:               

There is a separate class ``Haplotype`` which can hold any combination
of objects of variant classes above.

.. autoclass:: genomvar.variant.Haplotype
   :noindex:
      
Additionally module contains two technically driven types: 
Null (no variant at all, reference), 
Asterisk (for * in the ALT field of VCF files).

