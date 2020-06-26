Comparing two VCF files from command line
****************************************

For convenience a common operation of VCF comparison can
be done from command line like this::

  $ python3 -m genomvar compare_vcf vcf1.vcf vcf2.vcf > comparison.vcf
  Unit variants:
     first only: 1
     second only: 2
     both: 105

The numbers can be used for plotting Venn diagram and indicate
corresponding patch sizes.

File ``comparison.vcf`` will contain something like this::

  chr1    18353435        .       C       G       100     .       mt=SNP;whichVCF=both;ln=186;ln2=374
  chr1    18704102        .       T       A       100     .       mt=SNP;whichVCF=both;ln=187;ln2=375

INFO fields indicate which VCF contains the variant on the line. ``ln`` and ``ln2`` indicate line numbers in original files.

  
