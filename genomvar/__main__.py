import sys
import argparse
from genomvar.cmdutils import _cmp_vcf
from genomvar import OverlappingHaplotypeVars,UnsortedVariantFileError

tasks = ['compare_vcf']
options=[{'arg':'--match_partial','action':'store_true',
          'help':'Flag; if set common positions will be searched in MNPs',
          'default':False},
         {'arg':'vcf1','metavar':'VCF1','type':str,
          'help':'First VCF file'},
         {'arg':'vcf2','metavar':'VCF2','type':str,
          'help':'Second VCF file'},
         {'arg':'--out','metavar':'OUT','type':str,
          'help':'Output VCF file (STDOUT by default)'},
]
prsr=argparse.ArgumentParser(prog='compare_vcf',
             description='Output VCF containing variants'\
                                 +' from two VCFs compared')
for opt in options:
    arg = opt.pop('arg')
    prsr.add_argument(arg,**opt)

if len(sys.argv)<2 or (sys.argv[1] not in tasks):
    usage = """Usage: python3 -m genomvar compare_vcf [options]
    > python3 -m genomvar compare_vcf -h (--help)
    """.format(*tasks)
    print(usage)
    sys.exit()

args = prsr.parse_args(sys.argv[2:])
if args.out is None:
    out = sys.stdout
    file_opened = False
else:
    out = open(args.out,'w')
    file_opened = True


try:
    cnt = _cmp_vcf(f1=args.vcf1,f2=args.vcf2,
             out=out,match_partial=args.match_partial)
except UnsortedVariantFileError as exc:
    print('Error: Unsorted Variant file detected'\
          +'\nMessage was: '+exc.args[0],
          file=sys.stderr)
print('Unit variants:\n'\
      +'     first only: {}\n'.format(cnt[0])\
      +'     second only: {}\n'.format(cnt[1])\
      +'     both: {}'.format(cnt[2]),
      file=sys.stderr)
if file_opened:
    out.close()
