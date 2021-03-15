from unittest import TestCase
from genomvar import Reference
from pkg_resources import resource_filename as pkg_file
from genomvar import variant

class MyTestCase(TestCase):
    def setUp(self):
        self.chr15rgn = Reference(
            pkg_file(__name__,'data/chr15rgn.fna'))
        self.chr1rgn = Reference(
            pkg_file(__name__,'data/chr1rgn.fasta'))
        self.nvf = variant.VariantFactory(self.chr15rgn,normindel=True)
        # Factory not normalizing indels
        self.svf = variant.VariantFactory()

    def tearDown(self):
        self.chr15rgn.close()
        self.chr1rgn.close()
