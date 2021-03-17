from unittest import TestCase
from genomvar import Reference
from pkg_resources import resource_filename as pkg_file
from genomvar import variant

class MyTestCase(TestCase):
    def setUp(self):
        self.chr24 = Reference(
            pkg_file(__name__,'data/chr24.fna'))
        self.nvf = variant.VariantFactory(self.chr24,normindel=True)
        # Factory not normalizing indels
        self.svf = variant.VariantFactory()

    def tearDown(self):
        self.chr24.close()
