from setuptools import setup
import re

def readme():
    with open('README.rst') as f:
        desc = f.read()
        desc = re.sub(':class:`~?(.+)`',r'\1',desc)
        desc = re.sub(':meth:`~?(.+)`',r'\1',desc)
        return desc

setup(name='genomvar',
      version='0.3.0',
      description='Sequence variant analysis in Python',
      long_description=readme(),
      classifiers=[
          'Intended Audience :: Science/Research',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries'],
      keywords= ['Bioinformatics', 'Computational Biology'],
      url='http://relativity.nprog.local/mikpom/genomvar',
      author='Mikhail Pomaznoy',
      author_email='mikpom@mailbox.org',
      license='BSD',
      packages=['genomvar'],
      package_data={
          'genomvar':['tmpl/*']
      },
      install_requires=[
          'rbi-tree>=0.1.3.1',
          'jinja2'],
      test_suite='unittest',
      tests_require=['setuptools'],
      zip_safe=False)
