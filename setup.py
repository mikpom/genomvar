from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='genomvar',
      version='0.1.13',
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
          'bx-python',
          'jinja2',
          'pysam'],
      test_suite='unittest',
      tests_require=['setuptools'],
      zip_safe=False)
