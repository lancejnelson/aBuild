#!/usr/bin/env python
try:
    from setuptools import setup
    args = {}
except ImportError:
    from distutils.core import setup
    print("""\
*** WARNING: setuptools is not found.  Using distutils...
""")

from setuptools import setup
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    read_md = lambda f: open(f, 'r').read()

from os import path
setup(name='aBuild',
      version='1.0',
      description='Automator for f.p. calcs and fitting.',
      long_description= "" if not path.isfile("README.md") else read_md('README.md'),
      author='Lance J Nelson',
      author_email='lancejnelson@gmail.com',
      url='https://github.com/lancejnelson/aBuild',
      license='MIT',
      setup_requires=['pytest-runner',],
      tests_require=['pytest', 'numpy', 'phonopy'],
      install_requires=[
          "argparse",
          "ase",
          "pyparsing",
          "termcolor",
          "six",
          "numpy",
          "phonopy",
          "requests",
          "beautifulsoup4",
          "tqdm",
          "html5lib",
          "mpld3",
          "phenum",
          "h5py",
          "lazy_import",
          "seekpath"
      ],
      packages=['aBuild', 'aBuild.database', 'aBuild.fitting','aBuild.calculators'],
      scripts=['aBuild/scripts/builder.py'],
      package_data={'aBuild': []},
      include_package_data=True,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          'Natural Language :: English',
          'Operating System :: MacOS',
          'Operating System :: Microsoft :: Windows',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering',
      ],
     )
