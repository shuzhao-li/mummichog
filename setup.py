#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
  name='mummichog',
  version='2.3.2',

  author='Shuzhao Li, Andrei Todor',
  author_email='shuzhao.li@gmail.com',
  description='Pathway and network analysis for metabolomics data',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='https://github.com/shuzhao-li/mummichog',
  license='BSD 3-Clause',
  keywords='metabolomics analysis bioinformatics mass spectrometry systems biology',

  classifiers=[
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  packages=['mummichog', 'mummichog/tests'],
  include_package_data=True,
  zip_safe=True,
  entry_points = {
        'console_scripts': ['mummichog=mummichog.command_line:main'],
    },

  python_requires='>=3.4',
  install_requires=[
    'matplotlib',
    'networkx>=1,<2',
    'numpy',
    'scipy',
    'xlsxwriter',
  ],

)
