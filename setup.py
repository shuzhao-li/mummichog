#!/usr/bin/env python

from setuptools import setup

setup(
  name='mummichog',
  version='2.1.0',

  author='Shuzhao Li, Andrei Todor',
  author_email='shuzhao.li@gmail.com',
  description='Pathway and network analysis for metabolomics data',
  long_description=open('README.rst').read(),
  url='http://mummichog.org',
  license='BSD 3-Clause',


  keywords='metabolomics analysis bioinformatics mass spectrometry systems biology',

  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 2.7',
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

  install_requires=[
    'matplotlib',
    'networkx>=1,<2',
    'numpy',
    'scipy',
    'xlsxwriter',
  ],

)