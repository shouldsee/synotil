#!/usr/bin/env python
from distutils.core import setup

setup(name='synotil',
      version='1.0',
      description='Utility on Synology datasets',
      author='Feng Geng',
      author_email='shouldsee.gem@gmail.com',
      url='none',
      install_requires=[
            "biopython",
            #biograpy
            #hmmlearn
            "hmmlearn==0.2.2",
            "biograpy @ https://github.com/shouldsee/BioGraPy/tarball/2b31b12",
            #multiprocessing
            "pandas",
            "pyBigWig",
            "scipy",
            "sklearn",
            #textwrap
            "wrapt",
      ],
      
      package_dir = {'synotil': 'synotil'},
      packages=['synotil', ],
     )
