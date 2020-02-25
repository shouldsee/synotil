#!/usr/bin/env python
from distutils.core import setup

setup(name='synotil',
      version='1.0',
      description='Utility on Synology datasets',
      author='Feng Geng',
      author_email='shouldsee.gem@gmail.com',
      url='none',
      install_requires=[
      x.strip() for x in open("requirements.txt","r")
      if x.strip() and not x.strip().startswith("#")
      ],
      
      package_dir = {'synotil': 'synotil'},
      packages=['synotil', ],
     )
