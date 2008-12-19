#!/usr/bin/env python
from setuptools import setup, Extension

# Data files for simulation
data = ['wing.stl', 'body.stl', 'wing_info.txt', 'yan_1.txt', 'default_config.txt']
package_data = {'fmech.data':data}

setup( name='fmech',
       version='0.2',
       description='GUF Flight Mechanics Model',
       author='Will Dickson',
       author_email='wbd@caltech.edu',
       packages=['fmech', 'fmech.data', 'stl_tools'],
       package_data = package_data,
       ext_modules=[Extension("fmech.mirtich", ["src/mirtich.pyx","src/volInt.c"],)],
       )
