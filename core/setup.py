#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension

import sys
sys.path.append("/home/giorgi/soft/multilevel-daphne/core")

callblackscholes_module = Extension('_callblackscholes',
                           sources=['callblackscholes_wrap.cxx', 
                                    'callblackscholes.cpp', 
                                    'structuralparameters.cpp',
                                    'multilevelparameters.cpp',
                                    'functions.cpp'],
                           extra_compile_args=['-std=c++11'],
                           include_dirs = ['/usr/include/eigen']
                           )

setup (name = 'callblackscholes',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       swig_opts=['-c++'],
       ext_modules = [callblackscholes_module],
       py_modules = ["callblackscholes"],
       )

