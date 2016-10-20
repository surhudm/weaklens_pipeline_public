#!/usr/bin/env python

"""
setup.py file for wpipe
"""

from distutils.core import setup, Extension
import numpy

weaklens_module = Extension('_weaklens',
                           sources=['src/weaklens.i', 'src/weaklens.cpp', 'src/cosmology.cpp', 'src/haloes.cpp', 'src/powerspectrum.cpp', 'src/gauleg.cpp', 'src/kdtree2.cpp'],
                           swig_opts=["-c++"],
                           libraries=['m','gsl','gslcblas', 'boost_system'],
                           include_dirs=[numpy.get_include()],
                           )

tomography_module = Extension('_tomography',
                           sources=['src/tomography.i', 'src/tomography.cpp', 'src/cosmology.cpp', 'src/haloes.cpp', 'src/powerspectrum.cpp', 'src/gauleg.cpp', 'src/kdtree2.cpp'],
                           swig_opts=["-c++"],
                           libraries=['m','gsl','gslcblas', 'boost_system'],
                           include_dirs=[numpy.get_include()],
                           )

setup (name        = 'wpipe',
       version     = '0.01',
       author      = "Surhud More",
       url         = "http://member.ipmu.jp/surhud.more/research",
       author_email= "surhud.more@ipmu.jp",
       description = """Weak lensing pipeline tools""",
       ext_modules = [weaklens_module, tomography_module],
       license     = ['GPL'],
       py_modules  = ["weaklens", "tomography"],
       package_dir = { '':'src'},
       )
