#!/usr/bin/env python

"""
setup.py file for wpipe
"""

from distutils.core import setup, Extension

weaklens_module = Extension('_weaklens',
                           sources=['src/weaklens.i', 'src/weaklens.cpp', 'src/cosmology.cpp', 'src/haloes.cpp', 'src/powerspectrum.cpp', 'src/gauleg.cpp', 'src/kdtree2.cpp'],
                           swig_opts=["-c++"],
                           libraries=['m','gsl','gslcblas', 'boost_system'],
                           )

setup (name        = 'wpipe',
       version     = '0.01',
       author      = "Surhud More",
       url         = "http://member.ipmu.jp/surhud.more/research",
       author_email= "surhud.more@ipmu.jp",
       description = """Weak lensing pipeline tools""",
       ext_modules = [weaklens_module],
       license     = ['GPL'],
       py_modules  = ["weaklens"],
       package_dir = { '':'src'},
       )
