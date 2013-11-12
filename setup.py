from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc
import numpy as np


keporb_sources = [ 'planetc/keporb.c', 'planetc/keporb_src.c' ]
ma02_sources = [ 'planetc/ma02.c', 'planetc/ma02_quad_src.c', \
                 'planetc/ma02_nonlin_src.c', 'planetc/ma02_utils_src.c', \
                 'planetc/hyp2f1.c', 'planetc/hypappell.c', 'planetc/const.c', \
                 'planetc/fabs.c', 'planetc/gamma.c', 'planetc/mtherr.c', \
                 'planetc/polevl.c', 'planetc/psi.c', 'planetc/round.c' ]

python_inc = get_python_inc() # python headers
numpy_inc = np.get_include() # numpy headers

setup(

    name = 'planetc',
    version='0.0.1', 
    author = 'Tom Evans',
    author_email = 'tom.evans@astro.ox.ac.uk',
    description = 'The planetc package contains code for implementing the analytic Mandel & Agol (2002) equationsfor planetary transit lightcurves. To do this, it uses a C backend with Python/Cython wrappers. Routines for calculating radial velocities of Keplerian orbits will be added next.',
    packages = [ 'planetc' ],
    ext_modules =
    [
        Extension(
        'planetc.keporb', sources=keporb_sources,
        include_dirs=[ python_inc, numpy_inc ] ), 
        Extension(
        'planetc.ma02', sources = ma02_sources,
        libraries = [ 'gsl', 'gslcblas' ], 
        include_dirs = [ python_inc, numpy_inc ] )
    ],
    data_files = [ ( '', [ 'README.md', 'INSTALL', 'MANIFEST.in', 'Makefile' ] ) ]
    
    )

