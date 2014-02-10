"""gridmsm: Markov state models with a grid discretization

More description
"""
from __future__ import print_function, division
import sys
from distutils.core import setup, Extension
import numpy as np
from numpy.distutils import system_info
from Cython.Distutils import build_ext
DOCLINES = __doc__.split("\n")

lapack_info = system_info.get_info('lapack')
if len(lapack_info) == 0:
    print('LAPACK libraries could not be located.', file=sys.stderr)
    sys.exit(1)

map_transmat = Extension('gridmsm._map_transmat',
    language='c++',
    libraries=lapack_info['libraries'],
    library_dirs=lapack_info['library_dirs'],
    sources=['src/map_transmat.pyx', 'src/cycledeterminant.cpp'],
    include_dirs=[np.get_include(), 'include']
)

setup(name='gridmsm',
    version='1.0',
    author='Robert McGibbon',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    py_modules = ['gridmsm'],
    ext_modules=[map_transmat],
    cmdclass={'build_ext': build_ext})
    
      