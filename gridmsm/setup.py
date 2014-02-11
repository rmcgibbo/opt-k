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

try:
    import scipy
    if scipy.version.full_version < '0.11':
        raise ImportError()
except ImportError:
    print('scipy version 0.11 or better is required.')


def get_lapack():
    from collections import defaultdict
    lapack_info = defaultdict(lambda: [])
    lapack_info.update(system_info.get_info('lapack'))
    if len(lapack_info) == 0:
        try:
            from scipy.linalg import _flapack
            lapack_info['extra_link_args'] = [_flapack.__file__]
            return lapack_info
        except ImportError:
            pass
        print('LAPACK libraries could not be located.', file=sys.stderr)
        sys.exit(1)
    return lapack_info

lapack_info = get_lapack()
map_transmat = Extension('gridmsm._bayes_transmat',
    language='c++',
    libraries=lapack_info['libraries'],
    library_dirs=lapack_info['library_dirs'],
    sources=['src/bayes_transmat.pyx', 'src/cycledeterminant.cpp'],
    include_dirs=[np.get_include(), 'include'],
    extra_link_args=lapack_info['extra_link_args'],
)

setup(name='gridmsm',
    version='1.0',
    author='Robert McGibbon',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    packages = ['gridmsm'],
    ext_modules=[map_transmat],
    cmdclass={'build_ext': build_ext})
    
      
