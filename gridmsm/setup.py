"""gridmsm: Markov state models with a grid discretization

More description
"""

DOCLINES = __doc__.split("\n")
from distutils.core import setup, Extension
import numpy as np
from Cython.Distutils import build_ext

map_transmat = Extension('gridmsm._map_transmat',
    language='c++',
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
    
      