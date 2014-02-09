"""gridmsm: Markov state models with a grid discretization

More description
"""

DOCLINES = __doc__.split("\n")
from distutils.core import setup

setup(name='gridmsm',
    version='1.0',
    author='Robert McGibbon',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    py_modules = ['gridmsm'])
      