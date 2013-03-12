from setuptools import setup

setup(name='opt-k',
      version='0.1',
      description='Likelihood framework for choosing the optimal number of states in a markov state model',
      url='https://github.com/rmcgibbo/opt-k',
      packages=['opt-k'],
      install_requires=['scikit-learn>=0.13.1']
      )
