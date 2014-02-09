from __future__ import division, print_function
import os
import sys
import numpy as np
import scipy.special
from scipy.misc import logsumexp
from mixtape._reversibility import reversible_transmat_likelihood, logsymsumexp

cimport numpy as np
from libcpp.string cimport string

cdef extern from "cycledeterminant.hpp" namespace "CycleDeterminant":
    cdef cppclass CycleMatrixBuilder:
        CycleMatrixBuilder(const int n_states, string lapack_lite_lib) except +
        double logSqrtDetCycleMatrix(const double* u)
        
    
cdef class MAPTransmat:
    cdef int n_states
    cdef tuple triu_indices
    cdef np.ndarray loop_indices
    cdef np.ndarray edge_indices
    cdef CycleMatrixBuilder * cycleMatrixBuilder

    def __cinit__(self, n_states):
        self.n_states = int(n_states)
        cdef string lapack_lite_lib = os.path.join(np.linalg.__path__[0], 'lapack_lite.so')
        self.cycleMatrixBuilder = new CycleMatrixBuilder(n_states, lapack_lite_lib)

        self.triu_indices = np.triu_indices(n_states)
        # u[loop_indices] are entries of u corresponding to self transitions,
        # u[edge_indices] are entries of u corresponding to non-self transitions
        self.loop_indices = np.where(self.triu_indices[0] == self.triu_indices[1])[0]
        self.edge_indices = np.where(self.triu_indices[0] != self.triu_indices[1])[0]
        

    def fit(self, countsmat, prior):
        """Maximum a-posteriori reversible transition matrix given a set of
        directed transition counts and a symmetric edge-reinforced random walk
        prior.
    
        Parameters
        ----------
        countsmat : np.ndarray
            Dense matrix of directed transition counts
        prior : float
            Strength of the symmetric prior. This is the starting weight on each
            edge in the ERRW.
        """
        symcounts = (countsmat + countsmat.T)[self.triu_indices]
        symcounts[self.loop_indices] -= np.diag(countsmat)
        rowsums = np.sum(countsmat, axis=1)
        logrowsums = np.log(rowsums)
        u0 = np.log(symcounts)

    def reversible_transmat_log_prior(self, np.ndarray[ndim=1, dtype=double] u, double prior):
        """Negative log prior of the reversible transition matrix parameterized
        by the weights in u
    
        Parameters
        ----------
        u : np.array, ndim=1
            The free parameters. These are the log of the upper triangular
            counts matrix entries in symmetric storage.
        prior : float
            Strength of the symmetric prior. This is the starting weight on each
            edge in the ERRW.
        """
        u_edge = u[self.edge_indices]
        u_loop = u[self.loop_indices]

        # \log [ \product_{e \in E\E_loop} x_e^{a_e - 1/2} ]
        term1 = np.sum((prior - 0.5) * u_edge)

        # \log [ \product_{e \in E_loop} x_e^{a_e/2 - 1} ]
        term2 = np.sum((prior/2 - 1) * u_loop)

        u_v = logsymsumexp(u, self.n_states)
        a_v = prior * self.n_states
        # \log [ \product_V x_v^{(a_v + 1)/2} ]
        term3 = np.sum((a_v + 1)/2 * u_v)

        # The determinant of matrix A(x). Equation 12
        term4 = self.cycleMatrixBuilder.logSqrtDetCycleMatrix(&u[0]);
        
        return term4;
