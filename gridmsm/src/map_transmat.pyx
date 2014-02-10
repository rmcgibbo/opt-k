#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from __future__ import division, print_function
import os
import sys
import numpy as np
import scipy.special
import scipy.optimize
from mixtape._reversibility import reversible_transmat_likelihood, logsymsumexp

cimport cython
cimport numpy as np
from libcpp.string cimport string

#-----------------------------------------------------------------------------
# Globals
#-----------------------------------------------------------------------------

DTYPE = np.float64
ctypedef np.float64_t DTYPE_T

cdef extern from "cycledeterminant.hpp" namespace "CycleDeterminant":
    cdef cppclass CycleMatrixBuilder:
        CycleMatrixBuilder(const int n_states) except +
        double logSqrtDetCycleMatrix(const double* u) except +

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------


cdef class MAPTransmat:
    cdef int n_states
    cdef tuple triu_indices
    cdef np.ndarray loop_indices
    cdef np.ndarray edge_indices
    cdef CycleMatrixBuilder * cycleMatrixBuilder

    def __cinit__(self, n_states):
        self.n_states = int(n_states)
        self.cycleMatrixBuilder = new CycleMatrixBuilder(n_states)

        self.triu_indices = np.triu_indices(n_states)
        # u[loop_indices] are entries of u corresponding to self transitions,
        # u[edge_indices] are entries of u corresponding to non-self transitions
        self.loop_indices = np.where(self.triu_indices[0] == self.triu_indices[1])[0]
        self.edge_indices = np.where(self.triu_indices[0] != self.triu_indices[1])[0]

    def fit(self, np.ndarray[ndim=2, dtype=DTYPE_T] countsmat, DTYPE_T prior, method='cobyla', options={}):
        """Maximum a-posteriori reversible transition matrix given a set of
        directed transition counts and a symmetric edge-reinforced random walk
        prior.
    
        Parameters
        ----------
        countsmat : np.ndarray
            Dense matrix of (directed) transition counts. countsmat[i,j] gives
            the number of observed transitions from state `i` to state `j`.
        prior : float
            Strength of the symmetric prior. This is the starting weight on each
            edge in the ERRW, denoted by :math:`a` in Eq. 13 of Diaconis and 
            Rolles (2006). In that paper, :math:`a` is a vector with a potentially
            different value for every edge. This implementation only supports the
            symmetric prior, where every edge has the same prior strength. Note,
            if prior <= 0, we will just do the MLE estimate.
        method : str
            Optimization method. This string is passed directly to `scipy.optimize.minimize`,
            and can be any of the supported scipy optimization methods.  See
            the scipy documentation for details.
        options : dict
            Options for the optimizer. These are passed directly to `scipy.optimize.minimize`,
            and can be used to control extra printing, the maximum number of
            function evaluations, etc. See the scipy documentation for details.
            
        See Also
        --------
        scipy.optimize.minimize
        """
        cdef np.ndarray[ndim=1, dtype=DTYPE_T] symcounts, symcounts_double_loop
        cdef np.ndarray[ndim=1, dtype=DTYPE_T] rowsums
        cdef np.ndarray[ndim=1, dtype=DTYPE_T] logrowsums
        if countsmat.shape[0] != self.n_states or countsmat.shape[1] != self.n_states:
            raise TypeError('counsmat must be %d by %d' % (self.n_states, self.n_states))

        # symmetrized counts, with the loops counted in both directions
        symcounts_double_loop = (countsmat + countsmat.T)[self.triu_indices]
        
        # symmetrized counts, with the loops counted only once
        symcounts = symcounts_double_loop.copy()
        symcounts[self.loop_indices] -= np.diag(countsmat)
        rowsums = np.sum(countsmat, axis=1)
        logrowsums = np.log(rowsums)
        
        # If values start at precisely negative infinity, the initial
        # value of the objective function is NaN and then the optimization
        # can't go anywhere
        PADDING = 1e-10
        u0 = np.log(symcounts_double_loop + PADDING)

        result = scipy.optimize.minimize(self.objective, u0,
           args=(symcounts, rowsums, logrowsums, prior), method=method, options=options)

        # reconstruct the final counts from the weights. avoid
        # double-counting the diagonal
        reversible_counts = np.zeros((self.n_states, self.n_states))
        reversible_counts[self.triu_indices] = np.exp(result['x'])
        reversible_counts[np.diag_indices(self.n_states)] -= 0.5*np.diag(reversible_counts)
        reversible_counts = reversible_counts + reversible_counts.T
        populations = reversible_counts.sum(axis=0) / reversible_counts.sum()
        transmat = reversible_counts / np.sum(reversible_counts, axis=1)[:, np.newaxis]
        return transmat

    def objective(self, u, symcounts, rowsums, logrowsums, priorstrength):
        """Objective function
        """
        negative_loglikelihood = reversible_transmat_likelihood(u, symcounts, rowsums, logrowsums)
        if priorstrength > 0:
            negative_logprior = self.reversible_transmat_log_prior(u, priorstrength)
        else:
            negative_logprior = 0

        # print('loglikelihood=%f, logprior=%f, value=%f' %
        #     (-negative_loglikelihood, -negative_logprior,
        #     -(negative_loglikelihood + negative_logprior)))
        return negative_loglikelihood + negative_logprior

    def reversible_transmat_log_prior(self, np.ndarray[ndim=1, dtype=DTYPE_T] u, DTYPE_T prior):
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
        cdef np.ndarray[ndim=1, dtype=DTYPE_T] logx = logX_from_u(u, self.n_states)        
        cdef np.ndarray[ndim=1, dtype=DTYPE_T] logx_edge = logx[self.edge_indices]
        cdef np.ndarray[ndim=1, dtype=DTYPE_T] logx_loop = logx[self.loop_indices]

        # \log [ \product_{e \in E\E_loop} x_e^{a_e - 1/2} ]
        term1 = np.sum((prior - 0.5) * logx_edge)

        # \log [ \product_{e \in E_loop} x_e^{a_e/2 - 1} ]
        term2 = np.sum((prior/2 - 1) * logx_loop)

        logx_v = logsymsumexp(logx, self.n_states)
        a_v = prior * self.n_states
        # \log [ \product_V x_v^{(a_v + 1)/2} ]
        term3 = np.sum((a_v + 1)/2 * logx_v)

        # The log sqrt determinant of matrix A(x). Eq 12.
        try:
            term4 = self.cycleMatrixBuilder.logSqrtDetCycleMatrix(&logx[0]);
        except OverflowError:
            term4 = self.n_states * np.log(np.finfo(np.double).max)
        
        # Note that we're MISSING the normalization constant Z from Eq. 14. This
        # is alright since it's not a function of the weights, so it's just
        # an additive constant
        
        #print('logprior',  term1, term2, term3, term4)
        logprior = term1 + term2 - term3 + term4
        return -logprior


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[ndim=1, dtype=DTYPE_T] logX_from_u(np.ndarray[ndim=1, dtype=DTYPE_T] u, int n_states):
    """Calculate the log of the transition probabilities along each edge. These
    are the variables x in Diaconis and  Rolles (2006), and T in the MSMBuilder
    MLE notes.
    
    Parameters
    ----------
    u : np.array, ndim=1
        The free parameters. These are the log of the upper triangular
        counts matrix entries in symmetric storage.

    Returns
    -------
    logX : np.array, ndim=1
        The log of the transition probabilities, also in symmetric storage.
    """
    cdef np.ndarray[ndim=1, dtype=DTYPE_T] q = logsymsumexp(u, n_states)
    cdef np.ndarray[ndim=1, dtype=DTYPE_T] logX = np.zeros(len(u), dtype=DTYPE)
    
    cdef int i, j
    cdef int k = 0
    for i in range(n_states):
        for j in range(i, n_states):
            logX[k] = u[k] - q[i]
            k += 1

    return logX
    