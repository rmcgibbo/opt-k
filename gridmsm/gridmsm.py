import numpy as np
import scipy.special
from msmbuilder import MSMLib as msmlib
__all__ = ['GridMarkovStateModel']


def _log_transmat_evidence(countsmat, starting_state, prior):
    """Compute the log evidence of a Markov model given directed transition counts
    
    This function implements Eq. 51 from [1], using a symmetric prior (all
    :math:`\alpha_e`) are equal
    
    Parameters
    ----------
    countsmat : np.ndarray
        Dense matrix of directed transition counts
    starting_state : int
        The initial state of the chain
    prior : float
        Strength of the symmetric prior. This is the starting weight on each
        edge in the ERRW.
    
    References
    ----------
    .. [1] Diaconis, Persi, and Silke WW Rolles. "Bayesian analysis for
           reversible Markov chains." The Annals of Statistics 34.3 (2006):
           1270-1292.
    """
    
    def logpochhammer(a, n):
        "Natural logarithm of the Pochhammer symbol (a)_n"
        return scipy.special.gammaln(a + n) - scipy.special.gammaln(a)

    countsmat = np.asarray(countsmat)
    n_states = countsmat.shape[0]
    k_v = np.sum(countsmat, axis=1)
    k_e = countsmat + countsmat.T

    numerator1 = 0
    for e in zip(*np.triu_indices(n_states, 1)):
        # e in E\E_loop
        numerator1 += logpochhammer(prior, k_e[e])

    numerator2 = 0
    for k_e in np.diag(k_e):
        # e in E_loop
        numerator2 += (k_e/2.0) * np.log(2.0) + \
                      logpochhammer(prior / 2.0, k_e/2.0)

    denominator = 0
    for v in range(n_states):
        offset = int(v != starting_state)
        denominator += (k_v[v]) * np.log(2.0) + logpochhammer((
            prior * n_states + offset) / 2.0, k_v[v])

    return (numerator1 + numerator2 - denominator)


class GridMarkovStateModel(object):
    """Markov State Model using a grid-based discretization

    Parameters
    ----------
    n_bins : int
        Number of bins along per degree of freedom. The total of states for
        the Markov model will be :math:`n_{bins}^{n_{dofs}}`. Where
        :math:`n_{dofs}` is the number of degrees of freedom.
    min : float
        The bins along each degree of freedom run from min to max.
    max : float
        The bins along each degree of freedom run from min to max.
    prior : float
        Strength of the symmetric edge reinforced random walk conjugate prior
        on the (reversible) transition matrix. The interpretation of the prior
        strength is very analogous to a symmetric Dirichlet prior -- the prior
        acts like a "psuedocount" on each edge of the graph.

    Attributes
    ----------
    """
    def __init__(self, n_bins, min, max, transmat_prior=1.0):
        self.n_bins = n_bins
        self.min = min
        self.max = max
        self.transmat_prior = transmat_prior

        self.n_states = None
        self.n_features = None
        self.transmat_ = None
        self.grid = np.linspace(self.min, self.max, self.n_bins + 1)

    def _n_parameters(self):
        return (self.n_states * (self.n_states + 1) / 2) - 1

    def _initialize_sequences(self, X):
        if not isinstance(X, list):
            raise ValueError('X must be a list')

        if self.n_features is None:
            if not isinstance(X, np.ndarray) or X[0].ndim != 2:
                raise ValueError("Each element must be a 2d array")
            self.n_features = X[0].shape[1]

        self.n_states = self.n_bins**self.n_features
        for xx in X:
            if not isinstance(xx, np.ndarray) or xx.ndim != 2 or xx.shape[1] != self.n_features:
                raise ValueError("Each element must be a 2D array "
                                 "with of shape N by %d" % self.n_features)

    def _discretize(self, datasequence):
        """Discretize a real-valued data sequence, projecting each observation
        into an integer-indexed state

        Parameters
        ----------
        datasequence : np.ndarray, shape=(n_observations, n_features)
            A single trajectory: a timeseries of n_dim dimensional observations
            of the system moving through phase space

        Returns
        -------
        labels : np.ndarray
             The labels for each observation
        """
        assert self.n_features == datasequence.shape[1]
        labels = np.zeros(len(datasequence))
        for i in range(self.n_features):
            labels += i*(np.digitize(datasequence[:, i], self.grid) - 1)
        return labels

    def fit(self, X, method='mle'):
        """Fit the maximum a posteriori (MAP) or maximum likelihood (MLE)
        transition matrix to a collection of sequences

        This method sets the `transmat_` attribute on the model

        Parameters
        ----------
        X : list of np.ndarrays
            Each array is a 1D trajectory through phase space. This should
            already be subsampled to your discretization interval / lag time
            of choice.
        method : {'mle', 'map'}
            Use the maximum likelihood or maximum a posterori method to fit the
            transition matrix. If 'mle', the `transmat_prior` field is ignored.

        Returns
        -------
        self
        """
        if not method in ['mle', 'map']:
            raise ValueError('method must be one of ["mle", "map"]')
        self._initialize_sequences(X)

        countsmat = None
        for trajectory in X:
            labels = self._discretize(trajectory)
            thesecounts = msmlib.get_counts_from_traj(labels, n_states=self.n_states)
            if countsmat is None:
                countsmat = thesecounts
            else:
                countsmat = countsmat + thesecounts

        if method == 'mle':
            rc = msmlib.mle_reversible_count_matrix(
                scipy.sparse.csr_matrix(countsmat.todense() + 1e-20))
            self.transmat_ = msmlib.estimate_transition_matrix(rc)
        else:
            raise NotImplementedError('Still working on MAP')

    def loglikelihood(self, X, terms='all'):
        """Log-likelihood of the current model fit and the proposed data

        Parameters
        ----------
        X : list of np.ndarrays
            Each array is a 1D trajectory through phase space. This should
            already be subsampled to your discretization interval / lag time
            of choice.
        terms : str
            One of ['transmat', 'emission', 'all']. default='all'
        """
        if not terms in ['transmat', 'emission', 'all']:
            raise ValueError("terms must be one ['transmat', 'emission', 'all']")
        logtransmat = np.log(self.transmat_)
        transition_log_likelihood = 0
        emission_log_likelihood = 0

        for trajectory in X:
            labels = self._discretize(trajectory)
            transition_log_likelihood += np.multiply(msmlib.get_counts_from_traj(
                labels, n_states=self.n_states).todense(), logtransmat).sum()

            if np.any(np.logical_or(trajectory < self.min, trajectory > self.max)):
                emission_log_likelihood += -np.inf
            emission_log_likelihood += -self.n_features * \
                np.log(self.grid[1] - self.grid[0]) * len(labels)

        if terms == 'transmat':
            return transition_log_likelihood
        elif terms == 'emission':
            return emission_log_likelihood
        return transition_log_likelihood + emission_log_likelihood

    def bic(self, X):
        """Bayesian information criterion for the current model fit and the proposed data

        Parameters
        ----------
        X : list of np.ndarrays
            Each array is a 1D trajectory through phase space. This should
            already be subsampled to your discretization interval / lag time
            of choice.

        Returns
        -------
        bic: float (the lower the better)
        """
        return -2 * self.loglikelihood(X) + self._n_parameters * np.log(sum(len(xx) for xx in X))

    def aic(self, X):
        """Akaike information criterion for the current model fit and the proposed data

        Parameters
        ----------
        X : list of np.ndarrays
            Each array is a 1D trajectory through phase space. This should
            already be subsampled to your discretization interval / lag time
            of choice.

        Returns
        -------
        aic: float (the lower the better)
        """
        return -2 * self.loglikelihood(X) + 2 * self._n_parameters

    def logevidence(self, X, terms='all'):
        """Bayesian evidence for the model. This is the probability of the data
        given the model, integrated over all of the possible transition matrices
        (with respect to their prior).

        .. math ::

            P(X | M) = \integral_{\theta} P(X | \theta) P(\theta | M) d\theta

        Parameters
        ----------
        X : list of np.ndarrays
            Each array is a 1D trajectory through phase space. This should
            already be subsampled to your discretization interval / lag time
            of choice.
        terms : str
            One of ['transmat', 'emission', 'all']. default='all'

        Returns
        -------
        logevidence : float (the higher the better)
        """
        if not terms in ['transmat', 'emission', 'all']:
            raise ValueError("terms must be one ['transmat', 'emission', 'all']")
        self._initialize_sequences(X)

        logevidence_transmat = 0
        emission_log_likelihood = 0
        for trajectory in X:
            labels = self._discretize(trajectory)
            starting_state = labels[0]
            countsmat = msmlib.get_counts_from_traj(labels, n_states=self.n_states).todense()

            logevidence_transmat += _log_transmat_evidence(
                countsmat, starting_state, self.transmat_prior)

            # this factors out when we integrate over model parameters
            # since the we're NOT integrating over the state-space parameters
            # as they're fixed given the number of states
            if np.any(np.logical_or(trajectory < self.min, trajectory > self.max)):
                emission_log_likelihood += -np.inf
            emission_log_likelihood += -self.n_features * \
                np.log(self.grid[1] - self.grid[0]) * len(labels)

        if terms == 'transmat':
            return logevidence_transmat
        elif terms == 'emission':
            return emission_log_likelihood
        return logevidence_transmat + emission_log_likelihood



