import os
import IPython as ip
import numpy as np
from mdtraj import io
import scipy.sparse
import scipy.special
from msmbuilder import MSMLib as msmlib
from msmbuilder import  msm_analysis
from matplotlib import mlab
import matplotlib.pyplot as pp
colors = ['#348ABD', 'g', '#7A68A6', '#A60628', '#467821', '#CF4457', '#188487', '#E24A33']

DIFFUSION_CONST = 1
DT = 1e-3
SQRT_2D_DT = np.sqrt(2 * DIFFUSION_CONST * DT)
grad_potential = lambda x : -2 * np.sin(2 * x)


def reflect_pbc(x, min, max):
    if x > max:
        return 2*max - x
    if x < min:
        return 2*min - x
    return x


def propagate(n_steps):
    import time; start = time.time()
    n_steps = int(n_steps)

    max = np.pi
    min = -np.pi
    rand = np.random.randn(n_steps)
    x = np.zeros(n_steps+1)
    for i in range(n_steps):
        x_i_plus_1 = x[i] -DT * grad_potential(x[i]) + SQRT_2D_DT*rand[i]
        x[i+1] = reflect_pbc(x_i_plus_1, -np.pi, np.pi)
        # reflecting PBCs

    print '%d steps/s' % (n_steps / (time.time() - start))
    return x


def exact_solution(n_grid, plot_eigfunctions=False):
    ONE_OVER_SQRT_2PI = 1.0 / (np.sqrt(2*np.pi))
    normalpdf = lambda x : ONE_OVER_SQRT_2PI * np.exp(-0.5 * (x*x))

    grid = np.linspace(-np.pi, np.pi, n_grid)
    width = grid[1]-grid[0]
    transmat = np.zeros((n_grid, n_grid))
    for i, x_i in enumerate(grid):
        for offset in range(-(n_grid-1), n_grid):
            x_j = x_i + (offset * width)
            j = reflect_pbc(i+offset, 0, n_grid-1)

            # What is the probability of going from x_i to x_j in one step?
            diff = (x_j - x_i + DT * grad_potential(x_i)) / SQRT_2D_DT
            transmat[i, j] += normalpdf(diff)
        #print transmat[i, :]
        transmat[i, :] =  transmat[i, :] / np.sum(transmat[i, :])


    eigvalues, eigvectors = np.linalg.eig(transmat)
    eigsort = np.argsort(np.real(eigvalues))[::-1]
    eigvectors = eigvectors[:, eigsort]
    eigvalues = eigvalues[eigsort]
    
    if plot_eigfunctions:
        pp.title('Double well transfer operator 1st eigenfunction')
        eigvectors[:, 1] /= np.max(eigvectors[:, 1])
        
        pp.plot(grid, eigvectors[:, 1], label='Exact Soln.', lw=3)
        pp.plot([-np.pi, 0, 0, np.pi], [0.85, 0.85, -0.85, -0.85], label='2-state MSM', lw=2)

        xx = np.linspace(-np.pi, np.pi, 10+1)
        vv = np.zeros(10+1)
        for i in range(10):
            center = (xx[i] + xx[i+1]) / 2
            vv[i] = eigvectors[np.argmin((grid-center)**2), 1]
        for i in range(1, 10):
            pp.plot([xx[i], xx[i]], [vv[i-1], vv[i]], c='k', lw=2)
        for i in range(10):
            pp.plot([xx[i], xx[i+1]], [vv[i], vv[i]], c='k', lw=2)
        pp.plot([0,0], [0,0], c='k', lw=2, label='10-state MSM')
        #pp.yticks([])
        pp.legend()
        pp.ylim(-1.2, 1.2)
        pp.xlim(-np.pi, np.pi)
        print 'saving eigenfunctions'
        pp.savefig('doublewell-eigenfunctions-msm.png')
        exit(1)


    eigs = np.sort(np.real(np.linalg.eigvals(transmat)))
    timescales =  -1 / np.log(eigs)
    return timescales


class MarkovStateModel1D(object):
    def __init__(self, n_states, discretization='grid', min_x=-np.pi, max_x=np.pi, transmat_prior=1.0):
        self.n_states = n_states
        self.discretization = discretization
        self.min_x = min_x
        self.max_x = max_x
        self.transmat_prior = transmat_prior


        self.transmat_ = None
        self.populations_ = None
        self.timescales_ = None
        self.grid_ = None
        self.starting_states_ = None

        self._n_parameters = (n_states * (n_states + 1) / 2) - 1

    def fit(self, X):
        """

        Parameters
        ----------
        X : list of np.ndarrays
            Each array is a 1D trajectory through phase space. This should
            already be subsampled to your discretization interval / lag time
            of choice.
        """

        if self.discretization != 'grid':
            raise ValueError()

        self.grid_ = np.linspace(self.min_x, self.max_x, self.n_states + 1)
        
        self.starting_states_ = []
        self.countsmats_ = []
        for xx in X:
            labels = np.digitize(xx, self.grid_) - 1
            self.starting_states_.append(labels[0])
            countsmat = np.asarray(msmlib.get_counts_from_traj(labels, n_states=self.n_states).todense())
            self.countsmats_.append(countsmat)

        rc = msmlib.mle_reversible_count_matrix(scipy.sparse.csr_matrix(self.transmat_prior + np.sum(self.countsmats_, axis=0)))
        self.transmat_ = msmlib.estimate_transition_matrix(rc)
        self.timescales_ = -1 / msm_analysis.get_eigenvectors(self.transmat_, self.n_states-2)[0][1:]
        self.transmat_ = self.transmat_.todense()
    
    def loglikelihood(self, X, terms='all'):
        """

        Parameters
        ----------
        X : list of np.ndarrays
            Each array is a 1D trajectory through phase space. This should
            already be subsampled to your discretization interval / lag time
            of choice.
        terms : str
            One of ['transmat', 'emission', 'all']. default='all'
        """
        assert terms in ['transmat', 'emission', 'all']
        logtransmat = np.log(self.transmat_)
        transition_log_likelihood = 0
        emission_log_likelihood = 0

        for xx in X:
            labels = np.digitize(xx, self.grid_) - 1
            transition_log_likelihood += np.multiply(msmlib.get_counts_from_traj(
                labels, n_states=self.n_states).todense(), logtransmat).sum()
            
            #if np.any(np.logical_or(xx < self.min_x, xx > self.max_x)):
            #    emission_log_likelihood += -np.inf
            emission_log_likelihood += -np.log(self.grid_[1] - self.grid_[0]) * len(labels)

        if terms == 'transmat':
            return transition_log_likelihood
        elif terms ==  'emission':
            return emission_log_likelihood
        return transition_log_likelihood + emission_log_likelihood

    def bic(self, X):
        """
        Bayesian information criterion for the current model fit
        and the proposed data

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
        """Akaike information criterion for the current model fit
        and the proposed data

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

    def logevidence(self, terms='all'):

        assert terms in ['transmat', 'emission', 'all']
        def logpochhammer(a, n):
            "Natural logarithm of the Pochhammer symbol (a)_n"
            return scipy.special.gammaln(a + n) - scipy.special.gammaln(a)

        logevidence_transmat = 0
        for countsmat, starting_state in zip(self.countsmats_, self.starting_states_):
            k_v = np.sum(countsmat, axis=1)
            k_e = countsmat + countsmat.T

            numerator1 = 0
            for e in zip(*np.triu_indices(self.n_states, 1)):
                # e in E\E_loop
                numerator1 += logpochhammer(self.transmat_prior, k_e[e])
            
            numerator2 = 0
            for k_e in np.diag(k_e):
                # e in E_loop
                numerator2 += (k_e/2.0) * np.log(2.0) + \
                              logpochhammer(self.transmat_prior / 2.0, k_e/2.0)
    
            denominator = 0
            for v in range(self.n_states):
                offset = int(v != starting_state)
                denominator += (k_v[v]) * np.log(2.0) + logpochhammer((
                    self.transmat_prior * self.n_states + offset) / 2.0, k_v[v])

            logevidence_transmat += (numerator1 + numerator2 - denominator)

        
        # this factors out when we integrate over model parameters
        # since the we're NOT integrating over the state-space parameters
        # as they're fixed given the number of states
        if terms in ['emission', 'all']:
            emission_log_likelihood = - np.log(self.grid_[1]-self.grid_[0]) * \
                                      (np.sum(self.countsmats_) + len(self.countsmats_))


        if terms == 'transmat':
            return logevidence_transmat
        elif terms ==  'emission':
            return emission_log_likelihood
        return logevidence_transmat + emission_log_likelihood




def simulate():
    trajectories = []
    for i in range(10):
        fn = 'trajectory-%d.h5' % i
        if os.path.exists(fn):
            print 'loading %s' % fn
            trajectories.append(io.loadh(fn)['arr_0'])
        else:
            x = propagate(5e5)
            io.saveh(fn, x)
            print 'saving %s' % fn
            trajectories.append(x)
    return trajectories

def main():
    trajectories = simulate()
    model = MarkovStateModel1D(n_states=2)
    t = [trajectories[i][::100] for i in range(2)]
    model.fit(t)
    print 'loglikelihood', model.loglikelihood(t)
    print 'evidence', model.logevidence()
    print 'bic', -model.bic(t)/2
    print 'aic', -model.aic(t)/2


def test_logevidence():
    model = MarkovStateModel1D(n_states=4)
    model.countsmats_ = [np.array([
        [91, 160, 261, 108],
        [213, 351, 161, 249],
        [251, 224, 388, 201],
        [66, 239, 254, 152]])]
    model.starting_states_ = [3]

    result = model.logevidence('transmat')
    expected = np.log(2.166939224648291) - 1961 * np.log(10)
    np.testing.assert_almost_equal(expected, result)

if __name__ == '__main__':
    pp.ion()
    main()
    test_logevidence()
