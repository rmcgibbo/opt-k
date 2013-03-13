"""
This library is used to calculate the likelihood of a particular MSM.

The likelihood takes as input the transition matrix as well as the volume
of each state.
"""
from msmbuilder import MSMLib
import numpy as np

def get_log_likelihood(tProb, trajs, volumes, lagtime=1):
    """
    Function to calculate the log likelihood according to:

    .. math:: \log L = \sum_{i,j} C_{ij} \log T_{ij} - 
        \sum_{i,j} c_i \log V_i

    where :math:`C_{ij}` is the number of transitions from state i to state
    j, :math:`T_{ij}` is the probability of transistioning from state i to
    state j, :math:`c_i` is the number of samples of state i, and :math:`V_i` 
    is the volume of state i in phase space.

    Parameters
    ----------
    tProb : scipy.sparse matrix or np.ndarray
        transition probability matrix from msmbuilder
    trajs : np.ndarray
        trajectory of states
    volumes : np.ndarray
        volume for each state in the msm
    lagtime : int
        lag time used to build the msm. Will subsample the trajs accordingly

    Returns
    -------
    log_like : float
        log-likelihood as defined above
    """

    counts = MSMLib.get_count_matrix_from_assignments(trajs, lag_time=lagtime)


    tProb = tProb.tocsr()
    counts.eliminate_zeros()
    state_counts = counts.sum(axis=1)

    ij = counts.nonzero()

    log_like = np.sum(np.array(counts.data) * np.array(np.log(tProb[ij]))) - \
        np.sum(state_counts * np.log(volumes))
    
    return log_like
