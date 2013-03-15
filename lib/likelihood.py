"""
This library is used to calculate the likelihood of a particular MSM.

The likelihood takes as input the transition matrix as well as the volume
of each state.
"""
from msmbuilder import MSMLib
import numpy as np

def get_log_likelihood(tProb, trajs, volumes, lagtime=1, separate=False, 
    kellogg=False):
    """
    Function to calculate the log likelihood according to:

    .. math:: \log L = \sum_{i,j} C_{ij} \log T_{ij} - 
        \sum_{i} c_i \log V_i

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
    separate : bool, optional
        separate transition probability component from the volume component
    kellogg : bool, optional
        if true then use Kellogg's approximation

    Returns
    -------
    log_like : float
        log-likelihood as defined above
    """

    counts = MSMLib.get_count_matrix_from_assignments(trajs, lag_time=lagtime)

    counts = counts.tocsr()
    tProb = tProb.tocsr()

    ij = counts.nonzero()

    state_counts = np.array(counts.sum(axis=1)).flatten()
    nonzero_counts = np.array(counts.data).flatten()
    nonzero_tprobs = np.array(tProb[ij]).flatten()
    volumes = np.array(volumes).flatten()
    # scipy.sparse uses np.matrices which behave differently than arrays
    # so we cast everything to arrays and do the math with them
    
    #log_like = np.sum(np.array(counts.data) * np.array(np.log(tProb[ij]))) - \
    #    np.sum(state_counts * np.log(volumes))
    #log_like = np.sum(nonzero_counts * np.log(nonzero_tprobs)) - \
    #    np.sum(state_counts * np.log(volumes))

    t_like = np.sum(nonzero_counts * np.log(nonzero_tprobs))

    if not kellogg:
        v_like = - np.sum(state_counts * np.log(volumes))
    else:
        v_like = - np.sum(state_counts * np.log(state_counts))

    if separate:
        return t_like, v_like
    
    return t_like + v_like
