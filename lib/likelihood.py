"""
This library is used to calculate the likelihood of a particular MSM.

The likelihood takes as input the transition matrix as well as the volume
of each state.
"""
from msmbuilder import MSMLib
import numpy as np

def get_log_likelihood(tProb, trajs, volumes, lagtime=1, separate=False, 
    kellogg=False, mapping=None):
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

    if not mapping is None:
        MSMLib.apply_mapping_to_assignments(trajs, mapping)

    counts = MSMLib.get_count_matrix_from_assignments(trajs, lag_time=lagtime,
        n_states=tProb.shape[0])

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


def get_kmeans_log_likelihood(tprob, trajs, state_dists, num_dims, 
    lagtime=1, train_dists=None, train_ass=None, mapping=None, add_weights=False):
    """
    get the log likelihood of an MSM given that each state emits a conformation
    according to a user-input probability model. 

    Parameters
    ----------
    tprob : scipy.sparse matrix
        transition probability matrix for your MSM
    trajs : np.ndarray
        two-dim. array corresponding to the assignments of each trajectory
    state_dists : np.ndarray
        distance of a conformation to its generator
    num_dims : int
        number of dimensions in your vector space
    lagtime : int, optional
        lag time of the MSM (default: 1)
    train_dists : np.ndarray, optional
        the distance to the cluster center for the training data. If this
        is not None, then we assume you are doing cross-validation, 
        which means the state variance is estimated from these distances
        and the state_dists are used in the likelihood separately
    mapping : mapping to map counts to a transition matrix
    add_weights : bool, optional
        add the per-state weights corresponding to R_n / R from the X-means 
        paper (default: False)

    Returns
    -------
    log_likelihood : float
        log likelihood of the model

    Notes
    -----
    These results come from Pelleg, D. and Moore A. X-means: Extending K-means
        with Efficient Estimation of the Number of Clusters. 
        (http://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf)
    """

    cross = False
    if not train_dists is None and not train_ass is None:
        cross = True

    if not mapping is None:
        MSMLib.apply_mapping_to_assignments(trajs, mapping)

    counts = MSMLib.get_count_matrix_from_assignments(trajs, lag_time=lagtime,
        n_states=tprob.shape[0])

    counts = counts.tocsr()
    tprob = tprob.tocsr()

    ij = counts.nonzero()

    state_counts = np.array(counts.sum(axis=1)).flatten()
    nonzero_counts = np.array(counts.data).flatten()
    nonzero_tprobs = np.array(tprob[ij]).flatten()
    # scipy.sparse uses np.matrices which behave differently than arrays
    # so we cast everything to arrays and do the math with them
   
    #log_like = np.sum(np.array(counts.data) * np.array(np.log(tProb[ij]))) - \
    #    np.sum(state_counts * np.log(volumes))
    #log_like = np.sum(nonzero_counts * np.log(nonzero_tprobs)) - \
    #    np.sum(state_counts * np.log(volumes))

    t_like = np.sum(nonzero_counts * np.log(nonzero_tprobs))

    #R = np.sum(state_counts) # total number of points
    R = np.where(trajs != -1)[0].shape[0]
    if cross:
        orig_R = np.where(train_ass != -1)[0].shape[0]
        state_var = np.square(train_dists[np.where(train_ass != -1)]).sum() / (orig_R - tprob.shape[0])
    else:
        state_var = np.square(state_dists[np.where(trajs != -1)]).sum() / (R - tprob.shape[0])

    print 'state variance: %.4e' % state_var
    s_like = - num_dims / 2. * R * np.log(2 * np.pi * state_var) 
    if cross:
        s_like -= np.square(state_dists[np.where(trajs != -1)]).sum() / (2. * state_var)
    else:
        s_like -= (R - tprob.shape[0]) / 2.
        
    print 't-like: %.4e s-like: %.4e tot: %.4e' % (t_like, s_like, t_like + s_like)
    return t_like + s_like

