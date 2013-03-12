"""
This library is used to compute the volume of states.
"""

class VolumeEstimator(object):
    """
    This base class can be subclassed by implementations that use an algorithm
    to compute the volume of a state.

    Each implementation should use the distance metric and generators to return
    a list of volumes
    """
    def __init__(self, metric, generators):
    """
    Parameters
    ----------
    metric : msmbuilder.metric.Vectorized instance or subclass instance
        distance metric to use to calculate distances between points. This
        needs to be a Vectorized type because we assume that phase space is
        euclidean
    generators : msmbuilder.Trajectory instance
        each state is defined by the voronoi cell defined by the generators
    """

        if not isinstance(metric, msmbuilder.metrics.Vectorized):
            raise Exception("Metric must be a subclass of msmbuilder.metrics.Vectorized")

        self.metric = metric
        self.generators = generators

        self.prep_generators = self.metric.prepare_trajectory(self.generators)
        self.dimension = self.prep_generators.shape[1]
    
    def get_state_volumes(self, which_states=None):
    """
    Parameters
    ----------
    which_states : np.ndarray, optional
        Calculate volumes only for requested states. If None given, then 
        get volume for each state
    
    Returns
    -------
    volumes : np.ndarray
        volume for each state in same units as metric.prepare_trajectory 
    """

        if which_states is None:
            which_states = np.arange(len(self.generators))

        volumes = np.ones(len(which_states)) * -1
        for i, s in enumerate(which_states):
            print "Working on state %d" % i
            
            volumes[i] = self._get_volume(s)
            
        return volumes
    
    def _get_volume(state):
    """Internal function to get the volume of a particular state.
    """
        raise NotImplementedError
