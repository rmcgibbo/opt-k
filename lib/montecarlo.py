
from msmbuilder import metrics
import numpy as np
from volumes import VolumeEstimator

class MonteCarlo(VolumeEstimator):
    """
    This class uses a brute force Monte Carlo algorithm to compute the volume
    of a state by randomly drawing samples from R^d and assigning them to a 
    state.
    """
    
    def __init__(self, metric, generators, num_points=1E6, cushion=1.):
        """
        Compute the volume of states by using a simple monte carlo algorithm.

        Parameters
        ----------
        metric : msmbuilder.metrics.Vectorized
            metric to compute distances. Must be vectorized because we assume
            that phase space is euclidean
        generators : msmbuilder.Trajectory
            generators that define the voronoi cells
        num_points : int, optional
            number of points to simulate to calculate the volumes
        cushion : float, optional
            cushion to use to define the bounding box for the data. The box
            will be defined based on a hyberrectangle with corners in each 
            dimension equal to 
                [min_d(generators) - cushion, max_d(generators) + cushion]
        """
        super(MonteCarlo, self).__init__(metric, generators)

        self.num_points = int(num_points)

        self.mins = self.prep_generators.min(axis=0) - cushion
        self.maxs = self.prep_generators.max(axis=0) + cushion
        self.lengths = self.maxs - self.mins
        
        self.total_volume = np.prod(self.lengths)

    def get_state_volumes(self, which_states=None):
        """
        Compute the volume of each state.
        
        Parameters
        ----------
        which_states : Placeholder. This is not used.

        Returns
        -------
        volumes : np.ndarray
            array of volumes for each state

        """
        # This algorithm is more efficient if we do all states at the same time
        # so we will not use the same framework as the base class VolumeEstimator
    
        state_counts = np.zeros(len(self.prep_generators))

        for i in xrange(self.num_points):
        # It might be more efficient to write this vectorially, but I don't 
        # want to run out of memory
            random_sample = np.random.random(self.dimension)
            random_sample = random_sample * self.lengths + self.mins
            random_sample = random_sample.reshape((1, -1))
            
            state_dists = self.metric.one_to_all(random_sample, self.prep_generators, 0)

            state_counts[np.argmin(state_dists)] += 1

        points_per_state = state_counts.astype(float) / self.num_points

        volumes = points_per_state * self.total_volume

        return volumes
        
    
    def _get_volume(self, state):
        raise Exception("""
            This algorithm is more efficient if we do all states at the same time
            so we will not use the same framework as the base class VolumeEstimator
            """)
