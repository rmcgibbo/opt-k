
from msmbuilder import metrics
import numpy as np
from volumes import VolumeEstimator

class MonteCarlo(VolumeEstimator):
    """
    This class uses a brute force Monte Carlo algorithm to compute the volume
    of a state by randomly drawing samples from R^d and assigning them to a 
    state.
    """
    
    def __init__(self, metric, generators, space_generators=None,
         num_points=1E6, cushion=1.):
        """
        Compute the volume of states by using a simple monte carlo algorithm.

        Parameters
        ----------
        metric : msmbuilder.metrics.Vectorized
            metric to compute distances. Must be vectorized because we assume
            that phase space is euclidean
        generators : msmbuilder.Trajectory
            generators that define the voronoi cells
        space_generators : msmbuilder.Trajectory, optional
            conformations that define the full space. If None, then generators
            are used
        num_points : int, optional
            number of points to simulate to calculate the volumes
        cushion : float, optional
            cushion to use to define full space of the data. A point is in
            this space if it is within <cushion> of any of the space_generators
        """
        super(MonteCarlo, self).__init__(metric, generators)

        self.num_points = int(num_points)
        self.cushion = float(cushion)

        self.space_generators = space_metrics
        self.prep_space_generators = self.metric.prepare_trajectory(self.space_generators)

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

        while state_counts.sum() < num_points:
        # It might be more efficient to write this vectorially, but I don't 
        # want to run out of memory
            random_sample = np.random.random(self.dimension)
            random_sample = random_sample * self.lengths + self.mins
            random_sample = random_sample.reshape((1, -1))
            
            space_dists = self.metric.one_to_all(random_sample, self.prep_space_generators, 0)

            if np.min(space_dists) > self.cushion:
                continue
            
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
