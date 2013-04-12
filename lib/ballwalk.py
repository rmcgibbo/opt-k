import numpy as np

class BallWalk(object):
    """Monte Carlo "Ball-Walk" sampler for computing the relative volume
    of high-dimensional Voronoi cells.
    """
    def __init__(self, metric, generators, space_generators=None, cushion=1.,):
        """Create the BallWalk sampler
        
        Parameters
        ----------
        metric : msmbuilder.metrics.Vectorized
            A vectorized distance metric used to compute distances in this
            space.
        generators : msmbuilder.Trajectory
            The generators that define the centroid of each microstate.
        space_generators : msmbuilder.Trajectory
            This is a set of points in the space that define the envelope
            our data occupies in the space, basically setting the outer
            extent of the volumes. If not supplied, we'll just use the state
            centers (generators) for this purpose.
        cushon : float
            The outer envelope is defined as the set of points within this
            distance from any of the space_generators.
        """

        self.metric = metric
        self.cushion = float(cushion)

        self.generators = generators
        self.prep_generators = self.metric.prepare_trajectory(self.generators)

        if space_generators is None:
            self.space_generators = generators
            self.prep_space_generators = self.prep_generators
            self.space_generators_are_different = False
        else:
            self.space_generators = space_generators
            self.prep_space_generators = self.metric.prepare_trajectory(self.space_generators)
            self.space_generators_are_different = True

        assert self.prep_generators.ndim == 2, 'generators must be 2d'
        self.n_gens, self.n_dims = generators.shape

        # how do we coose this?
        self.ball_radius = 1
        
        # current state of the random walk
        self.counts = np.zeros(len(generators), dtype=np.float)
    
    def _random_ball_point(self):
        # generate a vector dx that is distributed uniformly inside
        # the self.n_dims-dimensional ball, with radius self.ball_radius
        # http://mathworld.wolfram.com/BallPointPicking.html
        normal = np.random.randn(self.n_dims)
        exponential = np.random.exponential()
        dx = self.ball_radius * normal / np.sqrt(exponential + np.sum(np.square(normal)))
                
        return dx
        
    def sample(self, n_steps):
        """Run the BallWalk sampler for a given number of steps
        
        As the sampler runs, it will add counts to self.counts, the normalized
        version of which gives the relative volume of each state.
        """
        
        last_assignment = 0  # the assignment corresponding to the current `state`
        state = self.generators[last_assignment]  # initial point for the walk

        for i in xrange(n_steps):
            # half of the time we just stay put. this ensures the 
            # analytic properties of the sampler I think.
            if np.random.randint(2):

                # new trial point. should we accept it?
                newstate = state + self._random_ball_point()                
                # we need to decide if it's inside the space...
                space_dists = self.metric.one_to_all(newstate.reshape(1,-1), self.prep_space_generators, 0)
                # by checking if its within `cushon` of any of the space_points
                if np.min(space_dists) < self.cushion:
                    # its inside the region. yay! update the state
                    state = newstate
                    
                    if self.space_generators_are_different:
                        # we might not need to recompute these distances if the 
                        # space_generators are the SAME as the regular generators
                        state_dists = self.metric.one_to_all(newstate.reshape(1,-1), self.prep_generators, 0)
                    else:
                        state_dists = space_dists
                    last_assignment = np.argmin(state_dists)

                else:
                    # the proposed move took us out of the space, so we
                    # don't do anything and last_assignments stays
                    # constant
                    pass 
            else:
                # we flipped the wrong coin, so we stay put and last_assignment
                # stays constaitn
                pass

            self.counts[last_assignment] += 1


    def get_state_volumes(self, which_states=None):
        return self.counts / np.sum(self.counts)
        
        
if __name__ == '__main__':
    from msmbuilder import metrics
    class Euclidean(metrics.Vectorized):
        def prepare_trajectory(self, traj):
            return traj

    b = BallWalk(Euclidean(), np.random.randn(10, 3))
    
    import matplotlib.pyplot as pp

    for i in range(10):
        b.sample(10000)
        pp.plot(b.get_state_volumes(), label=str(i))
    
    pp.legend()
    pp.show()