#!/usr/bin/env python

from msmbuilder import Trajectory
from msmbuilder import arglib
import numpy as np
from montecarlo import MonteCarlo
from msmbuilder import metrics
import copy

class DihedralWrap(metrics.Vectorized):
    def __init__(*args, **kwargs):
    
        args[0].dih_metric = kwargs.pop('dih_metric')

        metrics.Vectorized.__init__(*args, **kwargs)

    def prepare_trajectory(self, traj):
        ptraj = self.dih_metric.prepare_trajectory(traj)

        N = ptraj.shape[1]

        arcC = np.arccos(ptraj[:,:N/2])
        arcS = np.arcsin(ptraj[:,N/2:])

        arcC[np.where(arcS < 0)] *= -1

        print 'hi'
        return arcC

parser = arglib.ArgumentParser(get_metric=True)

parser.add_argument('generators', help='Generators filename')
parser.add_argument('output', help='output filename to save volumes.')
parser.add_argument('num_points', help='Number of points to sample.')
parser.add_argument('use_dih', action='store_true', default=False)

args, metric = parser.parse_args()

if args.use_dih:
    metric = DihedralWrap(dih_metric=copy.deepcopy(metric))

gens = Trajectory.load_from_lhdf(args.generators)

MC = MonteCarlo(metric, gens, num_points=args.num_points)

volumes = MC.get_state_volumes()

np.save(args.output, volumes)
