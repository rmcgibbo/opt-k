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

class Euclidean(metrics.Vectorized):
    def prepare_trajectory(self, traj):
        return traj

parser = arglib.ArgumentParser(get_metric=True)

parser.add_argument('generators', help='Generators filename')
parser.add_argument('space_gens', help='Generators to use to define the whole space')
parser.add_argument('cushion', type=float, help='Cushion to require points be within a generator to include it as a random sample.')
parser.add_argument('output', help='output filename to save volumes.')
parser.add_argument('num_points', help='Number of points to sample.')
parser.add_argument('use_dih', action='store_true', default=False)
parser.add_argument('use_euc', action='store_true', default=False)
#parser.add_argument('space_metric', default=None, help='metric to use to define which points are in the whole space. None or pickled distance metric')

args, metric = parser.parse_args()

if args.use_euc:
    metric = Euclidean()
elif args.use_dih:
    metric = DihedralWrap(dih_metric=copy.deepcopy(metric))

try: 
    gens = np.load(args.generators)
except: 
    gens = Trajectory.load_from_lhdf(args.generators)

try:
    space_gens = np.load(args.space_gens)
except:
    space_gens = Trajectory.load_from_lhdf(args.space_gens)
#
#if not args.space_metric is None:
#    space_metric = pickle.unpickle(args.space_metric)
#else:
#    space_metric = None

MC = MonteCarlo(metric, gens, space_gens, num_points=args.num_points, 
        cushion=args.cushion)
#, space_metric=space_metric)

volumes = MC.get_state_volumes()

np.save(args.output, volumes)
