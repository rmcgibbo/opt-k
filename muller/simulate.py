import random
from mullermsm import muller
mullerforce = muller.muller_force()
import scipy.linalg
from matplotlib.pyplot import *
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--dt', dest='dt', type=float, default=0.1)
parser.add_argument('-n', dest='num_frames', type=int, default=100000)
parser.add_argument('-o', dest='output', default='pos.npy')

args = parser.parse_args()

kT = 15.0
dt = args.dt
mGamma = 1000.0
traj_length = args.num_frames 
initial_x = [random.uniform(-1.5, 1.2), random.uniform(-0.2, 2)]
positions = muller.propagate(traj_length, initial_x, kT, dt, mGamma, mullerforce)

np.save(args.output, positions)
