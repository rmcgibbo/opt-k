
from msmbuilder import clustering, io
import numpy as np
from mullermsm.metric import EuclideanMetric
import os


from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-t', dest='data')
#parser.add_argument('-d', dest='dist_cut')
parser.add_argument('-k', dest='num_states', type=int)
parser.add_argument('-o', dest='output_dir')

args = parser.parse_args()

traj = np.load(args.data)

metric = EuclideanMetric()

gen_ids, ass_gen_ids, distances = clustering._kcenters(metric, traj, k=args.num_states, seed=np.random.randint(len(traj)))

ass_gen_ids = np.array([ass_gen_ids])

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

ass_contig = np.ones(ass_gen_ids.shape) * -1
ass_contig = ass_contig.astype(int)

for i, j in enumerate(gen_ids):
    ass_contig[np.where(ass_gen_ids == j)] = i

np.savetxt(os.path.join(args.output_dir, 'gen_ids.dat'), gen_ids)
np.save(os.path.join(args.output_dir, 'gens.npy'), traj[gen_ids])
io.saveh(os.path.join(args.output_dir, 'Assignments.h5'), ass_contig)
io.saveh(os.path.join(args.output_dir, 'Assignments.h5.distances'), np.array([distances]))
print "Saved output to %s" % args.output_dir
