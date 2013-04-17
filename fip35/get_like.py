import likelihood
import numpy as np
from msmbuilder import io
from scipy.io import mmread
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='msm_dirs', help='MSM directory containing, tProb.mtx, Mapping.dat, and Assignments.Fixed.h5', action='append')
parser.add_argument('-v', dest='vols', help='Volumes to use for each MSM.', action='append')
parser.add_argument('-o', dest='out_fn', default='likelihoods.npy', help='output filename (np.savetxt)')
parser.add_argument('-l', type=int, dest='lagtime', help='lag time (in frames). Right now it only makes sense to compare MSMs built at the same lag time. Though I suppose it is possible to change this in the future.')
args = parser.parse_args()

if len(args.msm_dirs) != len(args.vols):
    print "Need to input as many data directories as volume files."

likelihoods = []
states = []
likelihoods.append(0) # for a one state model the likelihood is 0
states.append(1)

for i in xrange(len(args.msm_dirs)):

    msm_folder = args.msm_dirs[i]
    
    print "Calculating likelihood for MSM in %s" % msm_folder

    t = mmread(os.path.join(msm_folder, 'tProb.mtx'))
    m = np.loadtxt(os.path.join(msm_folder, 'Mapping.dat'))
    v = np.load(args.vols[i])
    a = io.loadh(os.path.join(msm_folder, 'Assignments.Fixed.h5'))
    try:
        a = a['arr_0']
    except:
        a = a['Data']

    v = v[np.where(m != -1)] # remove trimmed states
    likelihoods.append(likelihood.get_log_likelihood(t, a, v, lagtime=args.lagtime, mapping=m))
    states.append(np.where(m != -1)[0].shape[0])

out_ary = np.array([states, likelihoods]).T

np.save(args.out_fn, out_ary)

