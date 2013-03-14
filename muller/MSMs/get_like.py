import likelihood
import numpy as np
from msmbuilder import io
from scipy.io import mmread


Ls = []
Ls.append(0) # for a one state model the likelihood is 0

for k in range(50,1050,50):
    print k
    t = mmread('k%d/Lag30/tProb.mtx' % k)
    v = np.load('k%d/vols.npy' % k)
    a = io.loadh('k%d/Assignments.h5' % k)
    a = a['arr_0']
    Ls.append(likelihood.get_log_likelihood(t, a, v, lagtime=30))

np.save('likes1000_s50.npy', np.array(Ls))
