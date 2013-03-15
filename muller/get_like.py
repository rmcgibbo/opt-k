import likelihood
import numpy as np
from msmbuilder import io
from scipy.io import mmread


Ls = []
Ls.append(0) # for a one state model the likelihood is 0

folder='KMeans'
for k in range(50, 1050, 50):
    print k
    t = mmread('%s/k%d/Lag30/tProb.mtx' % (folder, k))
    v = np.load('%s/k%d/vols1M.npy' % (folder, k))
    a = io.loadh('%s/k%d/Assignments.h5' % (folder, k))
    a = a['arr_0']
    Ls.append(likelihood.get_log_likelihood(t, a, v, lagtime=30))

np.save('%s.npy' % folder, np.array(Ls))
