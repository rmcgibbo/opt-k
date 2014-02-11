import numpy as np
from glob import glob
import mdtraj as md
import matplotlib.pyplot as pp
from gridmsm import GridMarkovStateModel

N_BINS = np.arange(2, 7)
DATA_GLOB = '/home/rmcgibbo/local/msmbuilder/Tutorial/XTC/RUN*/*.xtc'
TOP = '/home/rmcgibbo/local/msmbuilder/Tutorial/native.pdb'

X = []
for fn in glob(DATA_GLOB):
    t = md.load(fn, top=TOP)
    _, phi = md.compute_phi(t)
    _, psi = md.compute_psi(t)
    X.append(np.hstack((phi, psi)))


akaike = np.empty(len(N_BINS))
schwarz = np.empty(len(N_BINS))
logevidence = np.empty(len(N_BINS))
loglikelihood = np.empty(len(N_BINS))
n_states = np.empty(len(N_BINS))
for i, n in enumerate(N_BINS):
    model = GridMarkovStateModel(n_bins=n, min=-(np.pi+1e-5),
                                 max=np.pi+1e-5, prior=1)
    logevidence[i] = model.logevidence(X)
    model.fit(X)
    akaike[i] = -0.5*model.aic(X)
    schwarz[i] = -0.5*model.bic(X)
    loglikelihood[i] = model.loglikelihood(X)
    n_states[i] = model.n_states

    model.sample(X, 10)


pp.plot(n_states, logevidence, label='Log Evidence', lw=2)
pp.plot(n_states, loglikelihood, label='MLE Log Likelihood', lw=2)
pp.plot(n_states, schwarz, label='MLE Schwarz Criterion', lw=2)
pp.plot(n_states, akaike, label='MLE Akaike Criterion', lw=2)
pp.xlabel('Number of states')
pp.legend(loc=4)
pp.show()
