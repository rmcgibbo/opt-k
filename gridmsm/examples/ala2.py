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
    break

akaike = np.empty(len(N_BINS))
schwarz = np.empty(len(N_BINS))
logevidence = np.empty(len(N_BINS))
loglikelihood = np.empty(len(N_BINS))
n_states = np.empty(len(N_BINS))
timescale_1 = np.empty(len(N_BINS))
timescale_2 = np.empty(len(N_BINS))

for i, n in enumerate(N_BINS):
    model = GridMarkovStateModel(n_bins=n, min=-(np.pi+1e-5),
                                 max=np.pi+1e-5, prior=1)
    model.fit(X)
    n_states[i] = model.n_states
    logevidence[i] = model.logevidence(X)
    loglikelihood[i] = model.loglikelihood(X)
    akaike[i] = -0.5*model.aic(X)
    schwarz[i] = -0.5*model.bic(X)
    ts = -1/np.log(np.sort(np.linalg.eigvals(model.transmat_)))
    timescale_1[i] = ts[-2]
    timescale_2[i] = ts[-3]


pp.figure()
pp.subplot(211)
pp.plot(n_states, logevidence, 'b-o', label='Log Evidence', lw=2)
pp.plot([n_states[np.argmax(logevidence)]], [np.max(logevidence)], 'b^', ms=10, lw=2)

pp.plot(n_states, loglikelihood, 'g-o', label='MLE Log Likelihood', lw=2)
pp.plot([n_states[np.argmax(loglikelihood)]], [np.max(loglikelihood)], 'g^', ms=10, lw=2)

pp.plot(n_states, schwarz, 'c-o', label='MLE Schwarz Criterion', lw=2)
pp.plot([n_states[np.argmax(schwarz)]], [np.max(schwarz)], 'c^', ms=10, lw=2)

pp.plot(n_states, akaike, 'm-o', label='MLE Akaike Criterion', lw=2)
pp.plot([n_states[np.argmax(akaike)]], [np.max(akaike)], 'm^', ms=10, lw=2)

pp.xlabel('Number of states')
pp.legend(loc=4)

pp.title('ALA2 Phi/Psi-Grid MSM: Model Selection')
pp.subplot(212)
pp.plot(n_states, timescale_1, '--', lw=2, c='k', label='MLE ITS 1')
pp.plot(n_states, timescale_2, '--', lw=2, c='k', label='MLE ITS 2')
pp.ylabel('Timescale')
pp.xlabel('Number of States')
pp.legend(loc=4)


pp.show()
