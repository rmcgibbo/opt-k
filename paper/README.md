This file was automatically generated from latex using pandoc

Likelihood framework for choosing the optimal number of states in a
Markov State Model

Overview
========

A significant challenge in the automated construction of markov state
models is choosing the number of microstates. The number of microstates,
![equation](http://latex.codecogs.com/gif.latex?k) governs a bias-variance tradeoff – when ![equation](http://latex.codecogs.com/gif.latex?k) is large, the models
have a lot of parameters enabling it to fit the data with low bias, but
the simultaneous estimation of these parameters from finite data can
lead to high variance. On the other hand, when ![equation](http://latex.codecogs.com/gif.latex?k) is small, we have
fewer parameters and better statistics, but at the expense of model
flexibility.

Currently, most MSM construction pipelines involve manually selecting
the number of states. This is done primarily by human intuition, and is
a significant bottleneck.

Here, we introduce a likelihood framework enabling the automatic
selection of the number of states.

Within the discrete state space, the likelihood of an ensemble of
trajectories given the MSM is straightforward. We simply take the
product of the state to state transition probabilities along the path.
These state to state transition probabilities are our main central
parameters, the entries of the transition matrix. However, as we vary
the number of states, it is not permissible to simply compare these
likelihoods to select an optimal number of states. Doing so, we would
always chose the trivial model, with only one state, as the transition
probability in that model would be ![equation](http://latex.codecogs.com/gif.latex?1) and thus the likelihood of any
trajectory within that state space would be ![equation](http://latex.codecogs.com/gif.latex?1).

The proper likelihood is not of the trajectory within the discrete state
space; instead it’s the likelihood of the trajectory within the
continuous phase space, on which the discrete states are merely an
indicator function basis.

![equation](http://latex.codecogs.com/gif.latex?P%5Bx_%7B0...T-1%7D%5D%20dx%5EN%20%3D%20%5Cprod_%7Bi%3D0%7D%5E%7BT-1%7D%20T%28x_i%20%5Crightarrow%20x_%7Bi%2B1%7D%29%20%5Ccdot%20%5Cprod_%7Bi%3D0%7D%5E%7BT%7D%20p%28x_%7Bi%7D%20%7C%20%5Csigma%28x_%7Bi%7D%29%29)

With a discrete, non-overlapping state space, the likelihood of the
trajectory can be decomposed into a product over the trajectory of two
types of terms: the state to state transition probabilities and the
“emission” probabilities of each state, the probability of observing a
conformation at a given location in phase space given that the
conformation (![equation](http://latex.codecogs.com/gif.latex?x_t)) is within a certain state (![equation](http://latex.codecogs.com/gif.latex?%5Csigma%28x_t%29)).

Emission Distributions
======================

Up until this point, the emission probability distributions have not
been a model parameter for MSMs. Nevertheless, they are a critical
parameter from the point of view of evaluating trajectory likelihoods.
For example, consider two MSMs with the same transition probabilities.
However, in one case, the emission distributions for each state are
highly peaked at specific locations in phase space, whereas in the other
model the emission distributions are uniform over the volume of the
clusters. If the trajectory actually does go through the regions of high
likelihood in the first model, we would say that the first model more
accurately fits the data.

However, the long timescale behavior of the MSM, the quantity of most
interest scientifically, is independent of the choice of the emission
distributions. It’s only determined by the eigenspectrum of the
transition matrix elements. (The emission distributions determine the
estimated eigenfunctions, but not the eigenvectors.)

Therefore, the most appropriate emission distribution for discrete state
MSMs is that of the uniform distribution over the phase-space volume of
the state. That is, the likelihood of observing a conformation in phase
space given that the conformation is assigned to state ![equation](http://latex.codecogs.com/gif.latex?i) is ![equation](http://latex.codecogs.com/gif.latex?0) if the
conformation is outside of the bounding volume of the state and constant
if the conformation is within the volume. The constant is set so that
the distribution integrates to ![equation](http://latex.codecogs.com/gif.latex?1), and is thus the reciprocal volume of
the microstate.

![equation](http://latex.codecogs.com/gif.latex?P%5Bx_%7B0...T-1%7D%5D%20dx%5EN%20%3D%20%5Cprod_%7Bi%3D0%7D%5E%7BT-1%7D%20T%28x_i%20%5Crightarrow%20x_%7Bi%2B1%7D%29%20%5Ccdot%20%5Cprod_%7Bi%3D0%7D%5ET%20%5Cfrac%7B1%7D%7BV_%7B%5Csigma%28x_%7Bi%7D%29%7D%7D)

Algorithm
=========

To use this uniform distribution emission model, computationally, we
need to compute the volume of our MSM states, which are high-dimensional
Voronoi cells. While trivial in two or three dimensions, this
computational geometry task becomes challenging in large dimensional
settings. The computation of high dimensional volumes has occupied
significant attention in recent years in the computational geometry
literature, especially via randomized algorithms. [See issue \#1 on the
github]

One challenge is how to model the volume of states which are at the
“edge” of the model. Is the volume of these states infinite, extending
all the way out to infinity in some direction? This seems problematic.
Instead, it seems appropriate to assert that the volume is bounded by
the convex hull of the data. That is, the volume of a state might be
defined as the volume of the intersection of its Voronoi cell and the
convex hull of the dataset.

Instead, we use a slightly modified version of this definition that
adopts the same spirit. Instead of taking the outer bounding region to
be the convex hull of the data, we take it to be the set of all trail
points such that the nearest neighbor of the trial point within the
dataset is less than a certain cutoff, ![equation](http://latex.codecogs.com/gif.latex?R). This can be computed
relatively efficiently using a BallTree data structure. For further
efficiency, we might use only a random subsample of the dataset for this
nearest neighbor computation.

Instead of computing the volume of the Voronoi cells explicitly, we
instead compute the ratio of volume of each Voronoi cell to the volume
of the entire bounding region. Because the bounding region’s volume is
independent of the clustering parameters, its inclusion changes all of
the calculated likelihoods by a constant multiplicative factor, and can
thus be discarded from the perspective of model comparison.

To compute the fractional volumes of the Voronoi cells, we use the
following randomized algorithm. First, we generate points inside of the
bounding region using the lazy random walk Markov chain Monte Carlo
algorithm.\cite{Kannan97} This random walk converges to sampling from
the uniform distribution over the interior. After a sufficient burn in
period, for each sampled point we compute its nearest MSM state center,
assigning it to that state. As the number of randomly sampled points
goes to infinity, the fraction of the points assigned to each state
converges to being proportional to the state’s volume.

This algorithm is highly amenable to parallel computation. For a
euclidean distance metric, the assignment can be performed efficiently
by using a BallTree data structure for fast neighbor search.

Choosing the Optimal Number of States
=====================================

We choose both the clustering algorithm (k-centers, k-means, etc) and
the number of states by maximizing the BIC/AIC scores of the model,
using this likelihood.

Unfortunately, this doesn’t really help with picking the projection
operator to vectorize the conformations (e.g. dihedrals, contact maps,
etc). Also, it’s not going to work rigorously with RMSD.

Validation
==========

We being by using the procedure described above on the Müller potential.
