This file was automatically generated from latex using pandoc

Likelihood framework for choosing the optimal number of states in a
Markov State Model

Overview
========

A significant challenge in the automated construction of markov state
models is choosing the number of microstates. The number of microstates,
![equation](http://latex.codecogs.com/gif.latex?k governs a bias-variance tradeoff – when http://latex.codecogs.com/gif.latex?k is large, the models
)have a lot of parameters enabling it to fit the data with low bias, but
the simultaneous estimation of these parameters from finite data can
![equation](lead to high variance. On the other hand, when http://latex.codecogs.com/gif.latex?k is small, we have
)fewer parameters and better statistics, but at the expense of model
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
![equation](probability in that model would be http://latex.codecogs.com/gif.latex?1 and thus the likelihood of any
)![equation](trajectory within that state space would be http://latex.codecogs.com/gif.latex?1.
)
The proper likelihood is not of the trajectory within the discrete state
space; instead it’s the likelihood of the trajectory within the
continuous phase space, on which the discrete states are merely an
indicator function basis.

![equation](http://latex.codecogs.com/gif.latex?P%5Bx%28t%29%5D%20dx%5EN%20%3D%20O_%7Bs%3D0%7D%28x_0%29%20%5Cprod_%7Bi%3D0%7D%5E%7BN-1%7D%20T%28x_i%20%5Crightarrow%20x_%7Bi%2B1%7D%29%20%5Ccdot%20O_%7Bs%3Di%7D%28x_i%29
)
With a discrete, non-overlapping state space, the likelihood of the
trajectory can be decomposed into a product over the trajectory of two
types of terms: the state to state transition probabilities and the
“emission” probabilities of each state, the probability of observing a
conformation at a given location in phase space given that the
conformation is within a certain state.

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

Therefore, the most appropriate emmision distribution for discrete state
MSMs is that of the uniform distribution over the phase-space volume of
the state. That is, the likelihood of observing a conformation in phase
![equation](space given that the conformation is assigned to state http://latex.codecogs.com/gif.latex?i is http://latex.codecogs.com/gif.latex?i if the
)conformation is outside of the bounding volume of the state and constant
if the conformation is within the volume. The constant is set so that
![equation](the distribution integrates to http://latex.codecogs.com/gif.latex?1, and is thus the reciprocal volume of
)the microstate.

Algorithm
=========

To use this uniform distribution emission model, computationally, we
need to compute the volume of our MSM states, which are high-dimensional
Voronoi cells. While trivial in two or three dimensions, this
computational geometry task becomes challenging in large dimensional
settings.

One challenge is how to model the volume of states which are at the
“edge” of the model. Is the volume of these states infinite, extending
all the way out to infinity in some direction? This seems problematic.
Instead, it seems appropriate to assert that the volume is bounded by
the convex hull of the data. That is, the volume of a state is defined
as the volume of the intersection of its Voronoi cell and the convex
hell of the dataset.

Computing both N-dimensional Voronoi volumes and convex hulls is
difficult too. We adopt the following approximate algorithm. First, we
compute the axis-aligned bounding box of all of the conformations. Next,
we fit a multi-dimensional gaussian to the coordinates of each state.
(We should probably use a diagonal covariance matrix, or even a
spherical gaussian?).

Now, we sample random points from the uniform distribution over the
bounding box. We check these points against all of the gaussian PDFs. If
the point is farther than `3-sigma` away from every cluster center, we
reject it, reasoning that it is probably outside the convex hull. The
remaining points are now a representation of the uniform distribution
over the convex hull of phase space. For each point, we compute its
distance to each of the cluster centers, “assigning” it to the center it
is closest to. As the number of randomly sampled points goes to
infinity, the fraction of the points assigned to each state converges to
being proportional to the state’s volume.

This algorithm is highly amenable to parallel computation.

Choosing the Optimal Number of States
=====================================

Choose both the clustering algorithm (k-centers, k-means, etc) and the
number of states by maximizing the BIC score of the model, using this
likelihood.