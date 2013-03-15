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

![equation](http://latex.codecogs.com/gif.latex?)\label{eq:like}
P[x_{0...T-1}] dx^N = \prod_{i=0}^{T-1} T(x_i \rightarrow x_{i+1}) \cdot \prod_{i=0}^{T} p(x_{i} | \sigma(x_{i}))![equation](http://latex.codecogs.com/gif.latex?)

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

![equation](http://latex.codecogs.com/gif.latex?)\label{eq:like_vol}
P[x_{0...T-1}] dx^N = \prod_{i=0}^{T-1} T(x_i \rightarrow x_{i+1}) \cdot \prod_{i=0}^T \frac{1}{V_{\sigma(x_{i})}}![equation](http://latex.codecogs.com/gif.latex?)

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
“edge” of the model. Is the volume of these states unbounded, given that
they extend all the way out to infinity in some direction? This would be
problematic. It seems appropriate to assert that the volume of these
edge states is bounded in some way by the extent of our dataset. For
example, the volume of a state might be defined as the volume of the
intersection of its Voronoi cell and the convex hull of the whole
dataset.

We use a slightly modified version of this definition that adopts the
same spirit. Instead of taking the outer bounding region to be the
convex hull of the data, we take it to be the set of all trial points
such that the nearest data point to the trial point is closer than a
certain cutoff, ![equation](http://latex.codecogs.com/gif.latex?R). This can be computed relatively efficiently using a
BallTree data structure. For further efficiency, we might use only a
random subsample of the dataset for this nearest neighbor computation.

Instead of computing the volume of the Voronoi cells explicitly, we
instead compute the ratio of volume of each Voronoi cell to the volume
of the entire bounding region. Because the bounding region’s volume is
independent of the clustering parameters, its inclusion changes all of
the calculated likelihoods by a constant multiplicative factor, and can
thus be discarded from the perspective of model comparison.

To compute the fractional volumes of the Voronoi cells, we use the
following randomized algorithm. First, we generate points inside of the
bounding region using the lazy random walk Markov chain Monte Carlo
algorithm.@Kannan97 This random walk converges to sampling from the
uniform distribution over the interior. After a sufficient burn in
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

Langevin Dynamics on the Müller Potential
-----------------------------------------

We being by using the procedure described above on data generated by
simulating Langevin dynamics on the Müller potential (). This data was
clustered into models of different sizes using one of two algorithms.
The first, K-Centers, tends to make equal sized states as it works to
minimize the maximum distance between a conformation and its generator.
The second, Hybrid, will make varying size states because it works to
minimize the average distance between conformations and their
generators. In both cases, we built models ranging from 50 to 1,000
total states. Since this space is only two-dimensional, we used the
simple, brute-force MC algorithm for calculating volumes.

Briefly, we generated uniform random vectors within a bounding box whose
size was slightly larger than the data. This random sample was then
accepted if it was ’close enough’ to a generator in the 100 state model,
where close enough meant it was within 0.2 units from any generator in
that model. This was used to approximate the convex hull of the data.
Then the accepted points were assigned to a state. Then the volume of a
state is proportional to the number of samples assigned to that state
divided by the total number of samples. We note, that these volumes are
actually in units of the total acceptable volume. This affects the
absolute value of the likelihood, however, since the total volume was
constant in all cases, the relative likelihoods are the same.

[h!] ![The Müller potential consists of three metastable wells. We
simulated Langevin dynamics on this potential as a validation of the new
likelihood function.](figs/muller_pot.png "fig:") [fig:muller~p~ot]

We built an MSM using a fixed lag time for each state decomposition and
then calculated the likelihood according to . These likelihoods
increased rapidly with ![equation](http://latex.codecogs.com/gif.latex?k), but then plateaued at a value of around 300
states . At this plateau, we can conclude that we do not gain any new
information by adding new states to our model and so should pick ![equation](http://latex.codecogs.com/gif.latex?k) to
be this point. This is a heuristic argument that can be made more
mathematically rigorous with the introduction of functions like the
Bayesian Information Content (BIC), which is defined in .

![equation](http://latex.codecogs.com/gif.latex?)BIC = -2 \log L + m \log(n) 
\label{eq:bic}![equation](http://latex.codecogs.com/gif.latex?)

Here, ![equation](http://latex.codecogs.com/gif.latex?L) denotes the likelihood of a model, while ![equation](http://latex.codecogs.com/gif.latex?m) is the number of
parameters used in the model, and ![equation](http://latex.codecogs.com/gif.latex?n) is the number of data points.
Essentially this is a mathematical way of saying we would like to use a
model with the maximum likelihood, but without introducing too many
parameters. We found that the BIC of the models built using K-Centers
clustering has a minimum at ![equation](http://latex.codecogs.com/gif.latex?k%3D250). This minimum is the “optimal” model
according to the BIC of our likelihood. As a check of the usefulness of
this minimum, we plotted the three slowest implied timescales in the
models as a function of ![equation](http://latex.codecogs.com/gif.latex?k). The optimal model was close to the point at
which the implied timescales became constant.

![For models built between 50 and 1,000 states using the K-Centers
algorithm, the log likelihood function increased quickly and plateaued
at approximately 300 states. Heuristically, at this point, adding
another state does not add anything to the model. As a result, by using
a criterion, like the BIC, that penalizes using additional parameters to
fit the data, we can find an “optimal” model. Interestingly, the BIC
optimal model occurs at 250 states, which also corresponds to where the
three slowest timescales plateau.](figs/like_comp.pdf "fig:")
[fig:like~k~centers]
