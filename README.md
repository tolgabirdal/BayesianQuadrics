# BayesianQuadrics

## Introduction
The code in this repository was created to fit quadrics to point clouds which have large portions of noisy, or missing data. In this scenario fitting an unconstrained quadric could lead to an undesirable surface, for example, a hyperboloid, flattened ellipsoid or plane. One way to constrain the surface is to fit it using a probabilistic model with a Bayesian prior over the surface parameters.

More information can be found in the paper,
 * ["Fitting quadrics with a Bayesian prior" - Daniel Beale et al. Journal of Computational Visual Media, 2016.](http://link.springer.com/article/10.1007/s41095-016-0041-9)

**Quadric fitting examples**

The image below gives an example of some data extracted on the roof of a phone booth, on which no points were created in the reconstruction process. It is clear that the surface forms is drawn from an ellipsoid given that we know that the object is a phone booth, but fitting the quadric without a prior leads to the flattened ellipsoid shown in Fit 1. As we modify the prior, the fit becomes closer and closer to a sphere, and finally gives an accurate representation of the surface as shown in Fit 3.

<img alt="Quadric fitting examples" style="float:right" width="500px" src="doc/QuadricFitExamples.png" />

The library documentation can be found at [dabeale.github.io/BayesianQuadrics](http://dabeale.github.io/BayesianQuadrics).

## Installation
First clone the git repository using the command,
> git --recursive clone https://github.com/dabeale/BayesianQuadrics

Note that the '--recursive' flag is required to clone the Eigen dependency, which is included at a submodule. If Eigen is already installed on your machine, edit the 'INCLUDES' variable in 'makefile'.

GNU Octave must be installed in order to compile the software. It can be downloaded from [www.gnu.org/software/octave/](https://www.gnu.org/software/octave/), or installed using homebrew for mac. The statistics package is required (`pkg -forge install statistics`).

Enter the 'BayesianQuadrics' directory and run
> make

## Tests
From octave, enter the scripts directory and run,
> testQuadricFit

The code will generate an Output directory and a collection of ply files containing the points generated on the test quadric, and the fitted surfaces.

## Matlab compatibility
The code is compatible with matlab. If octave is installed you can run,
> make mex

to generate a mexfile for the c++ implementation.

