<!-- README.md is generated from README.Rmd. Please edit that file -->



# pulsar: Paralellized Utilities for lambda Selection Along a Regularization path.

This package is an implementation of the Stability Approach to Regularization Selection 
([StARS, H. Liu - 2010](https://papers.nips.cc/paper/3966-stability-approach-to-regularization-selection-stars-for-high-dimensional-graphical-models.pdf)),
with options for speedups, paralellizations and alternate/complementary metrics 
(most notably, graphlet stability). This package is intended for selection of L1-regularized graphical
models (e.g. glasso or Meinshausen-Buhlmann neighborhood selection), where the level of sparsity is
typically parametrized by *lambda*. This R package uses function passing to allow you to use your
favorite implementation (e.g. huge, QUIC, etc) for graphical model learning along a lambda
regularization path.The method is quite generic and could be used for other Lasso-type problems as
well.

Additionally, this package implements an option to find lower/upper bounds on the StARS-selected 
lambda, using some heuristics that work well in practice. This allows us to cut off a large chunk of
the lambda path after after only N=2 subsamples. 
This works well particularly for sparse graphs when the dimensionality is high (above 20 variables
or so).  Thus, computational burden is greatly reduced, especially when many subsamples are required
and even when running in embarassingly parallel (batch) mode. We also use these bounds to narrow the
range that we look for relative stable points over an induced-subgraph (graphlet) stability metric.

In batch computing systems, we use the [BatchJobs](https://cran.r-project.org/web/packages/BatchJobs/)
Map/Reduce strategy (for batch computing systems such as Torque, LSF, SLURM or SGE) which can 
significantly reduce the computation and memory burdens for StARS. This is useful for hpc users, 
when the number of processors available on a multicore machine only allows modest parallelization.

Please see the paper preprint on [arXiv](https://arxiv.org/abs/1605.07072).


## Installation

Installation of the current development version of pulsar requires the devtools package. Also we
recommend some other packages for installation as needed (such as `huge`).


```r
library(devtools)
install_github("zdk123/pulsar")
library(pulsar)
```

## Basic usage

For this readme, we will use synthetic data, generated from the `huge` package.


```r
library(huge)
set.seed(10010)
p <- 40 ; n <- 1200
dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
lams  <- getLamPath(.2, .01, len=40)
```

You can use the `pulsar` package to run StARS, serially, as a drop-in replacement for the 
`huge.select` function in the `huge` package. Pulsar differs in that we run the model selection step
first and then refit using arguments stored in the original call. Remove the `seed` argument for
real data (this seeds the pseudo-random number generator to fix subsampling for reproducing test
code).


```r
hugeargs <- list(lambda=lams, verbose=FALSE)
time1    <- system.time(
out.p    <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
                rep.num=20, criterion='stars', seed=10010))
fit.p    <- refit(out.p)
```

Inspect the output:

```r
out.p
# Mode: serial
# Path length: 40 
# Subsamples:  20 
# Graph dim:   40 
# Criterion:
#   stars... opt: index 15, lambda 0.132
fit.p
# Pulsar-selected refit of huge::huge 
# Path length: 40 
# Graph dim:   40 
# Criterion:
#   stars... sparsity 0.0325
```

Including the lower bound option `lb.stars` and upper bound option `ub.stars` can improve runtime
for same StARS result.


```r
time2 <- system.time(
out.b <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
                rep.num=20, criterion='stars', seed=10010,
                lb.stars=TRUE, ub.stars=TRUE))
```

Compare runtimes and StARS selected lambda index for each method.


```r
time2[[3]] < time1[[3]]
# [1] TRUE
opt.index(out.p, 'stars') == opt.index(out.b, 'stars')
# [1] TRUE
```

## Using a custom graphical model method

You can pass in an arbitrary graphical model estimation function to fun. The function has some
requirements: the first argument must be the nxp data matrix, and one argument must be named lambda,
which should be a decreasing numeric vector containing the lambda path. 
The output should be a list of adjacency matrices (which can be of sparse representation from the
`Matrix` package to save memory).
Here is an example from `QUIC`.


```r
library(QUIC)
quicr <- function(data, lambda) {
    S    <- cov(data)
    est  <- QUIC::QUIC(S, rho=1, path=lambda, msg=0, tol=1e-2)
    path <-  lapply(seq(length(lambda)), function(i) {
                tmp <- est$X[,,i]; diag(tmp) <- 0
                as(tmp!=0, "lgCMatrix")
    })
    est$path <- path
    est
}
```

We can use pulsar with a similar call. We can also parallelize this a bit for multi-processor
machines by specifying ncores (which wraps mclapply in the parallel package).

```r
quicargs <- list(lambda=lams)
out.q <- pulsar(dat$data, fun=quicr, fargs=quicargs,
                rep.num=100, criterion='stars', seed=10010,
                lb.stars=TRUE, ub.stars=TRUE, ncores=2)
```

## Graphlet stability

We can use the graphlet correlation distance as an additional stability criterion. We could call
pulsar again with a new criterion, or simply `update` the arguments for model we already used. Then,
we can use our default approach for selecting the optimal index, based on the gcd+StARS criterion:
choose the minimum gcd summary statistic between the upper and lower StARS bouds.


```r
out.q2 <- update(out.q, criterion=c('stars', 'gcd'))
opt.index(out.q2, 'gcd') <- get.opt.index(out.q2, 'gcd')
fit.q2 <- refit(out.q2)
```

Compare model error by relative hamming distances between refit adjacency matrices and the true
graph and visualize the results:

















