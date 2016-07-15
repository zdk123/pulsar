context("pulsar: serial mode")

p <- 20
## generate synthetic data
set.seed(10010)
dat <- huge::huge.generator(p*10, p, "random", prob=.25, verbose=FALSE, v=.1, u=.5)

source('pulsarfuns.R')
runtests(pulsar, "pulsar", dat, seed=10010)
