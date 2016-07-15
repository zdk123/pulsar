context("pulsar: batch mode")

options(BatchJobs.verbose=FALSE)
suppressPackageStartupMessages(library(BatchJobs))
p <- 20
## generate synthetic data
set.seed(10010)
dat <- huge::huge.generator(p*10, p, "random", prob=.25, verbose=FALSE, v=.1, u=.5)

source('pulsarfuns.R')
conffile <- file.path(system.file(package="pulsar"), "extdata", "BatchJobsSerialTest.R")
runtests(batch.pulsar, "batch.pulsar", dat, conffile=conffile, progressbars=FALSE, cleanup=TRUE, seed=10010)
