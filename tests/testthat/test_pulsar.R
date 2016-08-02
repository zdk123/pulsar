source('pulsarfuns.R')

rseed <-  10010 #sample.int(1000, 1) #
p <- 10
## generate synthetic data
set.seed(rseed)
dat <- huge::huge.generator(p*100, p, "hub", verbose=FALSE, v=.1, u=.4)

library(BatchJobs)
options(BatchJobs.verbose=FALSE)
suppressPackageStartupMessages(library(BatchJobs))
conffile <- file.path(system.file(package="pulsar"), "extdata", "BatchJobsSerialTest.R")
unlink(getTempDir(), recursive=TRUE)

######################################################
context("pulsar: huge, serial mode")
huge.serial <- runtests(pulsar, "pulsar", dat, fun=huge::huge, fargs=list(verbose=FALSE, scr=TRUE), seed=rseed)

######################################################
context("pulsar: quic, serial mode")
quic.serial <- runtests(pulsar, "pulsar", dat, fun=quicr, fargs=list(), seed=rseed)

######################################################
context("pulsar: huge, batch mode")
huge.batch <- runtests(batch.pulsar, "batch.pulsar", dat, fun=huge::huge, fargs=list(verbose=FALSE, scr=TRUE),
                       conffile=conffile, progressbars=FALSE, cleanup=TRUE, seed=rseed)

######################################################
context("pulsar: quic, batch mode")
quic.batch <- runtests(batch.pulsar, "batch.pulsar", dat, fun=quicr, fargs=list(),
                conffile=conffile, progressbars=FALSE, cleanup=TRUE, seed=rseed)

######################################################
context("pulsar: serial vs batch")
msg  <- "huge: serial and batch mode are equivilent: no bounds"
msg2 <- "huge: serial and batch mode are equivilent: lower bound"
runcomptest(msg,  huge.serial$out,   huge.batch$out)
runcomptest(msg2, huge.serial$outb, huge.batch$outb)

msg  <- "quic: serial and batch mode are equivilent: no bounds"
msg2 <- "quic: serial and batch mode are equivilent: lower bound"
runcomptest(msg,  quic.serial$out,   quic.batch$out)
runcomptest(msg2, quic.serial$outb, quic.batch$outb)

#######################################################
context("refit estimation function")

testrefit("refitting pulsar gives correct warnings & output",        huge.serial$outb)
testrefit("refitting batch pulsar gives correct warnings & output",  huge.batch$outb)

testrefit("refitting pulsar gives correct warnings & output",        quic.serial$outb)
testrefit("refitting batch pulsar gives correct warnings & output",  quic.batch$outb)
