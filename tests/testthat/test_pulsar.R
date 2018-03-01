context("setup")

suppressPackageStartupMessages(library(batchtools))
options(batchtools.verbose=FALSE)
options(batchtools.progress=FALSE)
source('pulsarfuns.R')

rseed <- 10010
p     <- 28
set.seed(rseed)
dat <- huge::huge.generator(p*100, p, "hub", verbose=FALSE, v=.4, u=.2)
set.seed(rseed)
dat$data <- MASS::mvrnorm(p*100, mu=rep(0,p), Sigma=dat$sigma, empirical=TRUE)

fp    <- batchtools:::fp
npath <- function(x) normalizePath(x, mustWork=FALSE, winslash="/")
homedir <- npath(fp("~"))

conffile <- ""

######################################################
context("pulsar: huge, serial mode")
huge.serial <- runtests(pulsar, "pulsar", dat, fun=huge::huge, fargs=list(verbose=FALSE, scr=TRUE), seed=rseed)

######################################################
context("pulsar: quic, serial mode")
quic.serial <- runtests(pulsar, "pulsar", dat, fun=quicr, fargs=list(), seed=rseed)

######################################################
context("pulsar: huge, batch mode")
huge.batch <- runtests(batch.pulsar, "batch.pulsar", dat, fun=huge::huge,
                fargs=list(verbose=FALSE, scr=TRUE), conffile=conffile,
                cleanup=TRUE, seed=rseed, wkdir=homedir)

# ######################################################
context("pulsar: quic, batch mode")
quic.batch <- runtests(batch.pulsar, "batch.pulsar", dat, fun=quicr,
                    fargs=list(), conffile=conffile, cleanup=TRUE,
                    seed=rseed, wkdir=homedir)

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
