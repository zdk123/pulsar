context("setup pulsar test")

suppressPackageStartupMessages(library(batchtools))
options(batchtools.verbose=FALSE)
options(batchtools.progress=FALSE)
source('pulsarfuns.R')

rseed <- 10010
p     <- 15
set.seed(rseed)
dat <- huge::huge.generator(p*30, p, "hub", verbose=FALSE, v=.4, u=.2)
set.seed(rseed)
dat$data <- MASS::mvrnorm(p*30, mu=rep(0,p), Sigma=dat$sigma, empirical=TRUE)

tmpdir <- fs::path_real(tempdir())
conffile <- ""
pargs <- list(verbose=FALSE, scr=TRUE)

cstr <- "pulsar: %s; mode: %s"
######################################################
context(sprintf(cstr, 'huge', 'serial'))
huge.serial       <- runtests(pulsar, "pulsar", dat, fun=huge::huge,
                        fargs=pargs, seed=rseed, refit=FALSE)
huge.serial.refit <- runtests(pulsar, "pulsar", dat, fun=huge::huge,
                        fargs=pargs, seed=rseed, refit=TRUE)

######################################################
context(sprintf(cstr, 'clime', 'serial'))
clime.serial       <- runtests(pulsar, "pulsar", dat, fun=climer,
                              fargs=list(seed=rseed), seed=rseed, refit=FALSE)
clime.serial.refit <- runtests(pulsar, "pulsar", dat, fun=climer,
                              fargs=list(seed=rseed), seed=rseed, refit=TRUE)

######################################################
context(sprintf(cstr, 'huge', 'batch'))
huge.batch       <- runtests(batch.pulsar, "batch.pulsar", dat, fun=huge::huge,
                       fargs=pargs, conffile=conffile, cleanup=TRUE,
                       seed=rseed, wkdir=tmpdir, refit=FALSE)
huge.batch.refit <- runtests(batch.pulsar, "batch.pulsar", dat, fun=huge::huge,
                       fargs=pargs, conffile=conffile, cleanup=TRUE,
                       seed=rseed, wkdir=tmpdir, refit=TRUE)

######################################################
context(sprintf(cstr, 'clime', 'batch', ''))
clime.batch       <- runtests(batch.pulsar, "batch.pulsar", dat, fun=climer,
                             fargs=list(seed=rseed), conffile=conffile, cleanup=TRUE,
                             seed=rseed, wkdir=tmpdir, refit=FALSE)
clime.batch.refit <- runtests(batch.pulsar, "batch.pulsar", dat, fun=climer,
                             fargs=list(seed=rseed), conffile=conffile, cleanup=TRUE,
                             seed=rseed, wkdir=tmpdir, refit=TRUE)

######################################################
context("pulsar: serial vs batch")
msg  <- "huge: serial and batch mode are equivilent: no bounds"
msg2 <- "huge: serial and batch mode are equivilent: lower bound"
runcomptest(msg,  huge.serial$out,   huge.batch$out)
runcomptest(msg2, huge.serial$outb, huge.batch$outb)

msg  <- "clime: serial and batch mode are equivilent: no bounds"
msg2 <- "clime: serial and batch mode are equivilent: lower bound"
runcomptest(msg,  clime.serial$out,   clime.batch$out)
runcomptest(msg2, clime.serial$outb, clime.batch$outb)

#######################################################
context("pulsar: refit mode")
msg  <- "huge: serial refit test: no bounds"
msg2 <- "huge: serial refit test: lower bound"
runcomptest(msg, huge.serial$out,  huge.serial.refit$out)
runcomptest(msg, huge.serial$outb, huge.serial.refit$outb)

msg  <- "huge: batch refit test: no bounds"
msg2 <- "huge: batch refit test: lower bound"
runcomptest(msg, huge.batch$out,  huge.batch.refit$out)
runcomptest(msg, huge.batch$outb, huge.batch.refit$outb)

msg  <- "clime: serial refit test: no bounds"
msg2 <- "clime: serial refit test: lower bound"
runcomptest(msg, clime.serial$out,  clime.serial.refit$out)
runcomptest(msg, clime.serial$outb, clime.serial.refit$outb)

msg  <- "clime: batch refit test: no bounds"
msg2 <- "clime: batch refit test: lower bound"
runcomptest(msg, clime.batch$out,  clime.batch.refit$out)
runcomptest(msg, clime.batch$outb, clime.batch.refit$outb)


#######################################################
context("refit estimation function")
testrefit("refitting pulsar gives correct warnings & output", huge.serial$outb)
testrefit("refitting batch pulsar gives correct warnings & output", huge.batch$outb)

testrefit("prefitting pulsar gives correct warnings & output", huge.serial.refit$outb)
testrefit("prefitting batch pulsar gives correct warnings & output", huge.batch.refit$outb)

testrefit("refitting pulsar gives correct warnings & output", clime.serial$outb)
testrefit("refitting batch pulsar gives correct warnings & output", clime.batch$outb)

testrefit("prefitting pulsar gives correct warnings & output", clime.serial.refit$outb)
testrefit("prefitting batch pulsar gives correct warnings & output", clime.batch.refit$outb)
