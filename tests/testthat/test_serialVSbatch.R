#context("pulsar: serial vs batch")
### Make sure serial and batch mode results in identical solutions when the same seed is used

#options(BatchJobs.verbose=FALSE)
#suppressPackageStartupMessages(library(BatchJobs))

#p <- 8
#set.seed(10010)
#dat <- huge::huge.generator(p*10, p, "random", prob=.25, verbose=FALSE, v=.1, u=.5)
#G <- dat$theta
#conffile <- file.path(system.file(package="pulsar"), "extdata", "BatchJobsSerialTest.R")
### get a random, consistent seed
#rseed <- sample.int(1000, 1)

#msg  <- "serial and batch mode are equivilent: no bounds"
#msg2 <- "serial and batch mode are equivilent: lower bound"

#runtest <- function(msg, ...) {
#    test_that(msg, {
#        ## run pulsar in serial mode
#        lams <- getLamPath(.5, .05, 5)

#        out       <- pulsar(dat$data, fargs=list(lambda=lams, verbose=FALSE, scr=TRUE),
#                         criterion=c("stars"), rep.num=4, seed=rseed, ...)
#        out.batch <- batch.pulsar(dat$data, fargs=list(lambda=lams, verbose=FALSE, scr=TRUE),
#                         criterion=c("stars"), rep.num=4, seed=rseed,
#                         conffile=conffile, progressbars=FALSE, cleanup=TRUE, ...)

#        expect_gt(max(out$stars$summary), 0) # make sure summary isn't trivally zero
#        expect_equivalent(out$stars$summary, out.batch$stars$summary)
#    })
#}

#runtest(msg)
#runtest(msg2, lb.stars=TRUE)
