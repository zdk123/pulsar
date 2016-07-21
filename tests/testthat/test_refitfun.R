
#p <- 8
#### generate synthetic data
#set.seed(10010)
#dat <- huge::huge.generator(p*20, p, "hub", verbose=FALSE, v=.25, u=.4)
#G <- dat$theta
#lams <- getLamPath(.4, .1, 5)
#lb <- TRUE
#hargs <- list(lambda=lams, verbose=FALSE, scr=TRUE)
#out.p <- pulsar(dat$data, fargs=hargs, rep.num=5, criterion=c('stars', 'gcd'), seed=10010, 
#                lb.stars=lb, ub.stars=TRUE)

#suppressPackageStartupMessages(library(BatchJobs))
#options(BatchJobs.verbose=FALSE)
#conffile <- file.path(system.file(package="pulsar"), "extdata", "BatchJobsSerialTest.R")
#out.bp <- batch.pulsar(dat$data, fargs=hargs,  rep.num=4, criterion=c('stars', 'gcd'), 
#                  lb.stars=lb, ub.stars=TRUE, conffile=conffile, progressbars=FALSE, 
#                  cleanup=TRUE, seed=10010)


#testrefit <- function(desc, out) {
#    test_that(desc, {
#        expect_warning(fit1 <- refit(out, "stars"), regexp = NA)
#        expect_warning(fit2 <- refit(out, "gcd"), "No optimal index")
#        expect_equal(names(fit1$refit), "stars")
#        expect_error(opt.index(out, 'gcd') <- -1, "Index value")
#        expect_error(opt.index(out, 'gcd') <- get.opt.index(out, 'gcd'), NA)
#        expect_equal(opt.index(out, 'gcd'), out$gcd$opt.index)
#        expect_equal(opt.index(out, 'gcd'), get.opt.index(out, 'gcd'))
#        expect_warning(fit3 <- refit(out), regexp = NA)

#        expect_gt(sum(fit3$refit$stars), 0)
#        expect_gt(sum(fit3$refit$gcd),   0)
#        expect_warning(fit4 <- refit(out, "foo"), "Unknown criterion")

#        out2 <- update(out, lb.stars=FALSE)
#        expect_error(get.opt.index(out2, 'gcd'), "Lower bound needed")
#    })
#}

#testrefit("refitting pulsar gives correct warnings & output",        out.p)
#testrefit("refitting batch pulsar gives correct warnings & output", out.bp)
