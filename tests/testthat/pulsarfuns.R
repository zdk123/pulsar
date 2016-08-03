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


runtests <- function(pfun, pclass, dat, fun, fargs, ...) {
    G <- dat$theta
    test_that("bad crits results in right errors", {
        lams <- getLamPath(.7, .1, 5)
        hargs <- c(fargs, list(lambda=lams))
        expect_error(out <- pfun(dat$data, fun=fun, fargs=hargs, rep.num=2, criterion=c("stars", "foo"),
                       ...), "foo")
        expect_error(out <- pfun(dat$data, fun=fun, fargs=hargs, rep.num=2, 
                       criterion=c("estrada", "sufficiency")))
    })


    test_that("weird lambda path results in correct error or warning", {
        lams <- seq(.5, .7, length.out=5)
        hargs <- c(fargs, list(lambda=lams))
        expect_warning(out <- pfun(dat$data, fun=fun, fargs=hargs, rep.num=2,
                       ...), "lambda path")
        expect_warning(out <- pfun(dat$data, fun=fun, fargs=c(list(lambda=lams[1]), fargs), rep.num=2,
                       ...),"1 value")
        expect_error(out <- pfun(dat$data, fun=fun, fargs=c(list(lams=lams), fargs), rep.num=2), "missing")
        expect_warning(out <- pfun(dat$data, fun=fun, fargs=c(list(lambda=lams[c(5,4)]), fargs), rep.num=2,
                       ...), "supplied values")
    })

    mlam <- signif(getMaxCov(dat$sigmahat)+.1, 3)
    lams  <- getLamPath(mlam, .05, 15)
    hargs <- c(fargs, list(lambda=lams))
    out   <- pfun(dat$data, fun=fun, fargs=hargs, criterion=c("stars"), rep.num=4, ...)
    outb  <- update(out, lb.stars=TRUE, ub.stars=TRUE, criterion=c("stars", "gcd"))


    test_that("pulsar w/ lambda path works for fun", {
        ## run pulsar in serial mode
        expect_is(out, pclass)
        expect_equal(out$stars$criterion, "stars.stability")
        # stars summary is monotonic increasing
        expect_equal(out$stars$summary, cummax(out$stars$summary))
        # merge objects dims match original graph, data
        expect_true(all(sapply(out$stars$merge, function(x) all(dim(x) == dim(G)))))
        expect_true(all(sapply(out$stars$merge, function(x) all(dim(x) == ncol(dat$data)))))

        ## run with gcd
##        out <- pfun(dat$data, fun=fun, fargs=hargs, criterion=c("stars", "gcd"), rep.num=2, ...)
    })

    test_that("pulsar bounds are consistent", {
        ## check lengths
##        expect_equal(length(out$gcd$summary), length(out$stars$summary))
        expect_equal(outb$gcd$criterion, "graphlet.stability")
        expect_error(fit <- refit(out, 'stars'), NA)
        expect_equal(outb$stars$opt.ind, out$stars$opt.ind) # same answer using bounds
        ## gcd computed between bounds
        expect_equal(length(outb$gcd$summary), outb$stars$lb.index-outb$stars$ub.index+1)
        expect_equal(opt.index(outb, 'stars'), opt.index(out, 'stars')) # same answer using bounds

        ## check F1 score is OK
        opt.index(outb, 'gcd') <- get.opt.index(outb, 'gcd')
        pdf(NULL)
        starsF1 <- huge::huge.roc(list(outb$stars$merge[[ opt.index(outb, 'stars') ]] > 0), G, verbose=FALSE)$F1
        gcdF1   <- huge::huge.roc(list(outb$stars$merge[[ opt.index(outb, 'gcd')   ]] > 0), G, verbose=FALSE)$F1
        dev.off()
        expect_gte(gcdF1, starsF1)
    })
    return(list(out=out, outb=outb))
}


runcomptest <- function(msg, out, out.batch, ...) {
    test_that(msg, {
        expect_gt(max(out$stars$summary), 0) # make sure summary isn't trivally zero
        expect_equivalent(out$stars$summary, out.batch$stars$summary)
    })
}

testrefit <- function(desc, outb) {
    test_that(desc, {
        expect_message(fit1 <- refit(outb, "stars"), regexp = NA)
        expect_message(fit2 <- refit(outb, "gcd"), "No optimal index")
        expect_equal(names(fit1$refit), "stars")
        expect_error(opt.index(outb, 'gcd') <- -1, "Index value")
        expect_error(opt.index(outb, 'gcd') <- get.opt.index(outb, 'gcd'), NA)
        expect_equal(opt.index(outb, 'gcd'), outb$gcd$opt.index)
        expect_equal(opt.index(outb, 'gcd'), get.opt.index(outb, 'gcd'))
        expect_warning(fit3 <- refit(outb), regexp = NA)

        expect_gt(sum(fit3$refit$stars), 0)
        expect_gt(sum(fit3$refit$gcd),   0)
        expect_warning(fit4 <- refit(outb, "foo"), "Unknown criterion")
    })
}
