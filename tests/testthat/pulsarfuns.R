runtests <- function(pfun, pclass, dat, ...) {
    G <- dat$theta
    test_that("weird lambda path results in correct error or warning", {
        lams <- seq(.5, .7, length.out=5)
        expect_warning(out <- pfun(dat$data, fargs=list(lambda=lams, verbose=FALSE), rep.num=2,
                       ...), "lambda path")
        expect_warning(out <- pfun(dat$data, fargs=list(lambda=lams[1], verbose=FALSE), rep.num=2,
                       ...),"1 value")
        expect_error(out <- pfun(dat$data, fargs=list(lams=lams, verbose=FALSE), rep.num=2),
                      "missing")
        expect_warning(out <- pfun(dat$data, fargs=list(lambda=lams[4:5], verbose=FALSE), rep.num=2,
                       ...), "supplied path")
    })


    test_that("pulsar w/ lambda path works for huge", {
        ## run pulsar in serial mode
        lams <- getLamPath(.7, .1, 5)
        out <- pfun(dat$data, fargs=list(lambda=lams, verbose=FALSE), rep.num=2, ...)
        expect_is(out, pclass)
        expect_equal(out$stars$criterion, "stars.stability")
        # stars summary is monotonic increasing
        expect_equal(out$stars$summary, cummax(out$stars$summary))
        # merge objects dims match original graph, data
        expect_true(all(sapply(out$stars$merge, function(x) all(dim(x) == dim(G)))))
        expect_true(all(sapply(out$stars$merge, function(x) all(dim(x) == ncol(dat$data)))))

        ## run with gcd
        out <- pfun(dat$data, fargs=list(lambda=lams, verbose=FALSE), criterion=c("stars", "gcd"),
                    rep.num=2, ...)
        expect_equal(length(out$gcd$summary), length(out$stars$summary))
        expect_equal(out$gcd$criterion, "graphlet.stability")
        expect_error(fit <- refit(out, 'stars'), NA)
    })


    test_that("pulsar bounds are consistent", {
        ## run pulsar in serial mode
        lams <- getLamPath(.7, .1, 10) #exp(seq(log(.7), log(.05), length.out=20))
        ## run with gcd
        out   <- pfun(dat$data, fargs=list(lambda=lams, verbose=FALSE), criterion=c("stars"),
                      rep.num=5, ...)

        outlb <- update(out, lb.stars=TRUE)
        expect_equal(outlb$stars$opt.ind, out$stars$opt.ind) # same answer using bounds

        outb <- update(outlb, ub.stars=TRUE, criterion=c('stars', 'gcd'))
        ## gcd computed between bounds
        expect_equal(length(outb$gcd$summary), outb$stars$lb.index-outb$stars$ub.index+1)
        expect_equal(opt.index(outb, 'stars'), opt.index(out, 'stars')) # same answer using bounds

        ## check F1 score is OK
        print(getAnywhere("opt.index<-"))
        opt.index(outb, 'gcd') <- get.opt.index(outb, 'gcd')
        pdf(NULL)
        starsF1 <- huge::huge.roc(list(outb$stars$merge[[ opt.index(outb, 'stars') ]] > 0), G, verbose=FALSE)$F1
        gcdF1   <- huge::huge.roc(list(outb$stars$merge[[ opt.index(outb, 'gcd')   ]] > 0), G, verbose=FALSE)$F1
        dev.off()
        expect_gte(gcdF1, starsF1)
    })
}
