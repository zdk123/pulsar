context("failed subsample fits")

## Correlation thresholding estimator that can report failures.
## A subsample "fails" at the smallest lambdas when its first column sums
## above zero, so some, but not all, subsamples fail.
p <- 8
set.seed(10010)
dat <- huge::huge.generator(200, p, "hub", verbose=FALSE, v=.4, u=.2)
X   <- dat$data
lams <- getLamPath(.9, .02, 15)

corr_path <- function(data, lambda) {
  S <- cor(data)
  lapply(lambda, function(lam) {
    tmp <- abs(S) > lam ; diag(tmp) <- FALSE
    Matrix::Matrix(tmp, sparse=TRUE)
  })
}
fails_at <- function(data, lambda, nfail) {
  ## tie failures to lambda values, so they are stable over sub-paths (lb.stars)
  if (sum(data[,1]) > 0) which(lambda %in% utils::tail(lams, nfail)) else integer(0)
}
## report via NULL path entries
fun_null <- function(data, lambda, nfail=3) {
  path <- corr_path(data, lambda)
  path[fails_at(data, lambda, nfail)] <- list(NULL)
  list(path=path)
}
## report via a 'failed' vector, returning an empty graph at failures
fun_flag <- function(data, lambda, nfail=3) {
  path   <- corr_path(data, lambda)
  failed <- seq_along(lambda) %in% fails_at(data, lambda, nfail)
  path[failed] <- lapply(path[failed], function(x) x * 0)
  list(path=path, failed=failed)
}
## silent failure: empty graph, not reported
fun_silent <- function(data, lambda, nfail=3) {
  path   <- corr_path(data, lambda)
  failed <- seq_along(lambda) %in% fails_at(data, lambda, nfail)
  path[failed] <- lapply(path[failed], function(x) x * 0)
  list(path=path)
}

quiet <- function(expr) suppressWarnings(expr)

test_that(".fitpath marks failures as NULL", {
  m <- Matrix::Matrix(diag(2), sparse=TRUE)
  expect_equal(pulsar:::.fitpath(list(path=list(m, m))), list(m, m))
  expect_equal(pulsar:::.fitpath(list(path=list(m, NULL, NA))), list(m, NULL, NULL))
  expect_equal(pulsar:::.fitpath(list(path=list(m, m), failed=c(FALSE, TRUE))),
               list(m, NULL))
  expect_equal(pulsar:::.fitpath(list(path=list(m, m), failed=c(NA, FALSE))),
               list(NULL, m))
  expect_error(pulsar:::.fitpath(list(path=list(m, m), failed=TRUE)), "same length")
})

test_that("stars frequencies use successful fits only", {
  full  <- Matrix::Matrix(1 - diag(3), sparse=TRUE)
  empty <- full * 0
  ## lambda 2: two complete graphs and two failures
  pm <- list(list(empty, empty, empty, empty),
             list(full, NULL, full, NULL))
  est <- quiet(pulsar:::stars.stability(pm, .1, 4, 3))
  expect_equal(est$n.failed, c(0, 2))
  expect_equal(est$n.fit, c(4, 2))
  expect_equal(as.matrix(est$merge[[2]]), as.matrix(full))
  expect_equal(est$summary, c(0, 0))
  ## counting failures as empty graphs gives maximal instability instead
  pm0 <- list(pm[[1]], list(full, empty, full, empty))
  expect_equal(quiet(pulsar:::stars.stability(pm0, .1, 4, 3))$summary, c(0, 1))
})

test_that("a lambda where every fit failed is never selected", {
  full  <- Matrix::Matrix(1 - diag(3), sparse=TRUE)
  empty <- full * 0
  pm <- list(list(empty, empty), list(full, full), list(NULL, NULL))
  expect_warning(est <- pulsar:::stars.stability(pm, .1, 2, 3), "All subsample fits failed")
  expect_equal(est$n.failed, c(0, 0, 2))
  expect_true(est$opt.index != 3)
  expect_error(pulsar:::stars.stability(list(list(NULL)), .1, 1, 3), "every lambda")
})

test_that(".starsind moves off unselectable lambdas", {
  s <- c(0, .01, .02, .05, .2, .3)
  expect_equal(pulsar:::.starsind(s, .1), 4)
  expect_equal(pulsar:::.starsind(s, .1, selectable=rep(TRUE, 6)), 4)
  expect_warning(i <- pulsar:::.starsind(s, .1, selectable=c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)),
                 "unselectable")
  expect_equal(i, 2)
  expect_warning(i <- pulsar:::.starsind(s, .1, selectable=c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)),
                 "unselectable")
  expect_equal(i, 5)
})

test_that("reported failures are dropped and counted in pulsar", {
  out.n <- quiet(pulsar(X, fun_null,   fargs=list(lambda=lams), rep.num=10, seed=1, refit=FALSE))
  out.f <- quiet(pulsar(X, fun_flag,   fargs=list(lambda=lams), rep.num=10, seed=1, refit=FALSE))
  out.s <- quiet(pulsar(X, fun_silent, fargs=list(lambda=lams), rep.num=10, seed=1, refit=FALSE))
  out.0 <- quiet(pulsar(X, fun_null,   fargs=list(lambda=lams, nfail=0), rep.num=10, seed=1, refit=FALSE))

  nf <- out.n$stars$n.failed
  expect_equal(length(nf), length(lams))
  expect_true(all(nf[1:12] == 0))
  expect_true(all(nf[13:15] > 0) && all(nf[13:15] < 10))
  expect_equal(nf + out.n$stars$n.fit, rep(10, length(lams)))
  ## both reporting styles give the same result
  expect_equal(out.n$stars$summary, out.f$stars$summary)
  expect_equal(out.n$stars$n.failed, out.f$stars$n.failed)
  ## unaffected lambdas match a run with no failures
  expect_equal(out.n$stars$merge[1:12], out.0$stars$merge[1:12])
  ## frequencies at failed lambdas are over successful fits
  for (i in 13:15) {
    k <- out.n$stars$merge[[i]] * out.n$stars$n.fit[i]
    expect_equal(as.numeric(k), round(as.numeric(k)))
    expect_true(max(out.n$stars$merge[[i]]) <= 1)
  }
  ## silent failures keep the old behaviour
  expect_true(all(out.s$stars$n.failed == 0))
  expect_false(isTRUE(all.equal(out.s$stars$merge[13:15], out.n$stars$merge[13:15])))
})

test_that("pulsar warns about failures next to the selected lambda", {
  out <- quiet(pulsar(X, fun_null, fargs=list(lambda=lams, nfail=0), rep.num=10, seed=1, refit=FALSE))
  k   <- out$stars$opt.index
  nfail <- length(lams) - k + 1
  expect_warning(out2 <- pulsar(X, fun_null, fargs=list(lambda=lams, nfail=nfail),
                                rep.num=10, seed=1, refit=FALSE),
                 "Subsample fits failed near")
  ## unselectable: chosen lambda must have no failures
  out3 <- quiet(pulsar(X, fun_null, fargs=list(lambda=lams, nfail=nfail), rep.num=10,
                       seed=1, refit=FALSE, fail.action="unselectable"))
  expect_equal(out3$stars$n.failed[out3$stars$opt.index], 0)
  expect_lt(out3$stars$opt.index, k)
  ## a loose threshold keeps those lambdas selectable
  out4 <- quiet(pulsar(X, fun_null, fargs=list(lambda=lams, nfail=nfail), rep.num=10,
                       seed=1, refit=FALSE, fail.action="unselectable", fail.thresh=.99))
  expect_equal(out4$stars$opt.index, out2$stars$opt.index)
  expect_error(pulsar(X, fun_null, fargs=list(lambda=lams), fail.thresh=1), "fail.thresh")
  expect_error(pulsar(X, fun_null, fargs=list(lambda=lams), fail.action="foo"))
})

test_that("failures are handled with bounded stars and other criteria", {
  out <- quiet(pulsar(X, fun_null, fargs=list(lambda=lams), rep.num=10, seed=1,
                      refit=FALSE, lb.stars=TRUE, ub.stars=TRUE,
                      criterion=c("stars", "diss", "nc", "sufficiency")))
  expect_equal(length(out$stars$n.failed), length(lams))
  expect_equal(length(out$stars$merge), length(lams))
  expect_equal(out$stars$n.failed + out$stars$n.fit,
               c(rep(2, out$stars$ub.index - 1),
                 rep(10, out$stars$lb.index - out$stars$ub.index + 1),
                 rep(2, length(lams) - out$stars$lb.index)))
  expect_false(any(is.na(out$diss$summary)))
  expect_false(any(is.na(out$nc$summary)))
})

test_that("graphlet distances are computed over successful fits", {
  gp <- list(list(c(1, 0), c(0, 1), NULL), list(NULL, NULL, c(1, 1)))
  gm <- pulsar:::.gcdmerge(gp)
  expect_equal(gm$n.failed, c(1, 2))
  expect_equal(as.numeric(gm$merge[[1]]), sqrt(2))
  expect_equal(length(gm$merge[[2]]), 0)
  est <- pulsar:::gcd.stability(NULL, NULL, 3, 2, 2, merge=gm$merge, n.failed=gm$n.failed)
  expect_equal(est$summary, c(sqrt(2), NA))
  expect_equal(est$n.failed, c(1, 2))
})

test_that("batch.pulsar drops reported failures like pulsar", {
  skip_if_not_installed("batchtools")
  options(batchtools.verbose=FALSE, batchtools.progress=FALSE)
  wk <- tempdir()
  out   <- quiet(pulsar(X, fun_null, fargs=list(lambda=lams), rep.num=6, seed=1, refit=FALSE))
  outbt <- quiet(batch.pulsar(X, fun_null, fargs=list(lambda=lams), rep.num=6, seed=1,
                              refit=FALSE, wkdir=wk, cleanup=TRUE))
  expect_equal(outbt$stars$n.failed, out$stars$n.failed)
  expect_equal(outbt$stars$n.fit, out$stars$n.fit)
  expect_equal(outbt$stars$summary, out$stars$summary)
  expect_equal(outbt$stars$opt.index, out$stars$opt.index)
  expect_true(any(outbt$stars$n.failed > 0))

  outb  <- quiet(pulsar(X, fun_null, fargs=list(lambda=lams), rep.num=6, seed=1, refit=FALSE,
                        lb.stars=TRUE, ub.stars=TRUE))
  outbb <- quiet(batch.pulsar(X, fun_null, fargs=list(lambda=lams), rep.num=6, seed=1,
                              refit=FALSE, wkdir=wk, cleanup=TRUE, lb.stars=TRUE, ub.stars=TRUE))
  expect_equal(outbb$stars$n.failed, outb$stars$n.failed)
  expect_equal(outbb$stars$summary, outb$stars$summary)
  expect_equal(outbb$stars$opt.index, outb$stars$opt.index)
})
