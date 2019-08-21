context("error reporting")


p <- 4
rseed <- 10010
set.seed(rseed)
dat <- huge::huge.generator(p*5, p, "hub", verbose=FALSE, v=.4, u=.2)
set.seed(rseed)
X <- MASS::mvrnorm(p*5, mu=rep(0,p), Sigma=dat$sigma, empirical=TRUE)
X[4,] <- NA ## sample 4 not in first two subsamples, given rseed

# TODO: test only a subset of jobs failing
huge_error <- function(data, lambda, f=stop, ...) {
  if (any(is.na(data))) {
    f('NA detected in data', call.=FALSE)
  }
  return(huge::huge(na.exclude(data), lambda=lambda, ...))
}

filter_warning <- function(expr, regxpr) {
  withCallingHandlers(
    force(expr),
  warning=function(w) {
     if (grepl(regxpr, w$message))
        invokeRestart("muffleWarning")
  } )
}

hargs <- list(lambda=getLamPath(1, .005, 10), verbose=FALSE, f=stop)
test_that("Job stops when every parallel job has errors", {
  expect_error(
    pulsar(X, fun=huge_error, hargs, ncores=1, subsample.ratio=1, seed=rseed,
          rep.num=8, refit=FALSE),
    regexp = "NA")

  expect_error(
    pulsar(X, fun=huge_error, hargs, ncores=1, subsample.ratio=1, seed=rseed,
          rep.num=8, refit=FALSE, lb.stars=TRUE),
    regexp = "NA")

  skip_on_os("windows")
  expect_error(
    pulsar(X, fun=huge_error, hargs, ncores=2, subsample.ratio=1, seed=rseed,
          rep.num=8, refit=FALSE),
    regexp = "NA")

  skip_on_os("windows")
  expect_error(
    pulsar(X, fun=huge_error, hargs, ncores=2, subsample.ratio=1, seed=rseed,
          rep.num=8, refit=FALSE, lb.stars=TRUE),
    regexp = "NA")

})


test_that("Job continues when subset of parallel jobs have errors/warnings", {
  hargs$f <- stop
  expect_error(
    pulsar(X, fun=huge_error, hargs, ncores=1, subsample.ratio=.5, seed=rseed,
          rep.num=8, refit=FALSE),
    regexp = "NA")

  expect_error(
    pulsar(X, fun=huge_error, hargs, ncores=1, subsample.ratio=.5, seed=rseed,
          rep.num=8, refit=FALSE, lb.stars=TRUE),
    regexp = "NA")

  skip_on_os("windows")
  expect_warning(
    filter_warning(
      out1 <- pulsar(X, fun=huge_error, hargs, ncores=2, subsample.ratio=.5,
            seed=rseed, rep.num=8, refit=FALSE), "Optimal lambda"),
    regexp = "NA")

  skip_on_os("windows")
  expect_warning(
    filter_warning(
      out2 <- pulsar(X, fun=huge_error, hargs, ncores=2, subsample.ratio=.5,
            seed=rseed, rep.num=8, refit=FALSE, lb.stars=TRUE),
      "Optimal lambda"),
    regexp = "NA")

  expect_equivalent(out1$stars$opt.ind, out2$stars$opt.ind)

  hargs$f <- warning
  expect_warning(
    filter_warning(
      out1 <- pulsar(X, fun=huge_error, hargs, ncores=1, subsample.ratio=.5,
             seed=rseed, rep.num=8, refit=FALSE), "Optimal lambda"),
    regexp = "NA")

  expect_warning(
    filter_warning(
      out2 <- pulsar(X, fun=huge_error, hargs, ncores=1, subsample.ratio=.5,
             seed=rseed, rep.num=8, refit=FALSE, lb.stars=TRUE),
      "Optimal lambda"),
    regexp = "NA")

  expect_equivalent(out1$stars$opt.ind, out2$stars$opt.ind)

})
