context("the raw StARS curve (summary.raw)")

p <- 10
set.seed(6)
dat  <- huge::huge.generator(80, p, "hub", verbose=FALSE)
lams <- getLamPath(getMaxCov(dat$data), .005, 12)
n <- nrow(dat$data)

## fits each lambda on its own; reports subsample fits at lambda `allfail`
## as failed (every subsample), and some at lambda `somefail`
coldfun <- function(allfail=0, somefail=0) {
  function(data, lambda) {
    ii <- match(lambda, lams)
    path <- lapply(lambda, function(l) huge::huge(data, lambda=l, verbose=FALSE)$path[[1]])
    sub <- nrow(data) < n
    list(path=path, failed=sub & (ii == allfail | (ii == somefail & data[1, 1] > 0)))
  }
}
go <- function(fun=coldfun(), ..., rep.num=8)
  suppressWarnings(pulsar(dat$data, fun=fun, fargs=list(lambda=lams), rep.num=rep.num,
                          seed=2, thresh=.05, ...))
## the statistic from merge: 4 theta (1 - theta) over pairs
frommerge <- function(st) vapply(st$merge, function(m) 4 * sum(m * (1 - m)) / (p * (p - 1)),
                                 numeric(1))

test_that("summary is the running maximum of summary.raw, which matches merge", {
  for (f in list(coldfun(), coldfun(somefail=4))) {
    st <- go(f)$stars
    expect_length(st$summary.raw, length(lams))
    expect_equal(st$summary, cummax(st$summary.raw))
    expect_equal(st$summary.raw, frommerge(st))
  }
  ## the curve should not be monotone already, or this proves little
  expect_false(isTRUE(all.equal(st$summary.raw, cummax(st$summary.raw))))
})

test_that("summary.raw is NA where every subsample fit failed", {
  st <- go(coldfun(allfail=3))$stars
  expect_true(is.na(st$summary.raw[3]))
  expect_false(anyNA(st$summary.raw[-3]))
  raw0 <- st$summary.raw ; raw0[is.na(raw0)] <- 0
  expect_equal(st$summary, cummax(raw0))
})

test_that("summary.raw is the same on futures, NA past an early stop, spliced with lb.stars", {
  ref <- go()$stars
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  expect_identical(go(backend="future")$stars$summary.raw, ref$summary.raw)

  es <- go(stop.chunk=2, stop.margin=0)$stars
  cut <- es$stop$index
  expect_true(cut < length(lams))
  expect_identical(es$summary.raw[seq_len(cut)], ref$summary.raw[seq_len(cut)])
  expect_true(all(is.na(es$summary.raw[-seq_len(cut)])))
  expect_equal(es$summary[seq_len(cut)], cummax(es$summary.raw[seq_len(cut)]))

  ## lb.stars: between the bounds, from all subsamples, monotonized from the
  ## upper bound; outside, from the first 2 (as summary is)
  lb <- go(lb.stars=TRUE, ub.stars=TRUE)$stars
  pind <- lb$ub.index:lb$lb.index
  expect_equal(lb$summary[pind], cummax(lb$summary.raw[pind]))
  expect_equal(lb$summary.raw, frommerge(lb))
})

test_that("batch.pulsar has summary.raw too", {
  skip_if_not_installed("batchtools")
  suppressPackageStartupMessages(library(batchtools))
  options(batchtools.verbose=FALSE, batchtools.progress=FALSE)
  wd <- fs::path_real(tempdir())
  for (lb in c(FALSE, TRUE)) {
    b <- suppressWarnings(batch.pulsar(dat$data, coldfun(), list(lambda=lams), rep.num=4,
                                       seed=2, thresh=.05, lb.stars=lb, ub.stars=lb,
                                       wkdir=wd, refit=FALSE, cleanup=TRUE))$stars
    expect_length(b$summary.raw, length(lams))
    expect_equal(b$summary.raw, frommerge(b), info=lb)
    pind <- if (lb) b$ub.index:b$lb.index else seq_along(lams)
    expect_equal(b$summary[pind], cummax(b$summary.raw[pind]), info=lb)
  }
})
