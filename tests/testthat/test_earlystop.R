context("early stopping for StARS")

p <- 10
set.seed(2)
dat  <- huge::huge.generator(80, p, "hub", verbose=FALSE)
lams <- getLamPath(getMaxCov(dat$data), .005, 16)
n <- nrow(dat$data)

## a fun that fits each lambda on its own (cold), so a fit doesn't depend on
## the other lambdas; it counts the lambdas it fits in `countfile`, and
## reports failed fits: on the full data from lambda `fullfail` on, and on
## subsamples where fail(data, lambda index) is TRUE
coldfun <- function(countfile=tempfile(), fullfail=Inf, fail=function(data, i) FALSE,
                    boom=function(data, lambda) FALSE) {
  function(data, lambda) {
    cat(strrep("x", length(lambda)), file=countfile, append=TRUE)
    if (boom(data, lambda)) stop("boom")
    i0 <- match(lambda[1], lams)
    path <- lapply(seq_along(lambda), function(j)
              huge::huge(data, lambda=lambda[j], verbose=FALSE)$path[[1]])
    ii <- i0 + seq_along(lambda) - 1L
    failed <- if (nrow(data) == n) ii >= fullfail else vapply(ii, function(i) fail(data, i), NA)
    list(path=path, failed=failed)
  }
}
nfits <- function(countfile) {
  if (!file.exists(countfile)) return(0L)
  k <- nchar(readLines(countfile, warn=FALSE)) ; unlink(countfile) ; k
}
go <- function(fun, ..., rep.num=8, lambda=lams)
  suppressWarnings(pulsar(dat$data, fun=fun, fargs=list(lambda=lambda), rep.num=rep.num,
                          seed=7, thresh=.05, ...))
same_pick <- function(a, b, cut=b$stars$stop$index) {
  expect_identical(b$stars$opt.index, a$stars$opt.index)
  expect_identical(b$stars$summary[seq_len(cut)], a$stars$summary[seq_len(cut)])
  expect_identical(b$stars$n.failed[seq_len(cut)], a$stars$n.failed[seq_len(cut)])
}

test_that("stopping at the crossing gives the whole path's pick, fitting fewer lambdas", {
  cf <- tempfile()
  ref <- go(coldfun(cf))
  full <- nfits(cf)
  cross <- which(ref$stars$summary >= .05)[1]
  expect_true(cross < length(lams) - 3)     # so there's something to skip
  for (chunk in c(1, 2, 5)) for (margin in 0:2) {
    out <- go(coldfun(cf), stop.chunk=chunk, stop.margin=margin)
    info <- sprintf("chunk %d, margin %d", chunk, margin)
    same_pick(ref, out)
    cut <- out$stars$stop$index
    expect_identical(cut, as.integer(min(cross + margin, length(lams))), info=info)
    expect_identical(out$stars$stop$reason, "crossing", info=info)
    expect_identical(out$stars$stop$crossing, cross, info=info)
    expect_true(all(is.na(out$stars$summary[-seq_len(cut)])), info=info)
    expect_equal(out$stars$n.run, c(rep(8, cut), rep(0, length(lams) - cut)), info=info)
    expect_identical(out$est, ref$est, info=info)
    expect_equal(nfits(cf), length(lams) + 8 * cut, info=info)
  }
  expect_equal(full, length(lams) * 9)
})

test_that("stop.full fits subsamples only before the full-data fit's first failure", {
  cf <- tempfile()
  out <- go(coldfun(cf, fullfail=6), stop.full=TRUE)
  expect_equal(nfits(cf), length(lams) + 8 * 5)
  expect_identical(out$stars$stop[c("index", "reason", "full.failed")],
                   list(index=5L, reason="full.failed", full.failed=6L))
  ## the same as a run on the lambdas before the failure
  ref <- go(coldfun(), lambda=lams[1:5])
  same_pick(ref, out)
  expect_true(all(is.na(out$stars$summary[6:16])))
  ## with chunks too: the cut is the earlier of the two
  out2 <- go(coldfun(cf, fullfail=12), stop.full=TRUE, stop.chunk=2)
  ref2 <- go(coldfun(), lambda=lams[1:11])
  same_pick(ref2, out2)
  cross <- which(ref2$stars$summary >= .05)[1]
  expect_identical(out2$stars$stop$index, as.integer(min(cross + 1L, 11L)))
  expect_error(go(coldfun(fullfail=1), stop.full=TRUE), "first lambda")
  ## without chunks, all lambdas before the failure are fitted, past the crossing
  out3 <- go(coldfun(fullfail=15), stop.full=TRUE)
  expect_true(out3$stars$stop$crossing < 14)
  expect_identical(out3$stars$stop[c("index", "reason")], list(index=14L, reason="full.failed"))
  out4 <- go(coldfun(), stop.full=TRUE)
  expect_identical(out4$stars$stop[c("index", "reason")], list(index=16L, reason="end"))
  expect_identical(out4$stars$opt.index, go(coldfun())$stars$opt.index)
})

test_that("subsample failures and fail.action are handled as with the whole path", {
  ## failures at and around the crossing, on some subsamples
  fl <- function(data, i) i >= 3 && i <= 6 && (sum(data[1, ]) > 0)
  for (fa in c("drop", "unselectable")) {
    ref <- go(coldfun(fail=fl), fail.action=fa)
    out <- go(coldfun(fail=fl), fail.action=fa, stop.chunk=1, stop.margin=0)
    same_pick(ref, out)
  }
  ## no lambda selectable up to the crossing: keep fitting until one is
  fl2 <- function(data, i) i <= 8
  ref <- go(coldfun(fail=fl2), fail.action="unselectable")
  out <- go(coldfun(fail=fl2), fail.action="unselectable", stop.chunk=1, stop.margin=0)
  same_pick(ref, out)
  expect_true(out$stars$stop$index >= 9)
})

test_that("a subsample that errors in a later chunk is dropped at every lambda", {
  ## (a serial run stops at any error; on futures a failed job is dropped)
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  bm <- function(data, lambda) nrow(data) < n && data[1, 1] > 0 && lambda[1] < lams[3]
  ref <- go(coldfun(boom=function(data, lambda) nrow(data) < n && data[1, 1] > 0), backend="future")
  out <- go(coldfun(boom=bm), stop.chunk=2, stop.margin=0, backend="future")
  same_pick(ref, out)
  expect_true(all(out$stars$n.run[out$stars$n.run > 0] < 8))
})

test_that("early stopping works on future backends and resumes from run.dir", {
  ref <- go(coldfun(), stop.chunk=2)
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  out <- go(coldfun(), stop.chunk=2, backend="future")
  same_pick(ref, out)
  expect_identical(out$stars$stop, ref$stars$stop)
  skip_on_cran()
  skip_if(requireNamespace("pkgload", quietly=TRUE) && pkgload::is_dev_package("pulsar"),
          "multisession workers need the installed pulsar")
  future::plan(future::multisession, workers=2)
  out <- go(coldfun(), stop.chunk=2, backend="future")
  same_pick(ref, out)
})

test_that("a stopped run resumes from run.dir", {
  cf <- tempfile() ; rd <- tempfile()
  ref <- go(coldfun(cf), stop.chunk=2, job.seed=1)
  total <- nfits(cf)
  out <- go(coldfun(cf), stop.chunk=2, job.seed=1, run.dir=rd)
  nfits(cf)
  expect_true(length(list.files(rd, "_chunk_")) > 0)
  ## lose the last chunk of two subsamples: only those are refitted
  last <- max(as.integer(sub(".*_chunk_([0-9]+).rds", "\\1", list.files(rd, "_chunk_"))))
  unlink(file.path(rd, sprintf("job_%06d_chunk_%03d.rds", 1:2, last)))
  expect_message(again <- go(coldfun(cf), stop.chunk=2, job.seed=1, run.dir=rd), "chunk results")
  expect_equal(nfits(cf), 2 * length(seq(2 * last - 1, ref$stars$stop$index)))
  same_pick(ref, again)
  expect_identical(again$stars$stop, ref$stars$stop)
})

test_that("early stopping checks its arguments", {
  f <- coldfun()
  expect_error(go(f, stop.chunk=0), "stop.chunk")
  expect_error(go(f, stop.chunk=1.5), "stop.chunk")
  expect_error(go(f, stop.chunk=2, stop.margin=-1), "stop.margin")
  expect_error(go(f, stop.full=NA), "stop.full")
  expect_error(go(f, stop.chunk=2, criterion=c("stars", "nc")), "criterion")
  expect_error(go(f, stop.chunk=2, lb.stars=TRUE), "lb.stars")
  expect_error(go(f, stop.full=TRUE, refit=FALSE), "refit")
})

test_that("without early stopping nothing changes", {
  ref <- go(coldfun())
  expect_null(ref$stars$stop)
  expect_null(ref$stars$n.run)
  expect_false(anyNA(ref$stars$summary))
})
