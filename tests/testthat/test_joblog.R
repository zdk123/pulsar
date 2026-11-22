context("the record of a run's jobs")

p <- 8
set.seed(3)
dat  <- huge::huge.generator(60, p, "hub", verbose=FALSE)
lams <- getLamPath(getMaxCov(dat$data), .01, 6)
n <- nrow(dat$data)
## reports a failed fit at the last lambda on the full data
hf <- function(data, lambda) {
  out <- huge::huge(data, lambda=lambda, verbose=FALSE)
  list(path=out$path, failed=nrow(data) == n & lambda == min(lambda))
}
go <- function(..., fun=hf, rep.num=4)
  suppressWarnings(pulsar(dat$data, fun=fun, fargs=list(lambda=lams), rep.num=rep.num,
                          seed=1, ...))
cols <- c("id", "kind", "stage", "chunk", "n.lambda", "lambda.max", "lambda.min",
          "status", "start", "secs", "host", "failed.fits", "error")

test_that("the result records each job", {
  out <- go()
  j <- getJobs(out)
  expect_identical(names(j), cols)
  expect_equal(j$id, 0:4)
  expect_equal(j$kind, c("full", rep("subsample", 4)))
  expect_true(all(j$status == "done"))
  expect_true(all(j$secs >= 0) && inherits(j$start, "POSIXct"))
  expect_equal(j$host, rep(Sys.info()[["nodename"]], 5))
  expect_equal(j$failed.fits, c(1, 0, 0, 0, 0))
  expect_true(all(j$n.lambda == 6 & j$lambda.max == max(lams)))
  ## bounded StARS: stage 2 jobs fit the narrowed path
  jb <- getJobs(go(lb.stars=TRUE, ub.stars=TRUE))
  expect_equal(jb$stage, c(1, 1, 1, 2, 2))
  expect_true(all(jb$n.lambda[jb$stage == 2] <= 6))
  ## early stopping: a row per chunk
  je <- getJobs(go(stop.chunk=2))
  expect_true(all(!is.na(je$chunk[je$kind == "subsample"])))
  expect_true(is.na(je$chunk[je$kind == "full"]))
  expect_null(getJobs(structure(list(), class="batch.pulsar")))
})

test_that("run.dir keeps the record in jobs.csv, and resumed jobs keep their times", {
  rd <- tempfile()
  a <- go(run.dir=rd)
  csv <- utils::read.csv(file.path(rd, "jobs.csv"), stringsAsFactors=FALSE)
  expect_identical(names(csv), cols)
  expect_equal(sort(csv$id), 0:4)
  expect_true(all(csv$status == "done"))
  expect_match(csv$start, "^\\d{4}-\\d\\d-\\d\\dT[0-9:.]+Z$")
  unlink(file.path(rd, "job_000002.rds"))
  b <- suppressMessages(go(run.dir=rd))
  jb <- getJobs(b)
  expect_equal(jb$status, c("loaded", "loaded", "done", "loaded", "loaded"))
  expect_equal(jb$secs[-3], getJobs(a)$secs[-3])
  ## each call appends its rows: the loaded ones too
  csv <- utils::read.csv(file.path(rd, "jobs.csv"), stringsAsFactors=FALSE)
  expect_equal(csv$status[6:10][order(csv$id[6:10])], jb$status)
})

test_that("a job whose error stops a serial run is in jobs.csv", {
  rd <- tempfile()
  bad <- function(data, lambda) {
    if (nrow(data) < n) stop("solver blew up, badly")
    huge::huge(data, lambda=lambda, verbose=FALSE)
  }
  expect_error(go(fun=bad, run.dir=rd), "solver blew up")
  csv <- utils::read.csv(file.path(rd, "jobs.csv"), stringsAsFactors=FALSE)
  expect_equal(csv$status, c("done", "error"))
  expect_equal(csv$error[2], "solver blew up, badly")
})

test_that("failed jobs are recorded with their error", {
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  rd <- tempfile()
  bad <- function(data, lambda) {
    if (nrow(data) < n && data[1, 1] > 0) stop("solver blew up")
    huge::huge(data, lambda=lambda, verbose=FALSE)
  }
  out <- go(fun=bad, backend="future", run.dir=rd, rep.num=8)
  j <- getJobs(out)
  expect_true(any(j$status == "error"))
  expect_true(all(grepl("solver blew up", j$error[j$status == "error"])))
  expect_true(all(is.na(j$secs[j$status == "error"])))
  csv <- utils::read.csv(file.path(rd, "jobs.csv"), stringsAsFactors=FALSE)
  expect_equal(sum(csv$status == "error"), sum(j$status == "error"))
})

test_that("progress is one bar per run, in lambdas fitted, with grouped messages", {
  skip_if_not_installed("progressr")
  ## progressr is off in non-interactive sessions unless enabled
  oopt <- options(progressr.enable=TRUE)
  on.exit(options(oopt), add=TRUE)
  watch <- function(expr) {
    rec <- new.env() ; rec$bars <- 0L ; rec$amount <- 0 ; rec$msgs <- character(0)
    h <- progressr::make_progression_handler("watch", reporter=list(
           initiate=function(config, state, progression, ...) {
             rec$bars <- rec$bars + 1L ; rec$steps <- config$max_steps },
           update=function(config, state, progression, ...) {
             rec$amount <- rec$amount + progression$amount
             if (length(progression$message)) rec$msgs <- c(rec$msgs, progression$message) }))
    progressr::with_progress(expr, handlers=h)
    as.list(rec)
  }
  w <- watch(go())
  expect_equal(w$bars, 1L)
  expect_equal(w$steps, 6 * 5)
  expect_equal(w$amount, 6 * 5)
  expect_true(all(c("full-data fit done", "subsample 4 done") %in% w$msgs))
  ## bounded StARS: the stages are named, and the skipped fits are credited
  w <- watch(go(lb.stars=TRUE, ub.stars=TRUE))
  expect_equal(c(w$bars, w$amount), c(1, 6 * 5))
  expect_true(any(grepl("^stage 1 \\(bounds, lambdas 1-6 of 6\\): subsample 2 done$", w$msgs)))
  expect_true(any(grepl("^stage 2 \\(lambdas? [0-9-]+ of 6\\): subsample 3 done$", w$msgs)))
  expect_true(any(grepl("^bounds: .*skipping", w$msgs)))
  ## early stopping: chunks are named, and the bar ends where the run stopped
  w <- watch(go(stop.chunk=2, thresh=.2))
  expect_equal(c(w$bars, w$amount), c(1, 6 * 5))
  expect_true(any(grepl("^chunk 1 \\(lambdas 1-2 of 6\\): subsample 1 done$", w$msgs)))
  expect_true(tail(w$msgs, 1) == "done" || grepl("^stopped at lambda", tail(w$msgs, 1)))
  ## resumed jobs and future backends advance the bar too
  rd <- tempfile() ; go(run.dir=rd)
  w <- suppressMessages(watch(go(run.dir=rd)))
  expect_true("subsample 1 loaded" %in% w$msgs)
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan), add=TRUE)
  w <- watch(go(backend="future"))
  expect_equal(c(w$bars, w$amount), c(1, 6 * 5))
})
