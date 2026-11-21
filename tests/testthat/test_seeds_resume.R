context("job seeds and resuming runs")

p <- 8
set.seed(1)
dat  <- huge::huge.generator(60, p, "hub", verbose=FALSE)
lams <- getLamPath(getMaxCov(dat$data), .01, 6)

## a fun that draws random numbers: it jitters the data before fitting
jitter_fun <- function(data, lambda) {
  data <- data + matrix(stats::rnorm(length(data), sd=.5), nrow(data))
  huge::huge(data, lambda=lambda, verbose=FALSE)
}
## call pulsar inside a function, so `fun` is not a top-level global (#43)
## (the jittered fits can put the optimum at the end of this short path)
run <- function(..., fun=jitter_fun, rep.num=6) {
  f <- fun
  withCallingHandlers(
    pulsar(dat$data, fun=f, fargs=list(lambda=lams), rep.num=rep.num, seed=10, ...),
    warning=function(w) if (grepl("Optimal lambda may be|Accurate lower bound", conditionMessage(w)))
                          invokeRestart("muffleWarning"))
}
strip <- function(x) { x$call <- NULL ; x$envir <- NULL ; x }
## PSOCK workers load the installed pulsar, which under devtools::test may
## be an older version
dev_loaded <- function() requireNamespace("pkgload", quietly=TRUE) &&
                         pkgload::is_dev_package("pulsar")

test_that("job streams depend only on job.seed and the job id", {
  s5  <- pulsar:::.jobstreams(42, 5L)
  s10 <- pulsar:::.jobstreams(42, 10L)
  expect_length(s5, 6L)
  expect_identical(s10[1:6], s5)
  expect_false(identical(s5[[1]], s5[[2]]))
  expect_false(identical(pulsar:::.jobstreams(43, 5L), s5))
  expect_error(pulsar:::.jobstreams("a", 1L), "job.seed")
})

test_that("planning and running seeded jobs leave the caller's RNG alone", {
  RNGkind("Mersenne-Twister")
  set.seed(5) ; before <- stats::runif(3)
  set.seed(5)
  jobs <- pulsar:::.planjobs(20, 4, .8, TRUE, FALSE, seed=NULL, job.seed=1)
  X <- matrix(stats::rnorm(40), 20)
  set.seed(5)
  pulsar:::.run_job(jobs[[2]], X, function(data, lambda) list(path=list(stats::runif(1))), list(lambda=1))
  expect_equal(stats::runif(3), before)
  expect_equal(RNGkind()[1], "Mersenne-Twister")
})

test_that("a seeded job gives the same result whatever RNG state it starts in", {
  spec <- pulsar:::.planjobs(20, 3, .8, FALSE, FALSE, seed=1, job.seed=9)[[2]]
  X <- matrix(stats::rnorm(40), 20)
  f <- function(data, lambda) list(path=list(stats::rnorm(2)))
  set.seed(1) ; a <- pulsar:::.run_job(spec, X, f, list(lambda=1))
  set.seed(2) ; b <- pulsar:::.run_job(spec, X, f, list(lambda=1))
  expect_identical(a, b)
})

test_that("job.seed makes a stochastic fun reproducible", {
  for (lb in c(FALSE, TRUE)) {
    expect_equal(strip(run(job.seed=3, lb.stars=lb)), strip(run(job.seed=3, lb.stars=lb)))
  }
  ## without job.seed the jitter differs between runs
  expect_false(isTRUE(all.equal(strip(run()), strip(run()))))
  expect_false(isTRUE(all.equal(strip(run(job.seed=3)), strip(run(job.seed=4)))))
})

test_that("job.seed gives the same result on every backend", {
  ref <- strip(run(job.seed=3, lb.stars=TRUE))
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  expect_equal(strip(run(job.seed=3, lb.stars=TRUE, backend="future")), ref)
  skip_on_cran()
  skip_if(dev_loaded(), "PSOCK workers need the installed pulsar")
  expect_equal(strip(run(job.seed=3, lb.stars=TRUE, ncores=2)), ref)
  future::plan(future::multisession, workers=2)
  expect_equal(strip(run(job.seed=3, lb.stars=TRUE, backend="future")), ref)
})

test_that("refit() of a refit=FALSE run uses the full-data job's stream", {
  a <- run(job.seed=3, refit=TRUE)
  b <- run(job.seed=3, refit=FALSE)
  expect_equal(refit(b)$est, a$est)
})

## ---- resuming ------------------------------------------------------------

## `fun` that fails when the file `stopfile` exists, as a session dying
## mid-run would stop the run, and counts its calls in `countfile`
resumable <- function(stopfile, countfile) {
  function(data, lambda) {
    cat("x", file=countfile, append=TRUE)
    if (file.exists(stopfile) && nchar(readLines(countfile, warn=FALSE)) > 3)
      stop("session died")
    huge::huge(data, lambda=lambda, verbose=FALSE)
  }
}
calls <- function(countfile) {
  if (!file.exists(countfile)) return(0L)
  n <- nchar(readLines(countfile, warn=FALSE)) ; unlink(countfile) ; n
}

test_that("a run stopped midway resumes from run.dir, running only the rest", {
  for (lb in c(FALSE, TRUE)) {
    info <- sprintf("lb.stars=%s", lb)
    stopfile <- tempfile() ; countfile <- tempfile() ; rd <- tempfile()
    f <- resumable(stopfile, countfile)
    ref <- strip(run(fun=f, lb.stars=lb))
    expect_equal(calls(countfile), 7L, info=info)

    ## serial: the 4th call stops the whole run; 3 jobs are saved by then
    file.create(stopfile)
    expect_error(run(fun=f, lb.stars=lb, run.dir=rd), "session died", info=info)
    calls(countfile)
    expect_length(list.files(rd, "^job_"), 3L)

    unlink(stopfile)
    expect_message(out <- run(fun=f, lb.stars=lb, run.dir=rd), "3 of 7", info=info)
    expect_equal(calls(countfile), 4L, info=info)
    expect_equal(strip(out), ref, info=info)

    ## everything saved: nothing runs again
    expect_message(again <- run(fun=f, lb.stars=lb, run.dir=rd), "7 of 7", info=info)
    expect_equal(calls(countfile), 0L, info=info)
    expect_equal(strip(again), ref, info=info)
  }
})

test_that("run.dir keeps the plan, so a resumed run needs no seed", {
  rd <- tempfile()
  hf <- function(data, lambda) huge::huge(data, lambda=lambda, verbose=FALSE)
  go <- function() pulsar(dat$data, fun=hf, fargs=list(lambda=lams), rep.num=4,
                          run.dir=rd)
  a <- go()
  unlink(file.path(rd, c("job_000002.rds", "job_000003.rds")))
  b <- suppressMessages(go())
  expect_equal(strip(b), strip(a))
})

test_that("resuming works across backends", {
  skip_if_not_installed("future")
  rd <- tempfile()
  ref <- strip(run(fun=jitter_fun, job.seed=2, run.dir=rd))
  unlink(file.path(rd, "job_000004.rds"))
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  expect_message(out <- run(fun=jitter_fun, job.seed=2, run.dir=rd, backend="future"),
                 "6 of 7")
  expect_equal(strip(out), ref)
})

test_that("run.dir refuses a different run", {
  rd <- tempfile()
  run(fun=jitter_fun, job.seed=1, run.dir=rd)
  expect_error(run(fun=jitter_fun, job.seed=2, run.dir=rd), "different run")
  expect_error(run(fun=jitter_fun, job.seed=1, rep.num=5, run.dir=rd), "different run")
  hf <- function(data, lambda) huge::huge(data, lambda=lambda, verbose=FALSE)
  expect_error(run(fun=hf, job.seed=1, run.dir=rd), "different run")
})

test_that("the run key sees captured values, not bytecode or later variables", {
  key <- function(fun, ...) pulsar:::.runkey(dat$data, fun, list(lambda=lams), ...)
  mk <- function(a) function(data, lambda) huge::huge(data * a, lambda=lambda, verbose=FALSE)
  expect_false(identical(key(mk(1)), key(mk(2))))
  f <- mk(1)
  k1 <- key(f)
  assign("unrelated", 1:10, envir=environment(f))
  expect_identical(key(f), k1)
  g <- function(data, lambda) huge::huge(data, lambda=lambda, verbose=FALSE)
  expect_identical(key(compiler::cmpfun(g)), key(g))
  expect_false(identical(pulsar:::.runkey(dat$data[-1, ], g, list(lambda=lams)), key(g)))
  expect_false(identical(pulsar:::.runkey(dat$data, g, list(lambda=lams[-1])), key(g)))
})

test_that("job streams don't depend on the caller's RNG kinds", {
  restore <- pulsar:::.saverng()
  on.exit(restore())
  ref <- pulsar:::.jobstreams(3, 2L)
  suppressWarnings(RNGkind("Mersenne-Twister", "Box-Muller", "Rounding"))
  expect_identical(pulsar:::.jobstreams(3, 2L), ref)
  spec <- pulsar:::.jobspec(1, "subsample", 1:3, seed=ref[[2]])
  f <- function(data, lambda) list(path=list(stats::rnorm(1)))
  expect_silent(pulsar:::.run_job(spec, matrix(0, 3, 2), f, list(lambda=1)))
  expect_equal(RNGkind(), c("Mersenne-Twister", "Box-Muller", "Rounding"))
})

test_that("a result that can't be saved stops the run", {
  rd <- tempfile()
  run(fun=jitter_fun, job.seed=1, run.dir=rd)
  unlink(file.path(rd, "job_000003.rds"))
  Sys.chmod(rd, "0555")
  on.exit(Sys.chmod(rd, "0755"))
  skip_if(file.access(rd, 2) == 0, "cannot make run.dir read-only here")
  expect_error(suppressMessages(run(fun=jitter_fun, job.seed=1, run.dir=rd)), "could not save")
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan), add=TRUE)
  expect_error(suppressMessages(run(fun=jitter_fun, job.seed=1, run.dir=rd, backend="future")),
               "could not save")
})

test_that("a saved stage-2 result for another lambda path is run again", {
  rd <- tempfile() ; countfile <- tempfile()
  f <- resumable(tempfile(), countfile)
  ref <- strip(run(fun=f, lb.stars=TRUE, run.dir=rd))
  calls(countfile)
  ## pretend job 5 was fitted on a different (stage 2) path
  path <- file.path(rd, "job_000005.rds")
  saved <- readRDS(path) ; saved$lambda <- saved$lambda[-1] ; saveRDS(saved, path)
  out <- suppressMessages(run(fun=f, lb.stars=TRUE, run.dir=rd))
  expect_equal(calls(countfile), 1L)
  expect_equal(strip(out), ref)
})
