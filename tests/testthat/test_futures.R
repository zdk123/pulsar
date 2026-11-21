context("future backends")

skip_if_not_installed("future")

p <- 8
set.seed(1)
dat  <- huge::huge.generator(60, p, "hub", verbose=FALSE)
lams <- getLamPath(getMaxCov(dat$data), .01, 6)

## call pulsar inside a function, so `fun` is not a top-level global (#43)
run <- function(...) {
  hfun <- function(data, lambda) huge::huge(data, lambda=lambda, verbose=FALSE)
  pulsar(dat$data, fun=hfun, fargs=list(lambda=lams), rep.num=6, seed=10, ...)
}
strip <- function(x) { x$call <- NULL ; x$envir <- NULL ; x }
combos <- expand.grid(refit=c(FALSE, TRUE), lb=c(FALSE, TRUE))

test_that("backend argument is checked", {
  expect_error(run(backend="slurm"), "backend")
  expect_error(run(backend=c("local", "future")), "backend")
})

test_that("future backend gives the same result as the local backend", {
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  for (k in seq_len(nrow(combos))) {
    refit <- combos$refit[k] ; lb <- combos$lb[k]
    expect_equal(strip(run(backend="future", refit=refit, lb.stars=lb)),
                 strip(run(refit=refit, lb.stars=lb)),
                 info=sprintf("refit=%s lb.stars=%s", refit, lb))
  }
})

test_that("future backend drops a job whose fun errors, with a warning", {
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  X <- dat$data ; X[3, ] <- NA
  err <- function(data, lambda) {
    if (anyNA(data)) stop("NA in data")
    huge::huge(data, lambda=lambda, verbose=FALSE)
  }
  ## row 3 is in some subsamples but not all, given this seed
  jobs <- pulsar:::.planjobs(nrow(X), 6, .8, FALSE, FALSE, 10)
  nbad <- sum(vapply(jobs, function(j) 3L %in% j$rows, logical(1)))
  skip_if(nbad == 0 || nbad == 6)
  expect_warning(out <- pulsar(X, fun=err, fargs=list(lambda=lams), rep.num=6,
                               seed=10, refit=FALSE, backend="future"),
                 "NA in data")
  expect_true(all(out$stars$n.fit == 6 - nbad))
})

## ---- Modal backend, with a local stand-in for the Python helper ----------

## The stand-in speaks the helper's protocol: it replies to a submission
## first, then evaluates the future, so R sees it running until its result
## file appears. Submission number `fail` (if any) is reported as a lost
## container; at submission `die` the helper exits without a result. Under
## devtools::test it loads this source tree; under R CMD check, the installed
## package.
fake_helper <- function(fail=0L, die=0L) {
  pkgdir <- if (requireNamespace("pkgload", quietly=TRUE) &&
                pkgload::is_dev_package("pulsar"))
              pkgload::pkg_path(testthat::test_path()) else ""
  script <- tempfile(fileext=".R")
  died <- tempfile()  # only the first helper dies
  writeLines(c(
    sprintf('pkgdir <- %s; fail <- %d; die <- %d; died <- %s; n <- 0L', deparse(pkgdir),
            as.integer(fail), as.integer(die), deparse(died)),
    'if (nzchar(pkgdir)) pkgload::load_all(pkgdir, quiet=TRUE) else library(pulsar)',
    'con <- file("stdin", "r")',
    'say <- function(...) { cat(..., "\\n", sep=""); flush(stdout()) }',
    'say("ready")',
    'while (length(line <- readLines(con, n=1L)) > 0) {',
    '  req <- strsplit(line, "\\t", fixed=TRUE)[[1]]',
    '  if (req[1] == "quit") { say("ok"); break }',
    '  if (req[1] == "cancel") { say("ok"); next }',
    '  n <- n + 1L ; out <- req[4]',
    '  inp <- tempfile() ; file.copy(req[3], inp)  # R deletes req[3] once answered',
    '  say("ok\\tfc-", n)',
    '  if (n == die && !file.exists(died)) { file.create(died); quit(status=1) }',
    '  if (n == fail) writeLines("container lost", paste0(out, ".err"))',
    '  else {',
    '    tmp <- paste0(out, ".part")',
    '    pulsar:::.modal_eval(inp, tmp)',
    '    file.rename(tmp, out)',
    '  }',
    '}'), script)
  c(file.path(R.home("bin"), "Rscript"), script)
}

with_modal <- function(code, fail=0L, die=0L, ...) {
  oopt  <- options(pulsar.modal.helper=fake_helper(fail, die), pulsar.modal.start.timeout=60)
  oplan <- future::plan(pulsar::modal, ref="test", ...)
  on.exit({ future::plan(oplan) ; options(oopt) })
  force(code)
}

test_that("modal backend gives the same result as the local backend", {
  skip_on_cran()
  skip_if_not_installed("processx")
  with_modal({
    for (k in seq_len(nrow(combos))) {
      refit <- combos$refit[k] ; lb <- combos$lb[k]
      expect_equal(strip(run(backend="future", refit=refit, lb.stars=lb)),
                   strip(run(refit=refit, lb.stars=lb)),
                   info=sprintf("refit=%s lb.stars=%s", refit, lb))
    }
  })
})

test_that("a lost Modal container drops its subsample with a warning", {
  skip_on_cran()
  skip_if_not_installed("processx")
  with_modal(fail=2L, {
    expect_warning(out <- run(backend="future", refit=FALSE), "container lost")
    expect_true(all(out$stars$n.fit == 5))
  })
  ## the full-data fit is submitted first: losing it stops the run
  with_modal(fail=1L, {
    expect_error(suppressWarnings(run(backend="future", refit=TRUE)),
                 "full data failed")
  })
})

test_that("more jobs than workers: each job runs once, failures are not rerun", {
  skip_on_cran()
  skip_if_not_installed("processx")
  ## workers=2 makes later jobs wait for a free worker; job 3 is lost
  with_modal(fail=3L, workers=2L, {
    expect_warning(out <- run(backend="future", refit=FALSE), "container lost")
    expect_true(all(out$stars$n.fit == 5))
    client <- get(ls(pulsar:::.modal_clients), pulsar:::.modal_clients)
    expect_equal(client$counter, 6L)
  })
})

test_that("a helper that dies drops its running job, and a new helper starts", {
  skip_on_cran()
  skip_if_not_installed("processx")
  with_modal(die=2L, workers=1L, {
    expect_warning(out <- run(backend="future", refit=FALSE), "helper process exited")
    expect_true(all(out$stars$n.fit == 5))
  })
})

test_that("changing the plan stops the Modal helper", {
  skip_on_cran()
  skip_if_not_installed("processx")
  with_modal({
    run(backend="future", refit=FALSE)
    expect_length(ls(pulsar:::.modal_clients), 1L)
    proc <- get(ls(pulsar:::.modal_clients), pulsar:::.modal_clients)$process
    future::plan(future::sequential)
    expect_length(ls(pulsar:::.modal_clients), 0L)
    expect_false(proc$is_alive())
  })
})

test_that("a helper that fails to start stops with its log", {
  skip_if_not_installed("processx")
  bad <- tempfile(fileext=".R")
  writeLines('cat("no Modal token\\n", file=stderr()); quit(status=1)', bad)
  oopt  <- options(pulsar.modal.helper=c(file.path(R.home("bin"), "Rscript"), bad))
  oplan <- future::plan(pulsar::modal, ref="test")
  on.exit({ future::plan(oplan) ; options(oopt) })
  expect_error(run(backend="future"), "failed to start.*no Modal token")
})

test_that("the Modal worker checks the pulsar version", {
  d <- tempfile() ; dir.create(d)
  saveRDS(list(data=NULL, pulsar="0.0.1"), file.path(d, "in.rds"))
  expect_error(pulsar:::.modal_eval(file.path(d, "in.rds"), file.path(d, "out.rds")),
               "0.0.1 on the controller")
})

test_that("modal settings reach the helper environment", {
  env <- pulsar:::.modalenv(list(packages=c("huge", "clime"), apt="libglpk-dev",
                                 cpu=2, memory=512L, timeout=60L, retries=1L,
                                 workers=10L, source=NULL, ref="abc123",
                                 r.version="4.6.2", app="pulsar"))
  expect_equal(env[["PULSAR_MODAL_R_PACKAGES"]], "huge,clime")
  expect_equal(env[["PULSAR_MODAL_CPU"]], "2")
  expect_equal(env[["PULSAR_MODAL_MAX_CONTAINERS"]], "10")
  expect_equal(env[["PULSAR_MODAL_REF"]], "abc123")
  expect_equal(env[["PULSAR_MODAL_SOURCE"]], "")
  oplan <- future::plan()
  on.exit(future::plan(oplan))
  expect_error(future::plan(pulsar::modal, source=tempdir()), "pulsar package directory")
})
