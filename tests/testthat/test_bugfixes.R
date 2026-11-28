context("Persistent local PSOCK pool and edge-case bug fixes")

test_that("persistent local PSOCK pool matches serial across full, lb.stars, and stop.chunk", {
  skip_on_cran()
  skip_if_not_installed("huge")
  set.seed(42)
  dat <- huge::huge.generator(n = 80, d = 15, graph = "hub", verbose = FALSE)
  lmax <- getMaxCov(dat$data)
  lams <- getLamPath(lmax, lmax * 0.1, len = 10)
  fargs <- list(lambda = lams, method = "glasso", verbose = FALSE)

  ser_full <- pulsar(dat$data, fun = huge::huge, fargs = fargs, rep.num = 6,
                     criterion = "stars", ncores = 1, seed = 101)
  par_full <- pulsar(dat$data, fun = huge::huge, fargs = fargs, rep.num = 6,
                     criterion = "stars", ncores = 2, seed = 101)
  expect_equal(par_full$stars$summary, ser_full$stars$summary)
  expect_identical(par_full$stars$opt.index, ser_full$stars$opt.index)

  ser_lb <- pulsar(dat$data, fun = huge::huge, fargs = fargs, rep.num = 6,
                   criterion = "stars", lb.stars = TRUE, ub.stars = TRUE,
                   ncores = 1, seed = 101)
  par_lb <- pulsar(dat$data, fun = huge::huge, fargs = fargs, rep.num = 6,
                   criterion = "stars", lb.stars = TRUE, ub.stars = TRUE,
                   ncores = 2, seed = 101)
  expect_equal(par_lb$stars$summary, ser_lb$stars$summary)
  expect_identical(par_lb$stars$opt.index, ser_lb$stars$opt.index)

  ser_sc <- pulsar(dat$data, fun = huge::huge, fargs = fargs, rep.num = 6,
                   criterion = "stars", stop.chunk = 2, stop.margin = 1,
                   ncores = 1, seed = 101)
  par_sc <- pulsar(dat$data, fun = huge::huge, fargs = fargs, rep.num = 6,
                   criterion = "stars", stop.chunk = 2, stop.margin = 1,
                   ncores = 2, seed = 101)
  expect_equal(par_sc$stars$summary, ser_sc$stars$summary)
  expect_identical(par_sc$stars$opt.index, ser_sc$stars$opt.index)
})

test_that("opt.index == 1 warning uses est$stars$thresh instead of hardcoded 0.1", {
  est <- list(stars = list(
    summary = c(0.05, 0.06, 0.07, 0.08),
    opt.index = 1L,
    thresh = 0.02,
    n.failed = rep(0, 4),
    n.fit = rep(10, 4)
  ))
  expect_warning(
    pulsar:::.pulsarout(est, NULL, "stars", quote(pulsar()), emptyenv()),
    "larger than the supplied values"
  )
})

test_that("print.pulsar reports Mode: future when backend = 'future'", {
  skip_if_not_installed("future")
  skip_if_not_installed("huge")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan), add = TRUE)
  set.seed(11)
  dat <- huge::huge.generator(n = 50, d = 8, graph = "hub", verbose = FALSE)
  lmax <- getMaxCov(dat$data)
  lams <- getLamPath(lmax, lmax * 0.2, len = 4)
  out <- suppressWarnings(pulsar(dat$data, fun = huge::huge,
                fargs = list(lambda = lams, method = "glasso", verbose = FALSE),
                rep.num = 4, criterion = "stars", backend = "future", seed = 1))
  txt <- paste(capture.output(print(out)), collapse = "\n")
  expect_match(txt, "Mode: future")
})

test_that("opt.index<- allows setting the last lambda index", {
  skip_if_not_installed("huge")
  set.seed(13)
  dat <- huge::huge.generator(n = 50, d = 8, graph = "hub", verbose = FALSE)
  lmax <- getMaxCov(dat$data)
  lams <- getLamPath(lmax, lmax * 0.2, len = 5)
  out <- pulsar(dat$data, fun = huge::huge,
                fargs = list(lambda = lams, method = "glasso", verbose = FALSE),
                rep.num = 4, criterion = "stars", seed = 1)
  opt.index(out, "stars") <- length(lams)
  expect_equal(opt.index(out, "stars"), length(lams))
  expect_error(opt.index(out, "stars") <- length(lams) + 1L, "Index value")
  expect_error(opt.index(out, "stars") <- NA_real_, "Index value")
})
