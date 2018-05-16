context("setup")

suppressPackageStartupMessages(library(batchtools))
options(batchtools.verbose=FALSE)
options(batchtools.progress=FALSE)

tmpdir <- fs::path_real(tempdir())


context("batchtools config file")

test_that("findConfFile", {
  d   <- fs::path_temp()
  fn  <- fs::path(d, "batchtools.conf.R")
  fs::file_create(fn)
  expect_equal(findConfFile(fn), fs::path_real(fn))

  fs::file_delete(fn)

  df <- system.file('config', "batchtools.conf.R", package='pulsar')
  expect_equal(findConfFile(), fs::path_real(df))
  df <- system.file('config', "batchtools.conf.snow.R", package='pulsar')
  expect_equal(findConfFile('snow'), fs::path_real(df))
  expect_equal(findConfFile('snow.R'), fs::path_real(df))
  expect_equal(findConfFile(df), fs::path_real(df))

})


context("batchtools template file")

test_that("default template files are in the search path", {
  df <- system.file('templates', "simpletorque.tmpl", package='pulsar')
  expect_equal(findTemplateFile("simpletorque"), fs::path_real(df))
})


context('batch.pulsar options')
rseed <- 10010
p     <- 12
## generate synthetic data
set.seed(rseed)
dat <- huge::huge.generator(p*100, p, "hub", verbose=FALSE, v=.1, u=.4)
lams  <- getLamPath(getMaxCov(cor(dat$data)), 1e-1, 5)
fargs <- list(lambda=lams)

test_that("batch.pulsar options", {

  ## set wd to home
  suppressWarnings(est <- batch.pulsar(dat$data, huge::huge, fargs, rep.num=3, lb.stars=TRUE, wkdir=tmpdir))
  expect_true(file.exists(est$reg$file.dir))
  expect_true(file.exists(est$init.reg$file.dir))
  expect_equivalent(est$id$job.id, 1)
  expect_equivalent(est$init.id$job.id, 1:2)
  expect_equal(est$init.reg$work.dir, fs::path_real(tmpdir))
  expect_equal(est$reg$work.dir, fs::path_real(tmpdir))


  suppressWarnings(est <- batch.pulsar(dat$data, huge::huge, fargs, rep.num=3, lb.stars=TRUE, cleanup=TRUE, wkdir=tmpdir))
  expect_false(file.exists(est$reg$file.dir))
  expect_false(file.exists(est$init.reg$file.dir))
  expect_equivalent(est$id$job.id, 1)
  expect_equivalent(est$init.id$job.id, 1:2)
  expect_equal(est$init.reg$work.dir, fs::path_real(tmpdir))
  expect_equal(est$reg$work.dir, fs::path_real(tmpdir))

})
