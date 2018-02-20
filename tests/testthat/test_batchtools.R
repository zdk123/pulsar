context("setup")

suppressPackageStartupMessages(library(batchtools))
options(batchtools.verbose=FALSE)
options(batchtools.progress=FALSE)

fp    <- batchtools:::fp
npath <- function(x) normalizePath(x, mustWork=FALSE, winslash="/")
homedir <- npath(fp("~"))


context("batchtools config file")

conftest <- function(conf, name) {
  expect_equivalent(basename(conf), name)
  expect_equivalent(dirname(conf), system.file('config', package='pulsar'))
}

conftestbt <- function(conf) {
  exit <- file.create(conf, showWarnings=FALSE, recursive = TRUE)
  if (exit) {
    expect_equivalent(batchtools:::findConfFile(), conf)
  }
}


test_that("default config files are in search path", {
  conf0 <- "%sbatchtools.conf%s.R"
  conf1 <- pulsar::findConfFile('')
  conf2 <- pulsar::findConfFile('parallel')
  conf3 <- pulsar::findConfFile('snow')
  conf4 <- pulsar::findConfFile('torque')

  conftest(conf1, sprintf(conf0, '', ''))
  conftest(conf2, sprintf(conf0, '', '.parallel'))
  conftest(conf3, sprintf(conf0, '', '.snow'))
  conftest(conf4, sprintf(conf0, '', '.torque'))
  expect_equivalent(length(pulsar::findConfFile('toque')), 0)

  conf5 <- npath(fp("~", sprintf(conf0, '.', '')))
  dir6  <- rappdirs::user_config_dir("batchtools", expand = FALSE)
  conf6 <- fp(dir6, "config.R")
  conf7 <- npath(fp(getwd(), npath(sprintf(conf0, '', ''))))

  expect_equivalent(length(batchtools:::findConfFile()), 0)

#  conftestbt(conf5)
  dir6exists <- dir.exists(dir6)
  dir.create(dir6, showWarnings = FALSE, recursive = TRUE)
  conftestbt(conf6)
  ## remove the directory unless it already existed
  if (!dir6exists)
    unlink(dir6, recursive=TRUE, force=TRUE)
  conftestbt(conf7)

  ## cleanup
  suppressWarnings(file.remove(conf5, conf7, conf6))
})

context("batchtools template file")
templtestbt <- function(templ, name='simpletorque') {
  exit <- file.create(templ, showWarnings = FALSE, recursive = TRUE)
  if (exit) {
    expect_equivalent(pulsar::findTemplateFile(name), templ)
  }
}

test_that("default template files are in the search path", {

  expect_error(batchtools:::findTemplateFile('foo'))
  expect_error(batchtools:::findTemplateFile())
  expect_error(pulsar::findTemplateFile())

  base <- 'simpletorque'
  tmpl1 <- pulsar::findTemplateFile(base)
  expect_equivalent(basename(tmpl1), paste0(base, ".tmpl"))
  expect_equivalent(dirname(tmpl1), system.file('templates', package='pulsar'))

  x  <- sprintf("batchtools.%s.tmpl", base)

  tmpl2 <- npath(fp(system.file("templates", package = "batchtools"), paste0(base, ".tmpl")))
  tmpl3 <- npath(fp("~", sprintf(".batchtools.%s.tmpl", base)))
  dir4  <- rappdirs::user_config_dir("batchtools", expand = FALSE)
  tmpl4 <- fp(dir4, paste0(base, ".tmpl"))
  tmpl5 <- npath(fp(getwd(), sprintf("batchtools.%s.tmpl", base)))
  tmpl6 <- npath(fp(getwd(), sprintf("batchtools.%s.tmpl", 'simpletorque')))

  templtestbt(tmpl2)
  templtestbt(tmpl3)

  dir4exists <- dir.exists(dir4)
  dir.create(dir4, showWarnings = FALSE, recursive = TRUE)
  templtestbt(tmpl4)
  if (!dir4exists)
    unlink(dir4, recursive=TRUE, force=TRUE)
  templtestbt(tmpl5)
  templtestbt(tmpl6)
  ## cleanup
  suppressWarnings(file.remove(tmpl2, tmpl3, tmpl4, tmpl5, tmpl6))
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
  suppressWarnings(est <- batch.pulsar(dat$data, huge::huge, fargs, rep.num=3, lb.stars=TRUE, wkdir=homedir))
  expect_true(file.exists(est$reg$file.dir))
  expect_true(file.exists(est$init.reg$file.dir))
  expect_equivalent(est$id$job.id, 1)
  expect_equivalent(est$init.id$job.id, 1:2)
  expect_equivalent(est$init.reg$work.dir, homedir)
  expect_equivalent(est$reg$work.dir, homedir)


  suppressWarnings(est <- batch.pulsar(dat$data, huge::huge, fargs, rep.num=3, lb.stars=TRUE, cleanup=TRUE, wkdir=homedir))
  expect_false(file.exists(est$reg$file.dir))
  expect_false(file.exists(est$init.reg$file.dir))
  expect_equivalent(est$id$job.id, 1)
  expect_equivalent(est$init.id$job.id, 1:2)
  expect_equivalent(est$init.reg$work.dir, homedir)
  expect_equivalent(est$reg$work.dir, homedir)

})
