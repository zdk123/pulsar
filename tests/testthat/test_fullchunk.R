context("splitting the full-data fit into chunks")

p <- 10
set.seed(4)
dat  <- huge::huge.generator(80, p, "hub", verbose=FALSE)
lams <- getLamPath(getMaxCov(dat$data), .005, 11)
n <- nrow(dat$data)

## fits each lambda on its own; besides path and failed, its output has
## per-lambda members (a vector and a list), a constant, a matrix and a
## per-lambda vector with a class, as a real fun might
coldfun <- function(countfile=tempfile(), fullfail=Inf, boom=function(lambda) FALSE) {
  function(data, lambda) {
    if (nrow(data) == n) cat(strrep("x", length(lambda)), file=countfile, append=TRUE)
    if (nrow(data) == n && boom(lambda)) stop("full fit blew up")
    ii <- match(lambda, lams)
    path <- lapply(lambda, function(l) huge::huge(data, lambda=l, verbose=FALSE)$path[[1]])
    ## per-lambda and constant attributes, as a fun keeping timings might
    attr(path, "steps") <- ii * 10L ; attr(path, "solver") <- "cd"
    list(path=path, failed=nrow(data) == n & ii >= fullfail,
         sparsity=vapply(path, function(g) sum(g) / (p * (p - 1)), numeric(1)),
         fits=lapply(ii, function(i) list(index=i)),
         method="cold", cov=stats::cov(data), k=1:4,
         df=vapply(path, function(g) as.numeric(rowSums(as.matrix(g))), numeric(p)),
         stamp=structure(ii, class="stamp"))
  }
}
nfits <- function(cf) { if (!file.exists(cf)) return(0L); k <- nchar(readLines(cf, warn=FALSE)); unlink(cf); k }
strip <- function(x) { x$call <- NULL ; x$envir <- NULL ; attr(x, "jobs") <- NULL ; x }
go <- function(fun=coldfun(), ..., rep.num=6, lambda=lams)
  suppressWarnings(pulsar(dat$data, fun=fun, fargs=list(lambda=lambda), rep.num=rep.num,
                          seed=3, job.seed=5, thresh=.05, ...))

test_that("a chunked full-data fit is the same as one fit on the whole path", {
  ref <- go()
  for (fc in c(2, 3, 4, 10, 11, 50)) {
    out <- go(full.chunk=fc)
    expect_identical(strip(out)$stars, strip(ref)$stars, info=fc)
    ## the per-lambda members are joined, the rest come from the first chunk
    expect_identical(out$est[c("path", "failed", "sparsity", "fits", "method")],
                     ref$est[c("path", "failed", "sparsity", "fits", "method")], info=fc)
    expect_identical(out$est[c("cov", "k", "df")], ref$est[c("cov", "k", "df")], info=fc)
    expect_identical(attributes(out$est$path), list(steps=1:11 * 10L, solver="cd"), info=fc)
  }
  ## a member with a class isn't joined, with a warning as it differs
  expect_warning(out <- pulsar(dat$data, fun=coldfun(), fargs=list(lambda=lams), rep.num=4,
                               seed=3, full.chunk=4), "'stamp' of fun's output differs")
  expect_identical(out$est$stamp, structure(1:4, class="stamp"))
  ## a path that divides evenly into chunks: the last lambda gets its own,
  ## so the constant k (of length 4, the chunk size) is not mistaken for a
  ## per-lambda member
  ref <- go(lambda=lams[1:8]) ; out <- go(lambda=lams[1:8], full.chunk=4)
  expect_identical(out$est$k, 1:4)
  expect_identical(out$est[c("path", "df", "sparsity")], ref$est[c("path", "df", "sparsity")])
  expect_equal(vapply(pulsar:::.splitfull(list(kind="full"), list(lambda=lams[1:8]), 4),
                      function(s) length(s$lambda), 1L), c(4, 3, 1))
})

test_that("full.chunk gives the same result with bounds, early stopping and futures", {
  for (args in list(list(lb.stars=TRUE, ub.stars=TRUE), list(stop.full=TRUE),
                    list(stop.chunk=2), list(stop.full=TRUE, stop.chunk=3))) {
    ref <- do.call(go, c(list(fun=coldfun(fullfail=9)), args))
    out <- do.call(go, c(list(fun=coldfun(fullfail=9), full.chunk=3), args))
    info <- paste(names(args), collapse="+")
    expect_identical(strip(out)$stars, strip(ref)$stars, info=info)
    expect_identical(out$est$path, ref$est$path, info=info)
    expect_identical(out$est$failed, ref$est$failed, info=info)
  }
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  ref <- go(stop.full=TRUE, stop.chunk=2)
  out <- go(stop.full=TRUE, stop.chunk=2, full.chunk=4, backend="future")
  expect_identical(strip(out)$stars, strip(ref)$stars)
  expect_identical(out$est$path, ref$est$path)
  skip_on_cran()
  skip_if(requireNamespace("pkgload", quietly=TRUE) && pkgload::is_dev_package("pulsar"),
          "multisession workers need the installed pulsar")
  future::plan(future::multisession, workers=3)
  out <- go(stop.full=TRUE, stop.chunk=2, full.chunk=4, backend="future")
  expect_identical(strip(out)$stars, strip(ref)$stars)
  expect_identical(out$est$path, ref$est$path)
})

test_that("each part is a job: recorded, resumable, and its failure stops the run", {
  cf <- tempfile() ; rd <- tempfile()
  out <- go(coldfun(cf), full.chunk=4, run.dir=rd)
  expect_equal(nfits(cf), 11)
  j <- getJobs(out)
  expect_equal(j$chunk[j$kind == "full"], 1:3)
  expect_equal(j$n.lambda[j$kind == "full"], c(4, 4, 3))
  ## lose part 2: only it is refitted
  unlink(file.path(rd, "job_000000_chunk_002.rds"))
  again <- suppressMessages(go(coldfun(cf), full.chunk=4, run.dir=rd))
  expect_equal(nfits(cf), 4)
  expect_identical(strip(again), strip(out))
  expect_error(go(coldfun(boom=function(l) min(l) < lams[6]), full.chunk=4),
               "full fit blew up|full data failed")
  skip_if_not_installed("future")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan))
  expect_error(go(coldfun(boom=function(l) min(l) < lams[6]), full.chunk=4, backend="future"),
               "full data failed")
})

test_that("full.chunk checks its argument", {
  expect_error(go(full.chunk=1), "full.chunk")
  expect_error(go(full.chunk=2.5), "full.chunk")
  expect_error(go(full.chunk=3, refit=FALSE), "refit")
})
