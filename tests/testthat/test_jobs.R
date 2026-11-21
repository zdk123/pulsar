context("job planning")

test_that(".planjobs draws the same subsamples as before and labels jobs", {
  n <- 50 ; ratio <- .8 ; rep.num <- 5
  ## the subsample draw pulsar used before jobs were explicit
  set.seed(7)
  old <- replicate(rep.num, sample(c(1L:n), floor(n * ratio), replace=FALSE), simplify=FALSE)
  for (refit in c(FALSE, TRUE)) for (lb in c(FALSE, TRUE)) {
    jobs <- pulsar:::.planjobs(n, rep.num, ratio, refit, lb, seed=7)
    info <- sprintf("refit=%s lb.stars=%s", refit, lb)
    expect_length(jobs, rep.num + refit)
    sub <- jobs[vapply(jobs, function(j) j$kind == "subsample", logical(1))]
    expect_equal(lapply(sub, `[[`, "rows"), old, info=info)
    expect_equal(vapply(sub, `[[`, integer(1), "id"), 1:rep.num, info=info)
    if (refit) {
      expect_equal(jobs[[1]]$kind, "full", info=info)
      expect_equal(jobs[[1]]$id, 0L, info=info)
      expect_equal(jobs[[1]]$rows, 1L:n, info=info)
      expect_equal(jobs[[1]]$stage, 1L, info=info)
    }
    ## bounded StARS: the first 2 subsamples (and the full fit) are stage 1
    st <- vapply(sub, `[[`, integer(1), "stage")
    expect_equal(st, if (lb) c(1L, 1L, rep(2L, rep.num - 2)) else rep(1L, rep.num), info=info)
    expect_equal(length(pulsar:::.stagejobs(jobs, 1L)) + length(pulsar:::.stagejobs(jobs, 2L)),
                 length(jobs), info=info)
  }
})

test_that(".run_job returns the full fit or the path with failures as NULL", {
  m <- Matrix::Matrix(diag(2), sparse=TRUE)
  f <- function(data, lambda) list(path=list(m, NULL), extra="x")
  full <- pulsar:::.jobspec(0, "full", 1:3)
  sub  <- pulsar:::.jobspec(1, "subsample", 1:2)
  X <- matrix(0, 3, 2)
  expect_equal(pulsar:::.run_job(full, X, f, list(lambda=1:2))$extra, "x")
  expect_equal(pulsar:::.run_job(sub, X, f, list(lambda=1:2)), list(m, NULL))
  expect_error(pulsar:::.run_job(sub, X, function(data, lambda) list(), list(lambda=1)), "path")
})
