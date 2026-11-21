## Job planning and execution (see #34)
##
## A job is one call of the graph estimator `fun`:
##   id:    0 for the full-data fit, 1..rep.num for subsamples
##   kind:  "full" (refit on all rows; returns the whole `fun` output) or
##          "subsample" (returns the path, with failed fits as NULL)
##   stage: 1 or 2. Bounded StARS (lb.stars) fits the full data and the first
##          2 subsamples in stage 1, and the rest on a reduced lambda path in
##          stage 2. Without lb.stars every job is in stage 1.
##   rows:  row indices of `data` the job is fitted on

#' @keywords internal
.jobspec <- function(id, kind, rows, stage=1L) {
  list(id=as.integer(id), kind=kind, stage=as.integer(stage), rows=rows)
}

## Plan the jobs for a pulsar run. Subsample rows are drawn exactly as before
## the job refactor (same RNG calls, in the same order), so a given `seed`
## selects the same subsamples.
#' @keywords internal
.planjobs <- function(n, rep.num, subsample.ratio, refit, lb.stars, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  rows <- replicate(rep.num,
            sample(c(1L:n), floor(n * subsample.ratio), replace = FALSE),
            simplify=FALSE)
  if (!is.null(seed)) set.seed(NULL)

  ## the bounds use the first 2 subsamples (see .lbcheck)
  nstage1 <- if (lb.stars) 2L else rep.num
  jobs <- lapply(seq_len(rep.num), function(i)
            .jobspec(i, "subsample", rows[[i]], stage=if (i <= nstage1) 1L else 2L))
  if (refit)
    jobs <- c(list(.jobspec(0L, "full", 1L:n, stage=1L)), jobs)
  jobs
}

#' @keywords internal
.stagejobs <- function(jobs, stage) {
  jobs[vapply(jobs, function(j) j$stage == stage, logical(1))]
}

## Run one job: fit `fun` on the job's rows. A full-data job returns the whole
## output of `fun`; a subsample job returns its path, with failed fits as NULL
## (see .fitpath).
#' @keywords internal
.run_job <- function(spec, data, fun, fargs) {
  tmp <- do.call(fun, c(fargs, list(data[spec$rows, ])))
  if (!('path' %in% names(tmp)))
    stop('Error: expected data stucture with \'path\' member')

  if (identical(spec$kind, "full"))
    return(tmp)
  else
    return(.fitpath(tmp))
}
