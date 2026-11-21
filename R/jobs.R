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
##   seed:  NULL, or the job's own L'Ecuyer-CMRG RNG stream (see .jobstreams)

#' @keywords internal
.jobspec <- function(id, kind, rows, stage=1L, seed=NULL) {
  list(id=as.integer(id), kind=kind, stage=as.integer(stage), rows=rows, seed=seed)
}

## Plan the jobs for a pulsar run. Subsample rows are drawn exactly as before
## the job refactor (same RNG calls, in the same order), so a given `seed`
## selects the same subsamples. With `job.seed`, each job also gets its own
## RNG stream for the random numbers `fun` draws.
#' @keywords internal
.planjobs <- function(n, rep.num, subsample.ratio, refit, lb.stars, seed=NULL,
                      job.seed=NULL) {
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
  if (!is.null(job.seed)) {
    streams <- .jobstreams(job.seed, max(vapply(jobs, `[[`, integer(1), "id")))
    for (i in seq_along(jobs)) jobs[[i]]$seed <- streams[[jobs[[i]]$id + 1L]]
  }
  jobs
}

## Independent L'Ecuyer-CMRG streams for job ids 0..maxid: the stream for job
## k is the (k+1)th after set.seed(job.seed). A job's stream depends only on
## job.seed and its id, not on the backend or the order jobs finish in. The
## caller's RNG state is left as it was.
#' @keywords internal
.jobstreams <- function(job.seed, maxid) {
  if (!is.numeric(job.seed) || length(job.seed) != 1L || is.na(job.seed))
    stop("job.seed must be NULL or a single number", call.=FALSE)
  restore <- .saverng()
  on.exit(restore())
  ## fixed normal and sample kinds, so the streams don't depend on the
  ## caller's settings
  suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"))
  set.seed(job.seed)
  s <- get(".Random.seed", envir=globalenv())
  lapply(0:maxid, function(k) s <<- parallel::nextRNGStream(s))
}

## Save the global RNG state (kind and seed); returns a function restoring it
#' @keywords internal
.saverng <- function() {
  kind <- RNGkind()
  had  <- exists(".Random.seed", envir=globalenv(), inherits=FALSE)
  seed <- if (had) get(".Random.seed", envir=globalenv(), inherits=FALSE)
  function() {
    suppressWarnings(RNGkind(kind[1], kind[2], kind[3]))  # "Rounding" warns
    if (had) assign(".Random.seed", seed, envir=globalenv())
    else if (exists(".Random.seed", envir=globalenv(), inherits=FALSE))
      rm(".Random.seed", envir=globalenv())
  }
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
  if (!is.null(spec$seed)) {
    ## run `fun` on the job's own stream, then put the caller's RNG back
    restore <- .saverng()
    on.exit(restore())
    assign(".Random.seed", spec$seed, envir=globalenv())
  }
  tmp <- do.call(fun, c(fargs, list(data[spec$rows, ])))
  if (!('path' %in% names(tmp)))
    stop('Error: expected data stucture with \'path\' member')

  if (identical(spec$kind, "full"))
    return(tmp)
  else
    return(.fitpath(tmp))
}
