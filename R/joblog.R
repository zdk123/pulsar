## The record of a run's jobs: one row per job (per chunk, with early
## stopping), kept as the "jobs" attribute of a pulsar result (see getJobs)
## and, with run.dir, appended to run.dir/jobs.csv as each job finishes (and
## as each saved job is loaded), so it outlives a session that dies mid-run.
##
##   id, kind, stage, chunk  the job (see jobs.R)
##   n.lambda, lambda.max, lambda.min
##                           the lambdas it fitted
##   status                  "done", "error" (fun raised an error, or the
##                           backend lost the worker) or "loaded" (from run.dir)
##   start, secs, host       when it started, how long it took, and where it
##                           ran (for loaded jobs, from the run that fitted it)
##   failed.fits             the number of lambdas at which fun reported a
##                           failed fit
##   error                   the error message, for failed jobs

#' @keywords internal
.jobrow <- function(spec, fargs, status, res=NULL, value=NULL, error=NA_character_) {
  lam <- .joblambda(spec, fargs)
  data.frame(id=spec$id, kind=spec$kind, stage=spec$stage,
             chunk=if (is.null(spec$chunk)) NA_integer_ else as.integer(spec$chunk),
             n.lambda=length(lam), lambda.max=max(lam), lambda.min=min(lam),
             status=status,
             start=if (is.null(res$start)) as.POSIXct(NA) else res$start,
             secs=if (is.null(res$end)) NA_real_
                  else as.numeric(difftime(res$end, res$start, units="secs")),
             host=if (is.null(res$host)) NA_character_ else res$host,
             failed.fits=if (is.null(value)) NA_integer_ else .nfailedfits(spec, value),
             error=error, stringsAsFactors=FALSE)
}

## the timing kept with a saved job
#' @keywords internal
.jobinfo <- function(r) r[c("start", "end", "host")]

#' @keywords internal
.nfailedfits <- function(spec, value) {
  path <- if (identical(spec$kind, "full")) tryCatch(.fitpath(value), error=function(e) NULL)
          else value
  if (is.null(path)) NA_integer_ else sum(vapply(path, is.null, logical(1)))
}

#' @keywords internal
.joblabel <- function(spec) {
  if (identical(spec$kind, "full")) "full-data fit" else paste("subsample", spec$id)
}

#' @keywords internal
.joblogpath <- function(run.dir) file.path(run.dir, "jobs.csv")

## Start run.dir/jobs.csv with its header, unless it exists
#' @keywords internal
.joblogopen <- function(run.dir) {
  f <- .joblogpath(run.dir)
  if (!file.exists(f)) {
    spec <- .jobspec(0L, "full", integer(0))
    try(utils::write.table(.jobrow(spec, list(lambda=1), "done")[0, ], f, sep=",",
                           row.names=FALSE, qmethod="double"), silent=TRUE)
  }
}

## Append a row to run.dir/jobs.csv. The log is a record, not a result: a
## failed write doesn't stop the run.
#' @keywords internal
.joblog <- function(run.dir, row) {
  if (is.null(run.dir)) return(invisible())
  row$start <- format(row$start, "%Y-%m-%dT%H:%M:%OS3Z", tz="UTC")
  row$error <- gsub("[\r\n]+", " ", row$error)
  try(utils::write.table(row, .joblogpath(run.dir), sep=",", append=TRUE,
                         col.names=FALSE, row.names=FALSE, qmethod="double"), silent=TRUE)
  invisible()
}

#' The record of a pulsar run's jobs
#'
#' Each pulsar run is a set of jobs: the fit on the full data, and one fit per subsample (per chunk of lambdas, with early stopping). \code{getJobs} gives a record of them, for checking where the time went and which jobs failed. With \code{run.dir}, the same record is appended to \code{jobs.csv} in \code{run.dir} as each job finishes, so it is kept even if the session dies mid-run (with times in UTC, and line breaks in error messages replaced by spaces). Each call of \code{pulsar} appends its own rows, including \code{"loaded"} rows for the jobs a resumed run reuses, and an \code{"error"} row for a job whose error stopped a serial run (which \code{getJobs} can't show, as there is no result).
#'
#' @param x a \code{pulsar} object
#' @return a data frame with one row per job, or \code{NULL} for a \code{batch.pulsar} object. The columns are:
#' \itemize{
#'   \item id: 0 for the full-data fit, 1 to \code{rep.num} for the subsamples
#'   \item kind: \code{"full"} or \code{"subsample"}
#'   \item stage: 1, or 2 for the subsamples fitted on the narrowed path with \code{lb.stars}
#'   \item chunk: the chunk of lambdas, with \code{stop.chunk} (else \code{NA})
#'   \item n.lambda, lambda.max, lambda.min: the lambdas the job fitted
#'   \item status: \code{"done"}, \code{"error"} (\code{fun} raised an error, or the backend lost the worker) or \code{"loaded"} (a result saved in \code{run.dir} by an earlier call)
#'   \item start, secs, host: when the job started, how many seconds it took and the machine it ran on (\code{NA} for failed jobs; for loaded jobs, from the call that fitted them)
#'   \item failed.fits: the number of lambdas at which \code{fun} reported a failed fit
#'   \item error: the error message of a failed job
#' }
#' @seealso \code{\link{pulsar}}
#' @export
getJobs <- function(x) attr(x, "jobs")
