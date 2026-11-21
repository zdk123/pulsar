## Running jobs (#34). backend="local" is the serial / PSOCK path; with
## backend="future" each job is a future, so the user's future::plan() picks
## where it runs: plan(multisession) on this machine, plan(pulsar::modal) on
## Modal, plan(future.batchtools::batchtools_slurm) on Slurm, or any other
## future backend.

#' @keywords internal
.backendcheck <- function(backend) {
  if (!is.character(backend) || length(backend) != 1L ||
      !(backend %in% c("local", "future")))
    stop('backend must be "local" or "future"', call.=FALSE)
  if (backend == "future" && !requireNamespace("future", quietly=TRUE))
    stop('backend="future" requires the \'future\' package', call.=FALSE)
  backend
}

## Run a set of jobs, returning their results in job order. Failed jobs are
## "try-error" strings, checked by .collectjobs (dropped with a warning when
## pass.errors, else the run stops). With run.dir, jobs already saved there
## are loaded instead of run, and each job is saved as it finishes.
#' @keywords internal
.execjobs <- function(jobs, data, fun, fargs, backend="local", ncores=1,
                      pass.errors=TRUE, run.dir=NULL) {
  out  <- vector("list", length(jobs))
  warn <- NULL
  todo <- rep(TRUE, length(jobs))
  if (!is.null(run.dir)) {
    saved <- lapply(jobs, .ckload, run.dir=run.dir, lambda=fargs$lambda)
    todo  <- !vapply(saved, `[[`, logical(1), "found")
    out[!todo] <- lapply(saved[!todo], `[[`, "value")
  }
  if (any(todo)) {
    res <- if (identical(backend, "future"))
             .futurejobs(jobs[todo], data, fun, fargs, run.dir)
           else
             .try_mclapply(jobs[todo], .estfun(data, fun, run.dir), fargs = fargs,
                           mc.cores = ncores, mc.preschedule = FALSE, collect = FALSE)
    out[todo] <- res$out
    warn <- res$warn
  }
  saveerr <- vapply(out, inherits, logical(1), what="pulsarSaveFailed")
  if (any(saveerr)) stop(out[[which(saveerr)[1]]]$message, call.=FALSE)
  .collectjobs(out, warn, pass.errors)
}

## The function a local (serial / PSOCK) worker runs for one job. Built here
## so that PSOCK workers get only data, fun and run.dir in its environment.
#' @keywords internal
.estfun <- function(data, fun, run.dir) {
  force(data) ; force(fun) ; force(run.dir)
  function(spec, fargs) {
    value <- .run_job(spec, data, fun, fargs)
    ## a failed save must not look like a failed fit (PSOCK workers turn
    ## errors into failed jobs), so it comes back as a marker
    tryCatch(.cksave(run.dir, spec, value, fargs$lambda),
             pulsarSaveError=function(e)
               structure(list(message=conditionMessage(e)), class="pulsarSaveFailed"))
  }
}

## One future per job. Each future carries only its own rows of the data;
## results are collected as they finish and returned in job order. Errors,
## from `fun` or from the backend (a lost worker), become failed jobs, as on a
## PSOCK cluster; warnings are relayed as in a serial run. Returns the raw
## results and warnings, for .execjobs to check.
#' @keywords internal
.futurejobs <- function(jobs, data, fun, fargs, run.dir=NULL) {
  run_job <- .run_job
  fs <- vector("list", length(jobs))
  ## on an error or interrupt, cancel the jobs still running
  collected <- logical(length(jobs))
  on.exit(for (i in which(!collected & !vapply(fs, is.null, logical(1))))
            try(future::cancel(fs[[i]]), silent=TRUE), add=TRUE)
  for (i in seq_along(jobs)) {
    spec <- jobs[[i]]
    rows <- spec$rows
    spec$rows <- seq_along(rows)
    ## seed=NULL: no seed and no RNG check, as with backend="local"
    fs[[i]] <- future::future(run_job(spec, X, fun, fargs),
                   globals=list(run_job=run_job, spec=spec,
                                X=data[rows, , drop=FALSE], fun=fun, fargs=fargs),
                   packages="pulsar", seed=NULL, lazy=FALSE,
                   label=sprintf("pulsar-job-%d", spec$id))
  }

  out  <- vector("list", length(fs))
  warn <- NULL
  pending <- seq_along(fs)
  while (length(pending)) {
    done <- pending[vapply(fs[pending], future::resolved, logical(1))]
    for (i in done) {
      out[[i]] <- withCallingHandlers(
        tryCatch(future::value(fs[[i]]),
                 error=function(e) structure(conditionMessage(e), class="try-error")),
        warning=function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      .cksave(run.dir, jobs[[i]], out[[i]], fargs$lambda)
    }
    collected[done] <- TRUE
    pending <- setdiff(pending, done)
    if (length(pending) && !length(done)) Sys.sleep(0.05)
  }
  list(out=out, warn=warn)
}
