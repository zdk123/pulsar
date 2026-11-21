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
## pass.errors, else the run stops).
#' @keywords internal
.execjobs <- function(jobs, data, fun, fargs, backend="local", ncores=1,
                      pass.errors=TRUE) {
  if (identical(backend, "future"))
    return(.futurejobs(jobs, data, fun, fargs, pass.errors))
  estFun <- function(spec, fargs) .run_job(spec, data, fun, fargs)
  .try_mclapply(jobs, estFun, fargs = fargs, mc.cores = ncores,
                mc.preschedule = FALSE, pass.errors = pass.errors)
}

## One future per job. Each future carries only its own rows of the data;
## results are collected as they finish and returned in job order. Errors,
## from `fun` or from the backend (a lost worker), become failed jobs, as on a
## PSOCK cluster; warnings are relayed as in a serial run.
#' @keywords internal
.futurejobs <- function(jobs, data, fun, fargs, pass.errors=TRUE) {
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
    for (i in done)
      out[[i]] <- withCallingHandlers(
        tryCatch(future::value(fs[[i]]),
                 error=function(e) structure(conditionMessage(e), class="try-error")),
        warning=function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
    collected[done] <- TRUE
    pending <- setdiff(pending, done)
    if (length(pending) && !length(done)) Sys.sleep(0.05)
  }
  .collectjobs(out, warn, pass.errors)
}
