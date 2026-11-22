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
## are loaded instead of run, and each job is saved as it finishes. The
## result's "jobs" attribute is the record of the jobs (see .jobrow), which
## with run.dir is also appended to run.dir/jobs.csv as each job finishes.
## Progress is reported through progressr, when it is installed.
#' @keywords internal
.execjobs <- function(jobs, data, fun, fargs, backend="local", ncores=1,
                      pass.errors=TRUE, run.dir=NULL) {
  out  <- vector("list", length(jobs))
  rows <- vector("list", length(jobs))
  warn <- NULL
  todo <- rep(TRUE, length(jobs))
  if (!is.null(run.dir)) {
    saved <- lapply(jobs, function(j) .ckload(run.dir, j, .joblambda(j, fargs)))
    todo  <- !vapply(saved, `[[`, logical(1), "found")
    for (i in which(!todo)) {
      out[[i]]  <- saved[[i]]$value
      rows[[i]] <- .jobrow(jobs[[i]], fargs, "loaded", saved[[i]]$info, out[[i]])
    }
  }
  if (any(todo)) {
    ## progressr's progressor belongs to the function that creates it
    prog <- if (requireNamespace("progressr", quietly=TRUE))
              progressr::progressor(steps=sum(todo)) else NULL
    if (identical(backend, "future")) {
      res <- .futurejobs(jobs[todo], data, fun, fargs, run.dir, prog)
    } else {
      efun <- .estfun(data, fun, run.dir)
      if (!is.null(prog) && ncores <= 1) {
        ef <- efun
        efun <- function(spec, fargs) { r <- ef(spec, fargs) ; prog(.joblabel(spec)) ; r }
      }
      res <- .try_mclapply(jobs[todo], efun, fargs = fargs, mc.cores = ncores,
                           mc.preschedule = FALSE, collect = FALSE)
      if (!is.null(prog) && ncores > 1) for (j in jobs[todo]) prog(.joblabel(j))
    }
    for (k in seq_along(res$out)) {
      i <- which(todo)[k] ; r <- res$out[[k]]
      if (inherits(r, "pulsarSaveFailed")) stop(r$message, call.=FALSE)
      if (inherits(r, "try-error")) {
        out[[i]]  <- r
        rows[[i]] <- .jobrow(jobs[[i]], fargs, "error", error=as.character(r))
        ## futures are logged as they are collected; local jobs by the worker
        if (!identical(backend, "future")) .joblog(run.dir, rows[[i]])
      } else {
        out[[i]]  <- r$value
        rows[[i]] <- .jobrow(jobs[[i]], fargs, "done", r, r$value)
      }
    }
    warn <- res$warn
  }
  out <- .collectjobs(out, warn, pass.errors)
  attr(out, "jobs") <- do.call(rbind, rows)
  out
}

## Run one job and time it: list(value, start, end, host)
#' @keywords internal
.timejob <- function(run_job, spec, data, fun, fargs) {
  start <- Sys.time()
  value <- run_job(spec, data, fun, fargs)
  list(value=value, start=start, end=Sys.time(), host=Sys.info()[["nodename"]])
}

## The function a local (serial / PSOCK) worker runs for one job. Built here
## so that PSOCK workers get only data, fun and run.dir in its environment.
#' @keywords internal
.estfun <- function(data, fun, run.dir) {
  force(data) ; force(fun) ; force(run.dir)
  function(spec, fargs) {
    r <- .timejob(.run_job, spec, data, fun, fargs)
    ## a failed save must not look like a failed fit (PSOCK workers turn
    ## errors into failed jobs), so it comes back as a marker
    tryCatch({
      .cksave(run.dir, spec, r$value, .joblambda(spec, fargs), .jobinfo(r))
      .joblog(run.dir, .jobrow(spec, fargs, "done", r, r$value))
      r
    }, pulsarSaveError=function(e)
         structure(list(message=conditionMessage(e)), class="pulsarSaveFailed"))
  }
}

## One future per job. Each future carries only its own rows of the data;
## results are collected as they finish and returned in job order. Errors,
## from `fun` or from the backend (a lost worker), become failed jobs, as on a
## PSOCK cluster; warnings are relayed as in a serial run. Returns the raw
## (timed) results and warnings, for .execjobs to check.
#' @keywords internal
.futurejobs <- function(jobs, data, fun, fargs, run.dir=NULL, prog=NULL) {
  run_job <- .run_job ; time_job <- .timejob
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
    fs[[i]] <- future::future(time_job(run_job, spec, X, fun, fargs),
                   globals=list(time_job=time_job, run_job=run_job, spec=spec,
                                X=data[rows, , drop=FALSE], fun=fun, fargs=fargs),
                   packages="pulsar", seed=NULL, lazy=FALSE,
                   label=sprintf("pulsar-%s", gsub(" ", "-", .joblabel(spec))))
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
      r <- out[[i]]
      if (inherits(r, "try-error")) {
        .joblog(run.dir, .jobrow(jobs[[i]], fargs, "error", error=as.character(r)))
      } else {
        .cksave(run.dir, jobs[[i]], r$value, .joblambda(jobs[[i]], fargs), .jobinfo(r))
        .joblog(run.dir, .jobrow(jobs[[i]], fargs, "done", r, r$value))
      }
      if (!is.null(prog)) prog(.joblabel(jobs[[i]]))
    }
    collected[done] <- TRUE
    pending <- setdiff(pending, done)
    if (length(pending) && !length(done)) Sys.sleep(0.05)
  }
  list(out=out, warn=warn)
}
