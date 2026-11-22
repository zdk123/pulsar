## Running jobs (#34). backend="local" is the serial / PSOCK path; with
## backend="future" each job is a future, so the user's future::plan() picks
## where it runs: plan(multisession) on this machine, plan(pulsar::modal) on
## Modal, plan(future.batchtools::batchtools_slurm) on Slurm, or any other
## future backend.

## `X` is each future's rows of the data, given in `globals` (.futurejobs)
utils::globalVariables("X")

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
## Progress goes to `tick` (see .ticker), one step per lambda fitted, with
## `phase` naming this set of jobs.
#' @keywords internal
.execjobs <- function(jobs, data, fun, fargs, backend="local", ncores=1,
                      pass.errors=TRUE, run.dir=NULL, tick=NULL, phase=NULL,
                      full.chunk=NULL) {
  ## a full-data fit split into chunks of lambdas (full.chunk) runs as several
  ## jobs, joined back into one result before the results are checked
  nparts <- 0L
  if (!is.null(full.chunk)) {
    isfull <- vapply(jobs, function(j) identical(j$kind, "full"), logical(1))
    if (any(isfull)) {
      parts  <- .splitfull(jobs[[which(isfull)]], fargs, full.chunk)
      nparts <- length(parts)
      jobs   <- c(parts, jobs[!isfull])
    }
  }
  step <- if (is.null(tick)) NULL else function(spec, status)
            tick(length(.joblambda(spec, fargs)),
                 ## the full-data fit is on the whole path, whatever the phase
                 paste0(if (length(phase) && !identical(spec$kind, "full")) paste0(phase, ": "),
                        .joblabel(spec), " ", status))
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
      .joblog(run.dir, rows[[i]])
      if (!is.null(step)) step(jobs[[i]], "loaded")
    }
  }
  if (any(todo)) {
    if (identical(backend, "future")) {
      res <- .futurejobs(jobs[todo], data, fun, fargs, run.dir, step)
    } else {
      efun <- .estfun(data, fun, run.dir)
      if (!is.null(step) && ncores <= 1) {
        ef <- efun
        efun <- function(spec, fargs) { r <- ef(spec, fargs) ; step(spec, "done") ; r }
      }
      res <- .try_mclapply(jobs[todo], efun, fargs = fargs, mc.cores = ncores,
                           mc.preschedule = FALSE, collect = FALSE)
      ## PSOCK results all come back at the end
      if (!is.null(step) && ncores > 1)
        for (k in seq_along(res$out))
          step(jobs[todo][[k]], if (inherits(res$out[[k]], "try-error")) "failed" else "done")
    }
    for (k in seq_along(res$out)) {
      i <- which(todo)[k] ; r <- res$out[[k]]
      if (inherits(r, "pulsarSaveFailed")) stop(r$message, call.=FALSE)
      if (inherits(r, "try-error")) {
        out[[i]]  <- r
        rows[[i]] <- .jobrow(jobs[[i]], fargs, "error", error=as.character(r))
      } else {
        out[[i]]  <- r$value
        rows[[i]] <- .jobrow(jobs[[i]], fargs, "done", r, r$value)
      }
    }
    warn <- res$warn
  }
  if (nparts > 1L) {
    bad <- vapply(out[1:nparts], inherits, logical(1), what="try-error")
    joined <- if (any(bad)) out[[which(bad)[1]]] else .joinfull(out[1:nparts])
    out <- c(list(joined), out[-(1:nparts)])
  }
  out <- .collectjobs(out, warn, pass.errors)
  attr(out, "jobs") <- do.call(rbind, rows)
  out
}

## The full-data job as jobs on chunks of `size` lambdas
#' @keywords internal
.splitfull <- function(spec, fargs, size) {
  lams <- .joblambda(spec, fargs)
  idx  <- split(seq_along(lams), ceiling(seq_along(lams) / size))
  if (length(idx) < 2L) return(list(spec))
  ## chunks of different sizes, so that .joinfull can tell a member with one
  ## element per lambda from a constant whose length is the chunk size
  if (length(lams) %% size == 0L) {
    last <- idx[[length(idx)]]
    idx[[length(idx)]] <- last[-length(last)]
    idx <- c(idx, list(last[length(last)]))
  }
  lapply(seq_along(idx), function(k) {
    s <- spec ; s$lambda <- lams[idx[[k]]] ; s$chunk <- k ; s
  })
}

## Join the chunked full-data fits into one output of `fun`: `path` and
## `failed`, and every other member with one element per lambda in each chunk
## (a vector or list without a class or dimensions), are concatenated, and so
## are their attributes with one element per lambda (such as timings kept on
## `path`); matrices with one column per lambda are joined by column. Other
## members are taken from the first chunk, with a warning if they differ
## between chunks (they may be per-lambda in a shape not joined here).
#' @keywords internal
.joinfull <- function(parts) {
  parts <- unname(parts)
  n   <- vapply(parts, function(x) length(x$path), integer(1))
  out <- parts[[1]]
  odd <- character(0)
  for (nm in names(out)) {
    vals <- lapply(parts, `[[`, nm)
    if (nm %in% c("path", "failed") || .perlambda(vals, n))
      out[[nm]] <- .joinvals(vals, n)
    else if (.percolumn(vals, n))
      out[[nm]] <- do.call(cbind, vals)
    else if (!all(vapply(vals[-1], identical, logical(1), vals[[1]])))
      odd <- c(odd, nm)
  }
  if (length(odd))
    warning(sprintf(paste0("full.chunk: member%s %s of fun's output differ%s between chunks ",
                           "and %s taken from the first chunk"),
                    if (length(odd) > 1) "s" else "", paste0("'", odd, "'", collapse=", "),
                    if (length(odd) > 1) "" else "s", if (length(odd) > 1) "are" else "is"),
            call.=FALSE)
  out
}

## matrices (base or Matrix) with one column per lambda in each chunk
#' @keywords internal
.percolumn <- function(vals, n) {
  dims <- lapply(vals, dim)
  all(lengths(dims) == 2L) && all(vapply(dims, `[`, integer(1), 2L) == n) &&
    length(unique(vapply(dims, `[`, integer(1), 1L))) == 1L &&
    !any(vapply(vals, is.data.frame, logical(1)))
}

#' @keywords internal
.perlambda <- function(vals, n) {
  all(lengths(vals) == n) &&
    all(vapply(vals, function(v) (is.atomic(v) || is.list(v)) && is.null(dim(v)) &&
                                 !is.object(v), logical(1)))
}

#' @keywords internal
.joinvals <- function(vals, n) {
  if (all(vapply(vals, is.null, logical(1)))) return(NULL)
  res <- do.call(c, unname(vals))
  for (a in setdiff(names(attributes(vals[[1]])), "names")) {
    av <- lapply(vals, attr, which=a, exact=TRUE)
    attr(res, a) <- if (.perlambda(av, n)) do.call(c, unname(av)) else av[[1]]
  }
  res
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
    ## an error is logged here, then goes on (a serial run stops; a PSOCK
    ## worker returns it as a failed job)
    r <- withCallingHandlers(.timejob(.run_job, spec, data, fun, fargs),
           error=function(e) .joblog(run.dir, .jobrow(spec, fargs, "error",
                                                      error=conditionMessage(e))))
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
.futurejobs <- function(jobs, data, fun, fargs, run.dir=NULL, step=NULL) {
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
      if (!is.null(step))
        step(jobs[[i]], if (inherits(r, "try-error")) "failed" else "done")
    }
    collected[done] <- TRUE
    pending <- setdiff(pending, done)
    if (length(pending) && !length(done)) Sys.sleep(0.05)
  }
  list(out=out, warn=warn)
}

## Progress of a run, in lambdas fitted: `total` counts every job on the whole
## path, and work that is skipped (bounded StARS, early stopping, dropped jobs)
## is credited with a message saying why, so the bar ends at 100%. `prog` is
## a progressr progressor, or NULL (no progress reported).
#' @keywords internal
.ticker <- function(prog, total) {
  done <- 0
  function(amount=0, msg=NULL, finish=FALSE) {
    if (finish) amount <- total - done
    amount <- max(0, min(amount, total - done))
    done <<- done + amount
    if (!is.null(prog)) prog(message=msg, amount=amount)
    invisible(done)
  }
}

#' @keywords internal
.lamrange <- function(from, to, of) {
  if (from == to) sprintf("lambda %d of %d", from, of)
  else sprintf("lambdas %d-%d of %d", from, to, of)
}
