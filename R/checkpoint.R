## Checkpoint and resume (#34). With run.dir, pulsar saves the run's job plan
## and each job's result as soon as it finishes. Calling pulsar again with
## the same inputs and run.dir reuses the saved plan (so the same subsamples,
## even without `seed`), loads the finished jobs and runs only the rest, on
## whatever backend is set now. Jobs that were running when the session died
## are run again.
##
## A job whose `fun` reports failed fits at some lambdas (#31) counts as
## finished: it is saved, and not run again on resume.
##
## Layout of run.dir:
##   manifest.rds     list(key, jobs, pulsar): the plan and the fingerprint of
##                    the inputs it was made for
##   job_<id>.rds     list(lambda, value, info): one finished job (info: when
##                    and where it ran); `lambda` is the
##                    path it was fitted on (bounded StARS narrows the path
##                    for stage 2, and a result is only reused for the same path)
##   job_<id>_chunk_<k>.rds
##                    the same, for one chunk of lambdas (early stopping with
##                    stop.chunk)
##   jobs.csv         the record of the jobs (see joblog.R), one row appended
##                    as each finishes

## Fingerprint of a run's inputs. Functions are compared by their code (not
## their bytecode, which changes when R compiles them), and a closure's
## captured values count too: make_fun(a=1) and make_fun(a=2) are different
## runs. The serialization header, which records the R version and locale, is
## left out, so a run can be resumed after an R upgrade.
#' @keywords internal
.runkey <- function(data, fun, fargs, ...) {
  x <- serialize(.keyable(list(data=data, fun=fun, fargs=fargs, ...)), NULL,
                 xdr=TRUE, version=2)
  f <- tempfile() ; on.exit(unlink(f))
  writeBin(x[-(1:14)], f)   # 14 bytes: "X\n" and three version integers
  unname(tools::md5sum(f))
}

## Replace functions by their code and, for closures defined inside another
## function, the captured values the code uses (other variables in that
## environment, such as a result assigned after the first run, don't count)
#' @keywords internal
.keyable <- function(x, depth=0L) {
  if (is.function(x)) {
    env <- environment(x)
    captured <- NULL
    if (depth < 3L && !is.null(env) && !.topenv(env)) {
      used <- intersect(codetools::findGlobals(x), ls(env, all.names=TRUE))
      captured <- lapply(mget(sort(used), envir=env), .keyable, depth=depth + 1L)
    }
    return(list(code=deparse(x), captured=captured))
  }
  if (is.list(x) && !is.data.frame(x)) {
    a <- attributes(x)
    x <- lapply(x, .keyable, depth=depth)
    attributes(x) <- a
  }
  x
}

#' @keywords internal
.topenv <- function(env) {
  identical(env, globalenv()) || identical(env, baseenv()) ||
    identical(env, emptyenv()) || isNamespace(env)
}

## Open run.dir: start it with this plan, or check it belongs to this run and
## return the saved plan
#' @keywords internal
.ckopen <- function(run.dir, key, jobs) {
  manifest <- file.path(run.dir, "manifest.rds")
  if (file.exists(manifest)) {
    m <- readRDS(manifest)
    if (!identical(m$key, key))
      stop("run.dir '", run.dir, "' holds a different run (other data, fun, ",
           "fargs or settings); use a new run.dir", call.=FALSE)
    ndone  <- length(list.files(run.dir, "^job_[0-9]+\\.rds$"))
    nchunk <- length(list.files(run.dir, "^job_[0-9]+_chunk_[0-9]+\\.rds$"))
    .joblogopen(run.dir)
    message(sprintf("Resuming the run in '%s': %d of %d job results saved%s",
                    run.dir, ndone, length(m$jobs),
                    if (nchunk) sprintf(", and %d results for chunks of lambdas", nchunk) else ""))
    return(m$jobs)
  }
  ## one run.dir is for one session at a time
  dir.create(run.dir, recursive=TRUE, showWarnings=FALSE)
  .ckwrite(list(key=key, jobs=jobs,
                pulsar=as.character(utils::packageVersion("pulsar"))), manifest)
  .joblogopen(run.dir)
  jobs
}

#' @keywords internal
.ckpath <- function(run.dir, spec) {
  file.path(run.dir, if (is.null(spec$chunk)) sprintf("job_%06d.rds", spec$id)
                     else sprintf("job_%06d_chunk_%03d.rds", spec$id, spec$chunk))
}

## write then rename, so a crash never leaves a partial file. A failed save
## stops the run: carrying on would lose results silently.
#' @keywords internal
.ckwrite <- function(object, path) {
  tmp <- tempfile(paste0(basename(path), "-"), tmpdir=dirname(path), fileext=".part")
  ok <- tryCatch({ saveRDS(object, tmp) ; file.rename(tmp, path) },
                 error=function(e) FALSE, warning=function(w) FALSE)
  if (!isTRUE(ok)) {
    unlink(tmp)
    stop(structure(class=c("pulsarSaveError", "error", "condition"),
                   list(message=sprintf("could not save '%s' in run.dir", path), call=NULL)))
  }
  invisible(path)
}

#' @keywords internal
.cksave <- function(run.dir, spec, value, lambda, info=NULL) {
  if (!is.null(run.dir) && !inherits(value, "try-error"))
    .ckwrite(list(lambda=lambda, value=value, info=info), .ckpath(run.dir, spec))
  invisible(value)
}

## list(found=, value=, info=): a saved result (info: its timing, see .jobinfo) for this job on this lambda path
#' @keywords internal
.ckload <- function(run.dir, spec, lambda) {
  path <- .ckpath(run.dir, spec)
  if (file.exists(path)) {
    saved <- readRDS(path)
    if (identical(saved$lambda, lambda))
      return(list(found=TRUE, value=saved$value, info=saved$info))
  }
  list(found=FALSE, value=NULL, info=NULL)
}
