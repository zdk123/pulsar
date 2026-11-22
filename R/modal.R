## Modal future backend (#34). pulsar runs its jobs as futures, so any future
## backend works (plan(multisession), future.batchtools on Slurm, ...). This
## file adds one for Modal (https://modal.com), which has no shared filesystem
## and no inbound connection to the controller: futures are submitted to, and
## collected from, a long-lived helper process (inst/modal/pulsar_modal.py)
## that holds an ephemeral Modal app for the R session. It needs nothing from
## pulsar beyond .modal_eval() on the worker, so it can move to its own
## package later.

#' Modal futures
#'
#' A \pkg{future} backend that evaluates each future in its own container on \href{https://modal.com}{Modal}. Use it to run pulsar's subsample fits in parallel on Modal: \code{future::plan(pulsar::modal, cpu=1, memory=2048)}, then \code{pulsar(..., backend="future")}. It works for other futures too.
#'
#' Requirements on the controlling machine: the \pkg{future} and \pkg{processx} packages, the \code{modal} Python package (\code{pip install modal}, or \code{uv tool install modal}), and a Modal token, set up once with \code{modal setup} (or the \code{MODAL_TOKEN_ID} and \code{MODAL_TOKEN_SECRET} environment variables). pulsar never reads or stores the token.
#'
#' Workers run a \code{rocker/r-ver} image (or \code{image} / \code{dockerfile}) with \pkg{future}, pulsar, and the packages that \code{fun} needs: those in \code{install}, from CRAN, Bioconductor or GitHub, and the local package directories in \code{install.local}. The image is built on first use and cached by Modal, so later sessions with the same settings start quickly. Changing a file in an \code{install.local} package rebuilds the image from that package on. A GitHub branch is cached at its first build: pin a commit (\code{"owner/repo@<sha>"}), or set \code{rebuild=TRUE} to pick up new commits. pulsar on the workers is built from \code{source}, a local package directory, or installed from GitHub at \code{ref}; with neither, the commit the controller's pulsar was installed from (with \code{remotes::install_github}) is used, and the worker checks that its pulsar version matches the controller's.
#'
#' The first future starts a helper process that opens an ephemeral Modal app; it is stopped by \code{plan()} changing to another backend, or when R exits. A future whose container fails (lost worker, timeout, out of memory) is retried up to \code{retries} times, then raises a \code{FutureError}; an error raised by the future's own code is returned as usual.
#'
#' @param install packages to install on the workers. (The name isn't \code{packages}, which \code{plan()} reserves for \code{future()}.) Each is a CRAN package name (installed from Posit Package Manager binaries), \code{"bioc::name"} for Bioconductor, or \code{"owner/repo@ref"} (also \code{"owner/repo/subdir@ref"}) for GitHub. GitHub and local packages are installed with \pkg{remotes}, along with their dependencies (\code{Depends}, \code{Imports}, \code{LinkingTo}, but not \code{Suggests}) from CRAN or Bioconductor. A private GitHub repository needs a token: put \code{GITHUB_PAT} in a Modal secret and name it in \code{secrets}.
#' @param install.local paths to local R package source directories (such as your own unpublished packages) to install on the workers. They are installed in dependency order, after pulsar and the GitHub packages (so they may depend on those, but not the reverse), and each worker checks it has the version in the directory's \code{DESCRIPTION}. Compiled objects (\code{*.o}, \code{*.so}), \code{.git}, \code{renv} and check output are not uploaded.
#' @param secrets names of \href{https://modal.com/docs/guide/secrets}{Modal secrets} whose environment variables are available while packages are installed (for example one holding \code{GITHUB_PAT}). They are not passed to the futures.
#' @param image a container image with R (and \code{Rscript}) to start from, instead of \code{rocker/r-ver:<r.version>}
#' @param dockerfile path to a Dockerfile building such an image, instead of \code{image}
#' @param rebuild rebuild the worker image from scratch (for example to pick up new commits on a GitHub branch)
#' @param apt system libraries to install on the workers with \code{apt-get} (the default, \code{libglpk-dev}, is needed by \pkg{igraph}, which \pkg{huge} imports)
#' @param cpu number of CPU cores per future
#' @param memory memory per future, in MiB
#' @param timeout time limit per future, in seconds
#' @param retries number of times to retry a future whose container failed
#' @param workers maximum number of futures running at once
#' @param source path to a local pulsar source directory to install on the workers
#' @param ref a GitHub ref (branch, tag or commit) of \code{zdk123/pulsar} to install on the workers
#' @param r.version the \code{rocker/r-ver} image tag
#' @param app name of the Modal app. Several R sessions can use the same name at once; when they start together, Modal's limit on how fast apps are created can refuse some, and those wait and retry (for up to a few minutes)
#' @param python path to a Python with the \code{modal} package. By default, the \code{PULSAR_MODAL_PYTHON} environment variable, else the Python that runs the \code{modal} command line tool.
#' @param verbose show Modal's progress output (image builds, container logs)
#' @param envir unused
#' @param \dots further arguments for \code{\link[future]{FutureBackend}}
#' @seealso \code{\link{pulsar}}, \code{\link[future]{plan}}
#' @examples
#' \dontrun{
#' ## requires a Modal account and the modal Python package
#' future::plan(pulsar::modal, install=c("huge", "bioc::SpiecEasi"), cpu=1, memory=1024)
#' ## your own package, and a private GitHub one (with a Modal secret holding GITHUB_PAT)
#' future::plan(pulsar::modal, install=c("huge", "me/privatepkg@<sha>"),
#'              install.local="~/src/mypkg", secrets="github")
#' p <- 40 ; dat <- MASS::mvrnorm(100, rep(0, p), diag(p))
#' lams <- getLamPath(getMaxCov(dat), .05, 10)
#' out <- pulsar(dat, fun=huge::huge, fargs=list(lambda=lams), rep.num=50,
#'               backend="future")
#' future::plan(future::sequential)   # stops the Modal app
#' }
#' @export
modal <- function(install="huge", install.local=NULL, apt="libglpk-dev", secrets=NULL,
                  cpu=1, memory=2048, timeout=3600, retries=0L, workers=100L,
                  source=NULL, ref=NULL, r.version="4.6.2", image=NULL,
                  dockerfile=NULL, rebuild=FALSE, app="pulsar", python=NULL,
                  verbose=FALSE, ..., envir=parent.frame()) {
  stop("INTERNAL ERROR: pulsar::modal() is a future plan; use future::plan(pulsar::modal)")
}
class(modal) <- c("modal", "multiprocess", "future", "function")
attr(modal, "init") <- TRUE

## The future backend made by plan(modal, ...)
#' @keywords internal
ModalFutureBackend <- function(install="huge", install.local=NULL, apt="libglpk-dev",
                               secrets=NULL, cpu=1, memory=2048, timeout=3600,
                               retries=0L, workers=100L, source=NULL, ref=NULL,
                               r.version="4.6.2", image=NULL, dockerfile=NULL,
                               rebuild=FALSE, app="pulsar", python=NULL,
                               verbose=FALSE, ..., envir=NULL) {
  if (!is.null(source)) {
    source <- normalizePath(source, mustWork=TRUE)
    if (!file.exists(file.path(source, "DESCRIPTION")))
      stop("modal: source must be a pulsar package directory", call.=FALSE)
  }
  if (!is.null(dockerfile)) dockerfile <- normalizePath(dockerfile, mustWork=TRUE)
  if (length(secrets) && !all(grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", secrets)))
    stop("modal: secret names may only use letters, digits, '_', '.' and '-'", call.=FALSE)
  if (!is.null(ref) && (!grepl(.SPEC_REF, ref) || grepl("..", ref, fixed=TRUE)))
    stop("modal: ref must be a branch, tag or commit name", call.=FALSE)
  settings <- list(packages=.modalspecs(install), local=.modallocal(install.local),
                   apt=as.character(apt), secrets=as.character(secrets), cpu=cpu,
                   memory=as.integer(memory), timeout=as.integer(timeout),
                   retries=as.integer(retries), workers=as.integer(workers),
                   source=source, ref=ref, r.version=r.version, image=image,
                   dockerfile=dockerfile, rebuild=isTRUE(rebuild), app=app,
                   python=python, verbose=verbose)
  core <- future::FutureBackend(reg=paste0("workers-modal-", .modalkey(settings)),
                                workers=as.integer(workers), modal=settings, ...,
                                timeout=getOption("future.wait.timeout", 30*24*60*60),
                                delta=getOption("future.wait.interval", 0.2),
                                alpha=getOption("future.wait.alpha", 1.01))
  core[["futureClasses"]] <- c("ModalFuture", "MultiprocessFuture", core[["futureClasses"]])
  structure(core, class=c("ModalFutureBackend", "MultiprocessFutureBackend",
                          "FutureBackend", class(core)))
}
attr(modal, "factory") <- ModalFutureBackend

## internal functions of 'future' that backends need (as future.mirai does)
#' @keywords internal
.future_internal <- function(name) get(name, envir=asNamespace("future"), inherits=FALSE)

## Names the helper for these settings. Local packages count by their file
## sizes and modification times, so editing one starts a new helper (and so
## a new image) instead of reusing the running one.
#' @keywords internal
.modalkey <- function(settings) {
  f <- tempfile() ; on.exit(unlink(f))
  stamp <- lapply(c(settings$source, settings$local$dirs), function(d) {
    files <- .modalfiles(d)
    info  <- file.info(file.path(d, files))[, c("size", "mtime")]
    paste(files, info$size, as.numeric(info$mtime))
  })
  if (!is.null(settings$dockerfile))
    stamp <- c(stamp, tools::md5sum(settings$dockerfile))
  writeLines(c(utils::capture.output(dput(settings[setdiff(names(settings), "verbose")])),
               unlist(stamp)), f)
  substr(unname(tools::md5sum(f)), 1, 10)
}

## The files of a package directory that are uploaded to the image: all but
## those pulsar_modal.py's _IGNORE skips (relative paths)
#' @keywords internal
.modalfiles <- function(d) {
  files <- list.files(d, recursive=TRUE, all.files=TRUE, no..=TRUE)
  top   <- sub("/.*", "", files)
  parts <- strsplit(files, "/", fixed=TRUE)
  skip  <- top %in% c(".git", "doc", "Meta", ".Rproj.user", "renv", ".Rhistory") |
           grepl("\\.Rcheck$", top) | (top == files & grepl("\\.tar\\.gz$", files)) |
           grepl("\\.(o|so|dll)$", files) | basename(files) == ".DS_Store" |
           vapply(parts, function(x) "__pycache__" %in% x, logical(1))
  files[!skip]
}

## Package specs are pasted into R code on the image builder, so only these
## shapes are allowed (the same patterns as in pulsar_modal.py)
.SPEC_NAME   <- "[A-Za-z][A-Za-z0-9.]*"
.SPEC_CRAN   <- paste0("^", .SPEC_NAME, "$")
.SPEC_BIOC   <- paste0("^bioc::", .SPEC_NAME, "$")
.SPEC_GITHUB <- "^[A-Za-z0-9][A-Za-z0-9-]*/[A-Za-z0-9_-][A-Za-z0-9_.-]*(/[A-Za-z0-9_][A-Za-z0-9_./-]*)?(@[A-Za-z0-9_][A-Za-z0-9_./-]*)?$"
.SPEC_REF    <- "^[A-Za-z0-9_][A-Za-z0-9_./-]*$"

#' @keywords internal
.modalspecs <- function(packages) {
  packages <- as.character(packages)
  ok <- (grepl(.SPEC_CRAN, packages) | grepl(.SPEC_BIOC, packages) |
         grepl(.SPEC_GITHUB, packages)) & !grepl("..", packages, fixed=TRUE)
  if (!all(ok))
    stop("modal: not a CRAN name, \"bioc::name\" or \"owner/repo@ref\": ",
         paste(sQuote(packages[!ok], FALSE), collapse=", "),
         ". Local package directories go in `install.local`.", call.=FALSE)
  packages
}

## Local package directories: their names and versions, in an order where
## each comes after the local packages it depends on
#' @keywords internal
.modallocal <- function(local) {
  if (!length(local)) return(NULL)
  if (any(grepl("[\n\r]", local)))
    stop("modal: install.local paths may not contain line breaks", call.=FALSE)
  dirs <- normalizePath(as.character(local), mustWork=TRUE)
  desc <- lapply(dirs, function(d) {
    f <- file.path(d, "DESCRIPTION")
    if (!file.exists(f))
      stop("modal: no DESCRIPTION in install.local directory ", d, call.=FALSE)
    read.dcf(f, fields=c("Package", "Version", "Depends", "Imports", "LinkingTo"))
  })
  pkgs <- vapply(desc, function(x) x[, "Package"], character(1))
  if (anyDuplicated(pkgs))
    stop("modal: the same package twice in `install.local`: ", pkgs[duplicated(pkgs)], call.=FALSE)
  needs <- lapply(desc, function(x) {
    deps <- unlist(strsplit(paste(stats::na.omit(x[, c("Depends", "Imports", "LinkingTo")]),
                                  collapse=","), ","))
    intersect(trimws(sub("\\(.*", "", deps)), pkgs)
  })
  order <- character(0)
  while (length(order) < length(pkgs)) {
    ready <- setdiff(pkgs[vapply(needs, function(n) all(n %in% order), logical(1))], order)
    if (!length(ready))
      stop("modal: circular dependencies among `install.local` packages", call.=FALSE)
    order <- c(order, ready)
  }
  i <- match(order, pkgs)
  list(dirs=dirs[i], packages=order,
       ## as packageVersion() reports them on the worker: "0.9-3" is "0.9.3"
       versions=stats::setNames(vapply(desc[i], function(x)
         as.character(package_version(x[, "Version"])), character(1)), order))
}

#' @exportS3Method future::nbrOfWorkers
nbrOfWorkers.ModalFutureBackend <- function(evaluator) evaluator[["workers"]]

#' @exportS3Method future::nbrOfFreeWorkers
nbrOfFreeWorkers.ModalFutureBackend <- function(evaluator, background=FALSE, ...) {
  max(0L, evaluator[["workers"]] - .modalrunning(evaluator))
}

## A launched future is done when its result (or failure) file exists, or
## when the helper that would write it has exited
#' @keywords internal
.modaldone <- function(m) {
  file.exists(m$outpath) || file.exists(paste0(m$outpath, ".err")) ||
    !m$client$process$is_alive()
}

## Number of futures still running on Modal. Finished ones stay registered
## until their value is collected, but no longer hold a worker.
#' @keywords internal
.modalrunning <- function(backend) {
  FutureRegistry <- .future_internal("FutureRegistry")
  futures <- FutureRegistry(backend[["reg"]], action="list", earlySignal=FALSE)
  sum(vapply(futures, function(f) !is.null(f[["modal"]]) && !.modaldone(f[["modal"]]),
             logical(1)))
}

#' @exportS3Method future::launchFuture
launchFuture.ModalFutureBackend <- function(backend, future, ...) {
  FutureRegistry <- .future_internal("FutureRegistry")
  getFutureData  <- .future_internal("getFutureData")
  .modalwait(backend)
  client <- .modalclient(backend[["modal"]])

  key <- sprintf("f%06d", client$counter <- client$counter + 1L)
  inpath  <- file.path(client$dir, paste0(key, ".in.rds"))
  outpath <- file.path(client$dir, paste0(key, ".out.rds"))
  saveRDS(list(data=getFutureData(future), pulsar=client$version,
               packages=client$expect), inpath)
  future[["backend"]] <- backend
  future[["modal"]]   <- list(key=key, inpath=inpath, outpath=outpath, client=client)
  future[["state"]]   <- "submitted"
  res <- .modalrequest(client, c("submit", key, inpath, outpath), timeout=600)
  unlink(inpath)
  if (res[1] != "ok")
    stop(future::FutureError(paste("Failed to submit future to Modal:",
                                   paste(res[-1], collapse=" ")), future=future))
  future[["modal"]]$call <- res[2]
  future[["state"]] <- "running"
  FutureRegistry(backend[["reg"]], action="add", future=future, earlySignal=FALSE)
  invisible(future)
}

## Wait until fewer than `workers` futures are running. This only watches
## for result files; it doesn't collect values (FutureRegistry's
## "collect-first" would reset a failed future, which then runs again).
#' @keywords internal
.modalwait <- function(backend) {
  interval <- backend[["delta"]]
  t0 <- Sys.time()
  while (.modalrunning(backend) >= backend[["workers"]]) {
    if (difftime(Sys.time(), t0, units="secs") > backend[["timeout"]])
      stop(future::FutureError(sprintf("TIMEOUT: all %d Modal workers still busy",
                                       backend[["workers"]])))
    Sys.sleep(interval)
    interval <- min(interval * backend[["alpha"]], 5)
  }
}

#' @exportS3Method future::resolved
resolved.ModalFuture <- function(x, ...) {
  if (NextMethod()) return(TRUE)
  if (x[["state"]] == "created") {
    future::run(x)
    return(FALSE)
  }
  .modaldone(x[["modal"]])
}

#' @exportS3Method future::result
result.ModalFuture <- function(future, ...) {
  FutureRegistry <- .future_internal("FutureRegistry")
  if (inherits(future[["result"]], "FutureResult")) return(future[["result"]])
  if (inherits(future[["result"]], "condition")) stop(future[["result"]])

  m <- future[["modal"]]
  errpath <- paste0(m$outpath, ".err")
  while (!file.exists(m$outpath) && !file.exists(errpath)) {
    if (!m$client$process$is_alive())
      break
    Sys.sleep(0.1)
  }
  reg <- future[["backend"]][["reg"]]
  on.exit(if (FutureRegistry(reg, action="contains", future=future))
            FutureRegistry(reg, action="remove", future=future), add=TRUE)

  if (file.exists(m$outpath)) {
    result <- readRDS(m$outpath)
    unlink(m$outpath)
    future[["result"]] <- result
    future[["state"]]  <- "finished"
    return(result)
  }
  msg <- if (file.exists(errpath)) paste(readLines(errpath, warn=FALSE), collapse="\n")
         else paste("the Modal helper process exited:", .modallog(m$client))
  unlink(errpath)
  future[["state"]]  <- "failed"
  future[["result"]] <- future::FutureError(
    sprintf("Modal future (%s) failed: %s", m$call %||% m$key, msg), future=future)
  stop(future[["result"]])
}

#' @exportS3Method future::interruptFuture
interruptFuture.ModalFutureBackend <- function(backend, future, ...) {
  FutureRegistry <- .future_internal("FutureRegistry")
  m <- future[["modal"]]
  if (!is.null(m) && !.modaldone(m))
    .modalrequest(m$client, c("cancel", m$key))
  future[["state"]]  <- "interrupted"
  future[["result"]] <- future::FutureInterruptError(
    sprintf("Modal future (%s) was interrupted", m$call %||% m$key), future=future)
  if (FutureRegistry(backend[["reg"]], action="contains", future=future))
    FutureRegistry(backend[["reg"]], action="remove", future=future)
  invisible(future)
}

#' @exportS3Method future::stopWorkers
stopWorkers.ModalFutureBackend <- function(backend, ...) {
  FutureRegistry <- .future_internal("FutureRegistry")
  ## cancel what is still running before the helper (and its app) stops
  for (f in FutureRegistry(backend[["reg"]], action="list", earlySignal=FALSE))
    if (!is.null(f[["modal"]]) && !.modaldone(f[["modal"]]))
      interruptFuture.ModalFutureBackend(backend, f)
  FutureRegistry(backend[["reg"]], action="reset")
  .modalstop(.modalkey(backend[["modal"]]))
  invisible(backend)
}

#' @exportS3Method future::tweak
tweak.modal <- function(strategy, ..., penvir=parent.frame()) {
  attr(strategy, "init") <- TRUE
  backend <- attr(strategy, "backend", exact=TRUE)
  if (!is.null(backend)) {
    future::stopWorkers(backend)
    attr(strategy, "backend") <- NULL
  }
  NextMethod("tweak")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

## ---- helper process ------------------------------------------------------

.modal_clients <- new.env(parent=emptyenv())

## The pulsar ref the workers install, when neither source nor ref is given
#' @keywords internal
.modalref <- function(settings) {
  if (!is.null(settings$source) || !is.null(settings$ref)) return(settings$ref)
  sha <- utils::packageDescription("pulsar")$RemoteSha
  if (is.null(sha))
    stop("modal: the installed pulsar did not come from GitHub; set `source` ",
         "(a local pulsar directory) or `ref` in plan(pulsar::modal, ...)", call.=FALSE)
  sha
}

#' @keywords internal
.modalenv <- function(settings) {
  ref <- .modalref(settings)
  c(PULSAR_MODAL_R_VERSION      = settings$r.version,
    PULSAR_MODAL_R_PACKAGES     = paste(settings$packages, collapse="\n"),
    PULSAR_MODAL_LOCAL          = paste(settings$local$dirs, collapse="\n"),
    PULSAR_MODAL_SECRETS        = paste(settings$secrets, collapse=","),
    PULSAR_MODAL_IMAGE          = settings$image %||% "",
    PULSAR_MODAL_DOCKERFILE     = settings$dockerfile %||% "",
    MODAL_FORCE_BUILD           = if (isTRUE(settings$rebuild)) "1" else "0",
    PULSAR_MODAL_APT            = paste(settings$apt, collapse=","),
    PULSAR_MODAL_SOURCE         = settings$source %||% "",
    PULSAR_MODAL_REF            = ref %||% "",
    PULSAR_MODAL_CPU            = as.character(settings$cpu),
    PULSAR_MODAL_MEMORY         = as.character(settings$memory),
    PULSAR_MODAL_TIMEOUT        = as.character(settings$timeout),
    PULSAR_MODAL_RETRIES        = as.character(settings$retries),
    PULSAR_MODAL_MAX_CONTAINERS = as.character(settings$workers),
    PULSAR_MODAL_APP            = settings$app)
}

## The Python that runs the modal CLI (pip, pipx and uv installs all put it in
## the script's shebang line)
#' @keywords internal
.modalpython <- function(python=NULL) {
  if (!is.null(python)) return(python)
  env <- Sys.getenv("PULSAR_MODAL_PYTHON")
  if (nzchar(env)) return(env)
  cli <- Sys.which("modal")
  if (nzchar(cli)) {
    first <- tryCatch(readLines(cli, n=1L, warn=FALSE), error=function(e) "")
    if (length(first) && grepl("^#!.*python", first))
      return(trimws(sub("^#!", "", first)))
  }
  py <- Sys.which("python3")
  if (!nzchar(py))
    stop("modal: no Python found; install the modal package (pip install modal) ",
         "or set PULSAR_MODAL_PYTHON", call.=FALSE)
  unname(py)
}

## The command that starts the helper; option pulsar.modal.helper overrides it
## (the tests use a local stand-in)
#' @keywords internal
.modalhelper <- function(settings) {
  getOption("pulsar.modal.helper",
            c(.modalpython(settings$python),
              system.file("modal", "pulsar_modal.py", package="pulsar", mustWork=TRUE),
              "serve"))
}

## Start (or reuse) the helper process for these settings
#' @keywords internal
.modalclient <- function(settings) {
  if (!requireNamespace("processx", quietly=TRUE))
    stop("modal: the 'processx' package is required", call.=FALSE)
  key <- .modalkey(settings)
  client <- .modal_clients[[key]]
  if (!is.null(client) && client$process$is_alive()) return(client)

  dir <- tempfile("pulsar-modal-") ; dir.create(dir)
  cmd <- .modalhelper(settings)
  log <- file.path(dir, "helper.log")
  proc <- processx::process$new(cmd[1], cmd[-1], stdin="|", stdout="|",
            stderr=if (isTRUE(settings$verbose)) "" else log,
            env=c("current", .modalenv(settings)), cleanup=TRUE)
  client <- new.env(parent=emptyenv())
  client$process <- proc ; client$dir <- dir ; client$log <- log
  client$counter <- 0L ; client$buffer <- character(0)
  client$version <- as.character(utils::packageVersion("pulsar"))
  client$expect  <- settings$local$versions
  reg.finalizer(client, function(e) .modalquit(e), onexit=TRUE)

  ## the first start builds the worker image, which can take a few minutes
  ready <- .modalreadline(client, timeout=getOption("pulsar.modal.start.timeout", 1800))
  if (!identical(ready, "ready")) {
    proc$kill()
    on.exit(unlink(dir, recursive=TRUE), add=TRUE)
    stop("modal: the helper process failed to start:\n", .modallog(client), call.=FALSE)
  }
  .modal_clients[[key]] <- client
  client
}

#' @keywords internal
.modalreadline <- function(client, timeout=60) {
  t0 <- Sys.time()
  repeat {
    if (length(client$buffer)) {
      line <- client$buffer[1] ; client$buffer <- client$buffer[-1]
      return(line)
    }
    left <- timeout - as.numeric(difftime(Sys.time(), t0, units="secs"))
    if (left <= 0) return(NULL)
    io <- client$process$poll_io(as.integer(min(left, 1) * 1000))
    lines <- if (io[["output"]] == "ready") client$process$read_output_lines()
             else character(0)
    ## at end of output poll_io keeps reporting "ready" with nothing to read
    if (!length(lines) && !client$process$is_alive()) return(NULL)
    client$buffer <- c(client$buffer, lines)
  }
}

#' @keywords internal
.modalrequest <- function(client, fields, timeout=120) {
  if (!client$process$is_alive())
    return(c("error", paste("the Modal helper process exited:", .modallog(client))))
  client$process$write_input(paste0(paste(fields, collapse="\t"), "\n"))
  line <- .modalreadline(client, timeout)
  if (is.null(line)) {
    ## a late reply would be read as the answer to the next request, so a
    ## helper that doesn't answer in time is stopped
    client$process$kill()
    return(c("error", "no reply from the Modal helper process"))
  }
  strsplit(line, "\t", fixed=TRUE)[[1]]
}

#' @keywords internal
.modallog <- function(client) {
  if (!file.exists(client$log))
    return("(no log: with verbose=TRUE the helper's output is printed above)")
  paste(utils::tail(readLines(client$log, warn=FALSE), 30), collapse="\n")
}

#' @keywords internal
.modalquit <- function(client) {
  if (!is.null(client$process) && client$process$is_alive()) {
    try(client$process$write_input("quit\n"), silent=TRUE)
    client$process$wait(30000)
    if (client$process$is_alive()) client$process$kill()
  }
  unlink(client$dir, recursive=TRUE)
}

#' @keywords internal
.modalstop <- function(key) {
  client <- .modal_clients[[key]]
  if (!is.null(client)) {
    .modalquit(client)
    rm(list=key, envir=.modal_clients)
  }
}

## ---- worker --------------------------------------------------------------

## Entry point in a Modal container: evaluate the future saved in `inpath`
## and save its FutureResult to `outpath`. The worker's pulsar must match the
## controller's, since pulsar's job code travels by reference to its namespace.
#' @keywords internal
.modal_eval <- function(inpath, outpath) {
  payload <- readRDS(inpath)
  here <- as.character(utils::packageVersion("pulsar"))
  if (!identical(here, payload$pulsar))
    stop(sprintf("pulsar %s on the Modal worker, %s on the controller", here, payload$pulsar))
  ## local packages: the worker must have the version in their DESCRIPTION
  for (pkg in names(payload$packages)) {
    have <- tryCatch(as.character(utils::packageVersion(pkg)), error=function(e) "none")
    if (!identical(have, payload$packages[[pkg]]))
      stop(sprintf("%s %s on the Modal worker, %s in its local directory",
                   pkg, have, payload$packages[[pkg]]))
  }
  evalFuture <- .future_internal("evalFuture")
  saveRDS(evalFuture(payload$data), outpath)
  invisible(NULL)
}
