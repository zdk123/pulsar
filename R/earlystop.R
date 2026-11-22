## Exact early stopping for StARS. Two opt-in modes, which combine:
##
## stop.full   fit the full data first; subsamples then fit only the lambdas
##             before the first one at which the full-data fit failed (see
##             .fitpath), since a lambda whose refit failed can't be used
## stop.chunk  fit the subsamples down the lambda path in chunks of this many
##             lambdas, and stop once the StARS summary has crossed `thresh`
##             and stop.margin more lambdas are fitted
##
## The StARS summary at a lambda depends only on the fits at that lambda, and
## its cummax up to the crossing only on the lambdas before it, so the selected
## index and the summary up to the cut are the same as with the whole path
## (same seed, job.seed, subsamples and failure rule), provided `fun`'s fit at
## a lambda doesn't depend on the other lambdas it is given, and it doesn't
## raise an error only at lambdas past the cut (see ?pulsar). With
## stop.full, lambdas from the first full-data failure on are not candidates.

#' @keywords internal
.stopcheck <- function(stop.full, stop.chunk, stop.margin, criterion, lb.stars, refit) {
  if (!is.logical(stop.full) || length(stop.full) != 1L || is.na(stop.full))
    stop("stop.full must be TRUE or FALSE", call.=FALSE)
  whole <- function(x, min) is.numeric(x) && length(x) == 1L && !is.na(x) &&
                            x >= min && x == round(x)
  if (!is.null(stop.chunk) && !whole(stop.chunk, 1))
    stop("stop.chunk must be NULL or a whole number of lambdas, at least 1", call.=FALSE)
  if (!whole(stop.margin, 0))
    stop("stop.margin must be a whole number of lambdas, at least 0", call.=FALSE)
  if (!stop.full && is.null(stop.chunk)) return(FALSE)
  if (!identical(as.character(criterion), "stars"))
    stop('early stopping (stop.full, stop.chunk) needs criterion = "stars" alone: ',
         'other criteria use every lambda', call.=FALSE)
  if (isTRUE(lb.stars))
    stop("early stopping (stop.full, stop.chunk) can't be combined with lb.stars", call.=FALSE)
  if (stop.full && !refit)
    stop("stop.full = TRUE needs refit = TRUE: the full-data fit sets the cut", call.=FALSE)
  TRUE
}

## Run the jobs with early stopping. Returns the full-data fit, the subsample
## paths up to the cut (one per job that didn't fail), and where and why the
## run stopped.
#' @keywords internal
.stoprun <- function(jobs, data, fun, fargs, backend, ncores, run.dir, p,
                     thresh, fail.action, fail.thresh, stop.full, stop.chunk,
                     stop.margin) {
  lams <- fargs$lambda
  last <- length(lams)
  full <- Filter(function(j) identical(j$kind, "full"), jobs)
  subs <- Filter(function(j) identical(j$kind, "subsample"), jobs)
  fullmodel <- NULL
  m <- NA_integer_
  run <- function(js) .execjobs(js, data, fun, fargs, backend, ncores,
                                pass.errors=TRUE, run.dir=run.dir)

  if (stop.full) {
    res <- .execjobs(full, data, fun, fargs, backend, ncores, pass.errors=FALSE,
                     run.dir=run.dir)
    fullmodel <- res[[1]]
    full <- list()
    failed <- vapply(.fitpath(fullmodel), is.null, logical(1))
    if (any(failed)) {
      m <- which(failed)[1]
      last <- m - 1L
      if (last < 1L)
        stop("The full-data fit failed at the first lambda: no lambda can be selected",
             call.=FALSE)
    }
  }

  chunk  <- if (is.null(stop.chunk)) last else as.integer(stop.chunk)
  alive  <- rep(TRUE, length(subs))      # jobs that haven't raised an error
  paths  <- rep(list(list()), length(subs))
  done   <- 0L
  target <- last
  cross  <- NA_integer_
  reason <- if (is.na(m)) "end" else "full.failed"
  k <- 0L
  while (done < target) {
    k   <- k + 1L
    idx <- (done + 1L):min(done + chunk, target)
    js  <- lapply(subs[alive], function(j) {
      j$lambda <- lams[idx]
      if (!is.null(stop.chunk)) j$chunk <- k
      j
    })
    res <- run(c(full, js))
    err <- attr(res, "errors")
    if (length(full)) {
      if (err[1]) stop("Fitting the full data failed", call.=FALSE)
      fullmodel <- res[[1]] ; res <- res[-1] ; err <- err[-1]
      full <- list()
    }
    ids <- which(alive)
    for (r in seq_along(ids[!err])) {
      i <- ids[!err][r]
      paths[[i]] <- c(paths[[i]], res[[r]])
    }
    ## a job that fails drops out at every lambda, as with the whole path
    alive[ids[err]] <- FALSE
    done <- max(idx)

    ## stop once the crossing, and the margin after it, are fitted, as long as
    ## some fitted lambda is selectable (else the pick may lie further on)
    n  <- sum(alive)
    st <- tryCatch(suppressWarnings(
            stars.stability(.tlist(paths[alive], done, n), thresh, n, p,
                            fail.action=fail.action, fail.thresh=fail.thresh)),
            error=function(e) NULL)
    cross  <- if (is.null(st)) NA_integer_ else which(st$summary >= thresh)[1]
    target <- last
    if (!is.na(cross) && any(.selectable(st$n.failed, n, fail.action, fail.thresh))) {
      target <- min(cross + as.integer(stop.margin), last)
      if (done >= target) {
        ## (a run that fitted all it could didn't stop early)
        if (done < last) reason <- "crossing"
        break
      }
    }
  }
  list(fullmodel=fullmodel, paths=paths[alive], n=sum(alive), cut=done,
       reason=reason, crossing=cross, full.failed=m, nlams=length(lams),
       chunk=if (is.null(stop.chunk)) NULL else as.integer(stop.chunk),
       margin=as.integer(stop.margin))
}

## The StARS result of an early-stopped run, on the whole lambda path: the
## summary, merge and fit counts past the cut are NA (NULL for merge)
#' @keywords internal
.stopstars <- function(es, thresh, p, fail.action, fail.thresh) {
  st <- stars.stability(.tlist(es$paths, es$cut, es$n), thresh, es$n, p,
                        fail.action=fail.action, fail.thresh=fail.thresh)
  rest <- es$nlams - es$cut
  st$summary  <- c(st$summary, rep(NA_real_, rest))
  st$merge    <- c(st$merge, vector("list", rest))
  st$n.failed <- c(st$n.failed, rep(NA_real_, rest))
  st$n.fit    <- c(st$n.fit, rep(NA_real_, rest))
  st$n.run    <- c(rep(es$n, es$cut), rep(0L, rest))
  st$stop     <- list(index=es$cut, reason=es$reason, crossing=es$crossing,
                      full.failed=es$full.failed, chunk=es$chunk, margin=es$margin)
  st
}
