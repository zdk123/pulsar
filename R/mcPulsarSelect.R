#' @keywords internal
.lamcheck <- function(lams) {
    if (is.null(lams)) {
        stop(paste('Error: missing members in fargs:',
             paste(c('lambda')[c(is.null(lams))])))
    } else {
        if (!all(lams == cummin(lams)))
            warning("Are you sure you don't want the lambda path to be monotonically decreasing")
        if (length(lams) < 2)
            warning("Only 1 value of lambda is given. Are you sure you want to do model selection?")
    }
}

#' @keywords internal
.ratcheck <- function(subsample.ratio, n) {
    if (is.null(subsample.ratio)) {
        if (n > 144)
            return(10 * sqrt(n)/n)
        else
            return(0.8)
    } else return(subsample.ratio)
}

#' @keywords internal
.lbcheck <- function(lb.stars, rep.num) {
    ## the bounds use the first 2 subsamples; the rest are fitted afterwards
    if (isTRUE(lb.stars) && rep.num <= 2)
        stop("lb.stars = TRUE needs rep.num > 2: the bounds are computed from the first 2 subsamples",
             call.=FALSE)
}

#' @keywords internal
.critcheck0 <- function(criterion, knowncrits) {
    if (!all(criterion %in% knowncrits)) {
       stop(paste('Unknown criterion', paste(criterion[!(criterion %in% knowncrits)],
                   collapse=", "), sep=": "))
    }
    starsdepend <- c("estrada", "sufficiency")
    if (any(starsdepend %in% knowncrits)) {
        if (any(starsdepend %in% criterion) && !("stars" %in% criterion)) {
             stop(paste('Criterion: ', paste(starsdepend[starsdepend %in% criterion],
                   collapse=", "), ' cannot be run unless stars is also a selected criterion', sep=""))
        }
    }

}

#' @importFrom parallel mclapply
#' @keywords internal
.tlist <- function(li, n, m) {
  if (missing(m)) m <- length(li)
  lapply(1L:n, function(i) {
    lapply(1L:m, function(j) li[[j]][[i]])
  })
}

#' @keywords internal
.try_mclapply <- function(X, FUN, mc.cores = getOption("mc.cores", 2L),
                          pass.errors=TRUE, mc.preschedule = TRUE, ...) {
  warn <- NULL
  env <- environment()

  if (mc.cores > 1L) {
    cl <- parallel::makePSOCKcluster(mc.cores)
    on.exit(parallel::stopCluster(cl))
    safe_fun <- function(x, ...) {
      tryCatch(FUN(x, ...), error = function(e) {
        structure(conditionMessage(e), class = "try-error")
      })
    }
    if (mc.preschedule) {
      out <- parallel::parLapply(cl, X, safe_fun, ...)
    } else {
      out <- parallel::parLapplyLB(cl, X, safe_fun, ...)
    }
  } else {
    withCallingHandlers({out <- mclapply(X, FUN, mc.cores=1L, ...)},
                        warning=function(w) {
                           assign('warn', c(warn, w$message), env)
                           invokeRestart("muffleWarning")
                         })
  }

  errors <- sapply(out, inherits, what="try-error")
  if (any(errors)) {
    error <- table(trimws(sapply(out[errors], '[', 1)))
    msg   <- paste(sprintf('\n%s job%s failed with: "%s"',
                  error, ifelse(error>1, 's', ''), names(error)), collapse="" )
    if (!pass.errors || all(errors)) {
      stop(msg, call.=FALSE)
    } else {
      warning(msg, call.=FALSE)
      out <- out[!errors]
    }
  } else {
    if (length(warn) > 0) {
      twarn <- table(trimws(sapply(warn, '[', 1)))
      msg   <- paste(sprintf('%s job%s had warning: "%s"',
                              twarn, ifelse(twarn>1, 's', ''), names(twarn)),
                      collapse="\n")
      warning(msg, call.=FALSE)
    }
  }

  attr(out, 'errors') <- errors
  return(out)
}


#' pulsar: serial or parallel mode
#'
#' Run pulsar using StARS' edge stability (or other criteria) to select an undirected graphical model over a lambda path.
#'
#' @param data A \eqn{n*p} matrix of data matrix input to solve for the \eqn{p*p} graphical model
#' @param fun pass in a function that returns a list representing \eqn{p*p} sparse, undirected graphical models along the desired regularization path. The expected inputs to this function are: a data matrix input and a sequence of decreasing lambdas and must return a list or S3 object with a member \emph{named} \code{path}. This should be a list of adjacency matrices for each value of \code{lambda}. See \code{\link{pulsar-function}} for more information.
#' @param fargs arguments to argument \code{fun}. Must be a named list and requires at least one member \code{lambda}, a numeric vector with values for the penalty parameter.
#' @param criterion A character vector of selection statistics. Multiple criteria can be supplied. Only StARS can be used to automatically select an optimal index for the lambda path. See details for additional statistics.
#' @param thresh threshold (referred to as scalar \eqn{\beta} in StARS publication) for selection criterion. Only implemented for StARS. \code{thresh=0.1} is recommended.
#' @param subsample.ratio determine the size of the subsamples (referred to as \eqn{b(n)/n}). Default is 10*sqrt(n)/n for n > 144 or 0.8 otherwise. Should be strictly less than 1.
#' @param rep.num number of random subsamples \eqn{N} to take for graph re-estimation. Default is \eqn{N=20}, but more is recommended for non-StARS criteria or if using edge frequencies as confidence scores.
#' @param seed A numeric seed to force predictable subsampling. Default is NULL. Use for testing purposes only.
#' @param lb.stars Should the lower bound be computed after the first \eqn{N=2} subsamples (should result in considerable speedup and only implemented if stars is selected). If this option is selected, other summary metrics will only be applied to the smaller lambda path.
#' @param ub.stars Should the upper bound be computed after the first \eqn{N=2} subsamples (should result in considerable speedup and only implemented if stars is selected). If this option is selected, other summary metrics will only be applied to the smaller lambda path. This option is ignored if the lb.stars flag is FALSE.
#' @param ncores number of cores to use for subsampling. When \code{ncores > 1}, a PSOCK cluster is used (via \code{\link[parallel]{parLapply}}) to avoid fork-safety issues with OpenMP libraries. Workers receive \code{fun} and \code{fargs}, but not other objects from your session, so a \code{fun} that calls your own helper functions should define them inside \code{fun} or load them from a package. See \code{batch.pulsar} for more parallelization options.
#' @param refit Boolean flag to refit on the full dataset after pulsar is run. (see also \code{\link{refit}})
#' @param fail.action how to treat subsample fits that \code{fun} reports as failed (see Details and \code{\link{pulsar-function}}). With \code{"drop"} (default), failed fits are left out and selection frequencies at each lambda are computed over the successful fits only. With \code{"unselectable"}, failed fits are also left out and, in addition, a lambda whose fraction of failed fits exceeds \code{fail.thresh} cannot be selected by StARS; the nearest selectable lambda with a larger penalty is chosen instead. A lambda at which every fit failed is never selectable.
#' @param fail.thresh threshold, in \eqn{[0, 1)}, on the fraction of failed subsample fits above which a lambda is unselectable when \code{fail.action="unselectable"}. The default, \eqn{0}, makes any lambda with a failed fit unselectable.
#'
#' @return an S3 object of class \code{pulsar} with a named member for each stability metric run. Within each of these are:
#' \itemize{
#'    \item summary: the summary statistic over \code{rep.num} graphs at each value of lambda
#'    \item criterion: the stability criterion used
#'    \item merge: the raw statistic over the \code{rep.num} graphs, prior to summarization
#'    \item opt.ind: index (along the path) of optimal lambda selected by the criterion at the desired threshold. Will return \eqn{0} if no optimum is found or \code{NULL} if selection for the criterion is not implemented.
#'   }
#' If \code{stars} is included as a criterion then additional arguments include
#' \itemize{
#'    \item lb.index: the lambda index of the lower bound at \eqn{N=2} samples if \code{lb.stars} flag is set to TRUE
#'    \item ub.index: the lambda index of the upper bound at \eqn{N=2} samples if \code{ub.stars} flag is set to TRUE
#'    \item n.failed: the number of subsample fits at each lambda that \code{fun} reported as failed
#'    \item n.fit: the number of successful subsample fits at each lambda (the denominator of the selection frequencies in \code{merge})
#'}
#' @return call: the original function call
#' @details
#' \strong{Failed fits.} If the graphical model function fails on a subsample at some lambda without raising an error (e.g. the solver does not converge), \pkg{pulsar} can only leave that fit out if \code{fun} reports the failure, by returning \code{NULL} or \code{NA} for that element of \code{path}, or by returning a logical vector \code{failed} alongside \code{path} (see \code{\link{pulsar-function}}). An unreported failure that returns an empty or diagonal estimate is counted as an empty graph, which biases the edge selection frequencies toward zero, and hence the StARS (and G-StARS) instability, at that lambda. Since failures tend to happen at small lambda, where StARS selects, \code{fun} should signal non-convergence whenever it can. A warning is given when reported failures occur at, or next to, the selected lambda.
#'

#' The options for \code{criterion} statistics are:
#' \itemize{
#'    \item stars (Stability approach to regularization selection)
#'    \item gcd   (Graphet correlation distance, requires the \pkg{orca} package) see \code{\link{gcvec}}
#'    \item diss  (Node-node dissimilarity) see \code{\link{graph.diss}}
#'    \item estrada (estrada class) see \code{\link{estrada.class}}
#'    \item nc  (natural connectivity) see \code{\link{natural.connectivity}}
#'    \item sufficiency (Tandon & Ravikumar's sufficiency statistic)
#' }
#' @examples
#'\dontrun{
#' ## Generate the data with huge:
#' library(huge)
#' p <- 40 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(getMaxCov(dat$data), .01, len=20)
#'
#' ## Run pulsar with huge
#' hugeargs <- list(lambda=lams, verbose=FALSE)
#' out.p <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                 rep.num=20, criterion='stars')
#'
#' ## Run pulsar in bounded stars mode and include gcd metric:
#' out.b <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                 rep.num=20, criterion=c('stars', 'gcd'),
#'                 lb.stars=TRUE, ub.stars=TRUE)
#' plot(out.b)
#' }
#' @importFrom Matrix mean triu
#' @importFrom parallel mclapply
#' @references Müller, C. L., Bonneau, R., & Kurtz, Z. (2016). Generalized Stability Approach for Regularized Graphical Models. arXiv. https://arxiv.org/abs/1605.07072
#' @references Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. Proceedings of the Twenty-Third Annual Conference on Neural Information Processing Systems (NIPS).
#' @references Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012). The huge Package for High-dimensional Undirected Graph Estimation in R. The Journal of Machine Learning Research, 13, 1059–1062.
#' @seealso \code{\link{batch.pulsar}} \code{\link{refit}}
#' @export
pulsar <- function(data, fun=huge::huge, fargs=list(),
                   criterion=c("stars"),
                   thresh = 0.1, subsample.ratio = NULL,
                   rep.num = 20, seed=NULL,
                   lb.stars=FALSE, ub.stars=FALSE,
                   ncores = 1, refit=TRUE,
                   fail.action = c("drop", "unselectable"), fail.thresh = 0)  {
  gcinfo(FALSE)
  ## evaluate fun now: parallel workers get the function itself, not a promise
  ## pointing at the caller's (unserialized) global environment (#43)
  force(fun)
  n <- nrow(data)
  p <- ncol(data)
  # min requirements for function args
  .lamcheck(fargs$lambda)
  nlams <- length(fargs$lambda)
  knowncrits <- c("stars", "diss", "estrada", "gcd", "nc", "sufficiency") #"vgraphlet", "egraphlet",
  .critcheck0(criterion, knowncrits)
  subsample.ratio <- .ratcheck(subsample.ratio, n)
  fail.action <- .failcheck(fail.action, fail.thresh)
  .lbcheck(lb.stars, rep.num)

  if (!is.null(seed)) set.seed(seed)
  ind.sample <- replicate(rep.num,
                  sample(c(1L:n), floor(n * subsample.ratio),
                  replace = FALSE), simplify=FALSE)
  if (refit) {
    tmp <- 1L:n
    attr(tmp, 'full') <- TRUE
    ind.sample <- c(list(tmp), ind.sample)
  }
  if (!is.null(seed)) set.seed(NULL)
  ## wrap estimator
  estFun <- function(ind.sample, fargs) {
    tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
    if (!('path' %in% names(tmp)))
      stop('Error: expected data stucture with \'path\' member')

    if (isTRUE(attr(ind.sample, 'full')))
      return(tmp)
    else
      return(.fitpath(tmp))
  }

  if (lb.stars) {
    if (!("stars" %in% criterion))
      stop('Lower/Upper bound method must be used with StARS')
    minN <- 2L + refit # Hard code for now
  } else minN <- rep.num + refit

  isamp <- ind.sample[1L:minN]
  ## don't pass on errors if lb.stars = TRUE
  premerge <- .try_mclapply(isamp, estFun, fargs = fargs, mc.cores = ncores,
                            mc.preschedule = FALSE, pass.errors = !lb.stars)
  errors <- attr(premerge, 'errors')
  # Adjust rep.num for failed jobs
  rep.num <- rep.num - ifelse(refit, sum(errors[-1]), sum(errors))

  if (refit && errors[1])
    stop('Fitting the full data failed', call.=FALSE)

  if (refit) {
    fullmodel <- premerge[[1]]
    premerge  <- premerge[-1]
    minN <- minN - 1
  } else fullmodel <- NULL

  if (lb.stars) {
    lb.premerge       <- premerge
    lb.premerge.reord <- lapply(1L:nlams, function(i)
                         lapply(1L:minN, function(j) lb.premerge[[j]][[i]]))

    lb.est <- stars.stability(lb.premerge.reord, thresh, minN, p)
    if (lb.est$opt.index == 1)
      warning("Accurate lower bound could not be determined with the first 2 subsamples")
    if (ub.stars) {
      # upper bound is determined by equivilent of MaxEnt of Poisson Binomial
      pmean      <- sapply(lb.est$merge, function(x) { sum(x)/(p*(p-1)) })
      ub.summary <- cummax(4*pmean*(1-pmean))
      tmpub      <- .starsind(ub.summary, thresh, 1)
      if (any(ub.summary == 0))  ## adjust upper bound to exclude empty graphs
        ub.index <- max(tmpub, max(which(ub.summary == 0))+1)
      else
        ub.index <- max(tmpub, 1)
    } else ub.index <- 1
    # reselect lambda between bounds
    fargs$lambda <- fargs$lambda[ub.index:lb.est$opt.index]
    nlams <- length(fargs$lambda)
    lb.premerge  <- lapply(lb.premerge,
                           function(ppm) ppm[ub.index:lb.est$opt.index])
    isamp <- ind.sample[-(1L:(minN+refit))]
#    tmp   <- mclapply(isamp, estFun, fargs=fargs, mc.cores=ncores, mc.preschedule = FALSE)
    tmp <- .try_mclapply(isamp, estFun, fargs = fargs, mc.cores = ncores,
                         mc.preschedule = FALSE)
    # Adjust rep.num for failed jobs
    rep.num <- rep.num - sum(attr(tmp, 'errors'))
    premerge <- c(lb.premerge, tmp)
  }

  premerge.reord <- .tlist(premerge, nlams, rep.num)
  rm(premerge) ; gc()
  est <- list()

  for (i in 1L:length(criterion)) {
    crit <- criterion[i]

    if (crit == "stars")
      est$stars <- stars.stability(premerge.reord, thresh, rep.num, p,
                                   fail.action=fail.action, fail.thresh=fail.thresh)

    else if (crit == "diss")
      est$diss <-  diss.stability(premerge.reord, thresh, rep.num, p, nlams)

    else if (crit == "estrada") {
      if (!("stars" %in% criterion))
        warning('Need StaRS for computing Estrada classes... not run')
      else
        est$estrada <- estrada.stability(est$stars$merge,thresh,rep.num,p,nlams)
    }

    else if (crit == "sufficiency") {
      if (!("stars" %in% criterion)) warning('Need StaRS for computing sufficiency... not run')
      else  est$sufficiency <- sufficiency(est$stars$merge, rep.num, p, nlams)
    }

#      else if (crit == "egraphlet")
#        est$egraphlet <- egraphlet.stability(premerge.reord, thresh, rep.num, p, nlams)

#      else if (crit == "vgraphlet")
#        est$vgraphlet <- vgraphlet.stability(premerge.reord, thresh, rep.num, p, nlams)

    else if (crit == "gcd") {
      est$gcd <- gcd.stability(premerge.reord, thresh, rep.num, p, nlams)

    } else if (crit == "nc")
      est$nc <- nc.stability(premerge.reord, thresh, rep.num, p, nlams)

  }

  if (lb.stars) {
    find <- 1:length(lb.est$summary)
    pind <- ub.index:lb.est$opt.index
    pinv <- setdiff(find, pind)
    tmpsumm <- vector('numeric', length(lb.est$summary))
    tmpsumm[pinv] <- lb.est$summary[pinv]
    tmpsumm[pind] <- est$stars$summary
    est$stars$summary   <- tmpsumm

    tmpmerg <- vector('list', length(lb.est$summary))
    tmpmerg[pinv]   <- lb.est$merge[pinv]
    tmpmerg[pind]   <- est$stars$merge
    est$stars$merge <- tmpmerg

    tmpfail <- lb.est$n.failed ; tmpfit <- lb.est$n.fit
    tmpfail[pind]      <- est$stars$n.failed
    tmpfit[pind]       <- est$stars$n.fit
    est$stars$n.failed <- tmpfail
    est$stars$n.fit    <- tmpfit

    est$stars$lb.index  <- lb.est$opt.index
    est$stars$ub.index  <- ub.index
    est$stars$opt.index <- est$stars$opt.index + ub.index - 1L
  }

  if ("stars" %in% criterion) {
      .failwarn(est$stars$n.failed, est$stars$n.failed + est$stars$n.fit,
                est$stars$opt.index)
      if (est$stars$opt.index == 1) {
          direction <- if (any(est$stars$summary >= .1)) "larger" else "smaller"
          warning(paste("Optimal lambda may be", direction, "than the supplied values"))
      }
  }
  est$call  <- match.call()
  est$envir <- parent.frame()
  est$est   <- fullmodel
  return(structure(est, class="pulsar"))
}

#' @keywords internal
.starsind <- function(summary, thresh, offset=1, selectable=NULL) {
  if(any(summary >= thresh)){
    ind <- max(which.max(summary >= thresh)[1] - offset, 1)
  } else {
    warning("Optimal lambda may be outside the supplied values")
    ind <- length(summary)
  }
  if (is.null(selectable) || selectable[ind]) return(ind)
  if (!any(selectable)) {
    warning("No lambda is selectable given the failed subsample fits; returning the StARS index anyway", call.=FALSE)
    return(ind)
  }
  cand  <- which(selectable)
  lower <- cand[cand < ind]
  new   <- if (length(lower)) max(lower) else min(cand)
  warning(sprintf("Lambda index %d is unselectable due to failed subsample fits; selecting index %d instead",
                  ind, new), call.=FALSE)
  new
}

#' @keywords internal
stars.stability <- function(premerge, stars.thresh, rep.num, p, merge=NULL,
                            n.failed=NULL, fail.action="drop", fail.thresh=0) {
    if (is.null(stars.thresh)) stars.thresh <- 0.1
    est <- list()

    if (is.null(merge)) {
      ## failed fits are NULL entries; drop them from the per-lambda sums
      n.failed  <- vapply(premerge, function(x) sum(vapply(x, is.null, logical(1))), numeric(1))
      est$merge <- lapply(premerge, function(x) {
        x <- .dropfailed(x)
        if (length(x)) Reduce("+", x) else NULL
      })
      gc() # flush
    } else {
      est$merge <- merge
      if (is.null(n.failed)) n.failed <- numeric(length(merge))
    }
    nlams   <- length(est$merge)
    n.fit   <- rep.num - n.failed
    est$summary <- rep(0, nlams)

    ok <- !vapply(est$merge, is.null, logical(1)) & n.fit > 0
    if (!any(ok)) stop("All subsample fits failed at every lambda", call.=FALSE)
    template <- est$merge[[which(ok)[1]]]
    for (i in 1:nlams) {
      if (!ok[i]) {
        ## no successful fits: frequencies are unknown, keep the dimensions but
        ## make sure this lambda can't be selected (summary is carried by cummax)
        est$merge[[i]] <- template * 0
        next
      }
      est$merge[[i]] <- est$merge[[i]]/n.fit[i]
      est$summary[i] <- 4 * sum(est$merge[[i]] * (1 - est$merge[[i]])) / (p * (p - 1))
    }
    if (any(!ok))
      warning(sprintf("All subsample fits failed at lambda index %s",
                      paste(which(!ok), collapse=", ")), call.=FALSE)
    ## monotonize variability
    est$summary   <- cummax(est$summary)
    selectable    <- .selectable(n.failed, rep.num, fail.action, fail.thresh)
    est$opt.index <- .starsind(est$summary, stars.thresh, selectable=selectable)
    est$criterion <- "stars.stability"
    est$thresh    <- stars.thresh
    est$n.failed  <- n.failed
    est$n.fit     <- n.fit
    return(est)
}

#' @keywords internal
sufficiency <- function(merge, rep.num, p, nlams) {
## Merge solution from StARS
  est <- list()
  est$merge <- sapply(merge, function(x) apply(x*(1-x), 2, max))
  est$summary <- colMeans(est$merge)
  est$criterion <- 'sufficiency'
  return(est)
}

#' @keywords internal
.sumsq <- function(x2,y) x2 + y^2

#' @keywords internal
diss.stability <- function(premerge, diss.thresh, rep.num, p, nlams) {
    est <- list()
    premerge  <- lapply(premerge, .dropfailed)
    n.fit     <- vapply(premerge, length, numeric(1))
    disslist  <- lapply(premerge, function(pm) lapply(pm, graph.diss))
    est$merge <- lapply(seq_along(disslist), function(i)
                        Reduce("+", disslist[[i]])/n.fit[i])
    mergesq   <- lapply(seq_along(disslist), function(i)
                        Reduce(.sumsq, disslist[[i]][-1], init=disslist[[i]][[1]]^2)/n.fit[i])

    gc() # flush
    est$summary <- rep(0, length(est$merge))
    for (i in 1:length(est$merge)) {
        if (n.fit[i] == 0) { est$summary[i] <- NA ; next }
        tmp <- mergesq[[i]] - est$merge[[i]]^2
        est$summary[i] <- sum(triu(tmp))  / (p * (p - 1))
    }
    est$mergesq <- mergesq
    est$criterion <- "diss.stability"
    est$n.failed  <- rep.num - n.fit
    return(est)
}

#estrada.stability <- function(premerge, thresh, rep.num, p, nlams) {
#    est <- list()
#    estrlist  <- lapply(premerge, function(pm) lapply(pm, estrada.class))
#    est$merge <- lapply(estrlist, function(x) table(unlist(x)))

##    gc() # flush
#    est$summary <- rep(0, length(est$merge))
#    for (i in 1:length(est$merge)) {
#        est$summary[i] <- 1-max(est$merge[[i]])/rep.num
#    }
#    ## monotonize variability
##    est$summary <- cummax(est$summary)
#    if (!is.null(thresh))
#      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
#    else
#      est$opt.index <- 0

#    est$criterion <- "estrada.stability"
#    return(est)
#}

#' @keywords internal
estrada.stability <- function(merge, thresh, rep.num, p, nlams) {
    est <- list()
    est$summary <- unlist(lapply(merge, function(x) estrada.class(x >= .05)))
    if (!is.null(thresh))
      est$opt.index <- max(which.max(est$summary >= thresh)[1] - 1, 1)
    else
      est$opt.index <- 0

    est$criterion <- "estrada.stability"
    return(est)
}

#' @keywords internal
nc.stability <- function(premerge, thresh, rep.num, p, nlams) {
    est <- list()
    est$merge <- sapply(premerge, function(x) vapply(x, function(g)
                    if (is.null(g)) NA_real_ else natural.connectivity(g), numeric(1)))
    est$merge <- matrix(est$merge, ncol=length(premerge))
    est$summary   <- colMeans(est$merge, na.rm=TRUE)
    est$criterion <- "nc.stability"
    est$n.failed  <- colSums(is.na(est$merge))
    return(est)
}

#' @importFrom stats dist
#' @keywords internal
gcd.stability <- function(premerge, thresh, rep.num, p, nlams, merge=NULL,
                          n.failed=NULL) {
    est <- list()
    if (is.null(merge)) {
        gcdpremerge <- lapply(premerge, function(pm) lapply(pm, function(g)
                              if (is.null(g)) NULL else gcvec(g)))
        tmp         <- .gcdmerge(gcdpremerge)
        est$merge   <- tmp$merge
        n.failed    <- tmp$n.failed
    } else
        est$merge <- merge
    if (is.null(n.failed)) n.failed <- numeric(nlams)

    est$summary <- vector('numeric', nlams)
    ## need at least two successful fits for a pairwise distance
    for (i in 1:nlams) est$summary[i] <-
        if (length(est$merge[[i]])) mean(est$merge[[i]]) else NA
    est$criterion <- "graphlet.stability"
    est$n.failed  <- n.failed
    return(est)
}

## Pairwise graphlet correlation distances over successful fits
##
## @param gcdpremerge list (over lambda) of lists (over subsamples) of graphlet
##        correlation vectors, with NULL for failed fits
#' @keywords internal
.gcdmerge <- function(gcdpremerge) {
    n.failed <- vapply(gcdpremerge, function(x) sum(vapply(x, is.null, logical(1))), numeric(1))
    merge    <- lapply(gcdpremerge, function(x) {
        x <- .dropfailed(x)
        if (length(x)) dist(do.call(rbind, x)) else dist(matrix(0, 0, 0))
    })
    list(merge=merge, n.failed=n.failed)
}

#egraphlet.stability <- function(premerge, thresh, rep.num, p, nlams, norbs=12) {
#    est <- list()
##    estrlist    <- lapply(premerge.reord, function(pm) lapply(pm, estrada))
##    est$merge <- lapply(estrlist, function(estrvec) Reduce("+", estrvec)/rep.num)
#    collect <- lapply(premerge, function(pm) lapply(pm, egraphletlist))
#    collect.reord <- lapply(collect, function(colam) lapply(1:norbs, function(i)
#                      lapply(1:rep.num, function(j) colam[[j]][[i]])))
#    rm(list=c('collect')) ; gc()
#    est$merge <- vector('list', nlams)
#    est$summary <- Matrix(0, length(est$merge), norbs)
#    for (i in 1:nlams) {
#      est$merge[[i]] <- vector('list', norbs)
#      for (j in 1:norbs) {
#        collam <- collect.reord[[i]][[j]]
#        EX2 <- Reduce(.sumsq, collam[-1], init=collam[[1]]^2)/rep.num
#        EX  <- Reduce('+', collam)/rep.num
#        est$merge[[i]][[j]] <- EX2 - EX^2
#        est$summary[i,j] <- 2 * sum(est$merge[[i]][[j]]) / (p * (p - 1))
#      }
#    }

##    if (!is.null(thresh))
##      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
##    else
##      est$opt.index <- 0

#    est$criterion <- "egraphlet.stability"
#    return(est)
#}


#vgraphlet.stability <- function(premerge, thresh, rep.num, p, nlams, norbs=15) {

#    est <- list()
#    collect <- lapply(premerge, function(pm) lapply(pm, vgraphletlist))
##    rm(list=c('collect')) ; gc()
#    est$merge <- vector('list', nlams)
#    est$summary <- matrix(0, nlams, norbs)
#    for (i in 1:nlams) {
##      est$merge[[i]] <- vector('list', norbs)
##      for (j in 1:norbs) {
##        collam <- collect.reord[[i]][[j]]
#        EX2 <- Reduce(.sumsq, collect[[i]][-1], init=collect[[i]][[1]]^2)/rep.num
#        EX  <- Reduce('+', collect[[i]])/rep.num
#        est$merge[[i]] <- EX2 - EX^2
#        est$summary[i,] <- 2*Matrix::colMeans(est$merge[[i]])
##      }
#    }
##    if (!is.null(thresh))
##      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
##    else
##      est$opt.index <- 0
#    est$criterion <- "graphlet.stability"
#    return(est)
#}
