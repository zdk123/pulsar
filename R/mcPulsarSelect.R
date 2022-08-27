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
                          pass.errors=TRUE, ...) {
  ## capture errors/warnings from mclapply
  warn <- NULL
  env <- environment()
  withCallingHandlers({out <- mclapply(X, FUN, mc.cores=mc.cores, ...)},
                      warning=function(w) {
                         # why doesn't assignment mode work when ncores>1
                         assign('warn', c(warn, w$message), env)
                         invokeRestart("muffleWarning")
                       })
  ## handle errors
  errors <- sapply(out, inherits, what="try-error")
  if (any(errors)) {
    error <- table(trimws(sapply(out[errors], '[', 1)))
    msg   <- paste(sprintf('\n%s job%s failed with: "%s"',
                  error, ifelse(error>1, 's', ''), names(error)), collapse="" )
    if (!pass.errors || all(errors)) {
      stop(msg, call.=FALSE)
    } else {
      warning(msg, call.=FALSE)
      ## continue with successful jobs
      out <- out[!errors]
    }
  } else {
    ## will only invoke if ncore=1, otherwise mclapply suppresses warnings
    if (length(warn) > 0) {
      twarn <- table(trimws(sapply(warn, '[', 1)))
      msg   <- paste(sprintf('%s job%s had warning: "%s"',
                              twarn, ifelse(twarn>1, 's', ''), names(twarn)),
                      collapse="\n")
      warning(msg, call.=FALSE)
    }
  }## no errors detected, continue

  # reset previous warn option
  # options(warn=warn.opt)
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
#' @param ncores number of cores to use for subsampling. See \code{batch.pulsar} for more parallelization options.
#' @param refit Boolean flag to refit on the full dataset after pulsar is run. (see also \code{\link{refit}})
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
#'}
#' @return call: the original function call
#' @details
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
                   ncores = 1, refit=TRUE)  {
  gcinfo(FALSE)
  n <- nrow(data)
  p <- ncol(data)
  # min requirements for function args
  .lamcheck(fargs$lambda)
  nlams <- length(fargs$lambda)
  knowncrits <- c("stars", "diss", "estrada", "gcd", "nc", "sufficiency") #"vgraphlet", "egraphlet",
  .critcheck0(criterion, knowncrits)
  subsample.ratio <- .ratcheck(subsample.ratio, n)

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
      return(tmp$path)
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
      est$stars <- stars.stability(premerge.reord, thresh, rep.num, p)

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

    est$stars$lb.index  <- lb.est$opt.index
    est$stars$ub.index  <- ub.index
    est$stars$opt.index <- est$stars$opt.index + ub.index - 1L
  }

  if ("stars" %in% criterion) {
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
.starsind <- function(summary, thresh, offset=1) {
    max(which.max(summary >= thresh)[1] - offset, 1)
}

#' @keywords internal
stars.stability <- function(premerge, stars.thresh, rep.num, p, merge=NULL) {
    if (is.null(stars.thresh)) stars.thresh <- 0.1
    est <- list()


    # do.call(batchtools::reduceResults,
    #                  c(list(reg=reg, fun=starsaggfun), reduceargs))

    if (is.null(merge)) {
      est$merge <- lapply(premerge, function(x) Reduce("+", x))
      gc() # flush
    } else est$merge <- merge
    est$summary <- rep(0, length(est$merge))

    for (i in 1:length(est$merge)) {
      est$merge[[i]] <- est$merge[[i]]/rep.num
      est$summary[i] <- 4 * sum(est$merge[[i]] * (1 - est$merge[[i]])) / (p * (p - 1))
    }
    ## monotonize variability
    est$summary   <- cummax(est$summary)
    est$opt.index <- .starsind(est$summary, stars.thresh)
    est$criterion <- "stars.stability"
    est$thresh    <- stars.thresh
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
    disslist  <- lapply(premerge, function(pm) lapply(pm, graph.diss))
    est$merge <- lapply(disslist, function(dissmat) Reduce("+", dissmat)/rep.num)
    mergesq   <- lapply(disslist, function(dissmat)
                        Reduce(.sumsq, dissmat[-1], init=dissmat[[1]]^2)/rep.num)

    gc() # flush
    est$summary <- rep(0, length(est$merge))
    for (i in 1:length(est$merge)) {
        tmp <- mergesq[[i]] - est$merge[[i]]^2
        est$summary[i] <- sum(triu(tmp))  / (p * (p - 1))
    }
    est$mergesq <- mergesq
    est$criterion <- "diss.stability"
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
    est$merge <- sapply(premerge, function(x) sapply(x, natural.connectivity))
    est$summary <- colMeans(est$merge)
    est$criterion <- "nc.stability"
    return(est)
}

#' @importFrom stats dist
#' @keywords internal
gcd.stability <- function(premerge, thresh, rep.num, p, nlams, merge=NULL) {
    est <- list()
    if (is.null(merge))
        est$merge <- lapply(premerge, function(pm) dist(t(sapply(pm, gcvec))))
    else
        est$merge <- merge

    est$summary <- vector('numeric', nlams)
    for (i in 1:nlams) est$summary[i] <- mean(est$merge[[i]])
    est$criterion <- "graphlet.stability"
    return(est)
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
