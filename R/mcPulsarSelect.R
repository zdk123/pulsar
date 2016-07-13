#' pulsar: serial or parallel mode
#'
#' Run pulsar stability selection, or other criteria, to select the sparsity of an undirected gaussian
#' graphical model. 
#'
#' @param data A \eqn{n*p} matrix of data matrix input to solve for the \eqn{p*p} graphical model
#' @param fun pass in a function that returns a list of \eqn{p*p} graphical models along the desired regularization path. The expected inputs to this function are based on the \code{huge} function in the \pkg{huge} package: a raw data input (as opposed to a covariance/correlation matrix) and must return a list or S3 object with a \emph{named} \code{path} list. This should be a list of matrices (or preferable sparse Matrix - see the \pkg{Matrix} package) of adjacency matrices without weights or signs representing the undirected graphical model at each value of \code{lambda}. See details for more information.
#' @param fargs arguments to argument \code{fun}. Must be a named list and requires at least one member \code{lambda}, a numeric vector with values for the penality parameter.
#' @param criterion a character vector of selection statistics. Multiple criteria can be supplied. By default, StARS is run. Currently, there are no selection criterion available for summary statistics, aside for stars, so the entire path and summary is returned. The options are:
#' \itemize{
#'    \item stars (Stability approach to regularization selection)
#'    \item gcd   (Graphet correlation distance, requires the \pkg{orca} package)
#'    \item diss  (Node dissimilarity stability)
#'    \item estrada (estrada class)
#'    \item vgraphlet (vertex orbit/graphlet frequency, requires the \pkg{orca} package)
#'    \item egraphlet (edge orbit/graphlet frequency, requires the \pkg{orca} package)
#'    \item nc  (natural connectivity)
#'    \item sufficiency (Ravikumar's sufficiency statistic)
#' }
#' @param thresh threshold for selection criterion. Only implemented for StARS. thresh = 0.1 is recommended.
#' @param subsample.ratio determine the size of the subsamples. Default is 10*sqrt(n)/n for n > 144 or 0.8 otherwise. Should be strictly less than 1.
#' @param rep.num number of random subsamples to take for graph re-estimation. Default is 20, but more is recommended if using other summary metrics or if using edge frequencies as confidence scores.
#' @param seed A numeric seed to force predictable subsampling. Default is NULL. Use for testing purposes only.
#' @param lb.stars Should the lower bound be computed after N=2 subsamples (should result in considerable speedup and only implemented if stars is selected). If this option is selected, other summary metrics will only be applied to the smaller lambda path.
#' @param ub.stars Should the upper bound be computed after N=2 subsamples (should result in considerable speedup and only implemented if stars is selected). If this option is selected, other summary metrics will only be applied to the smaller lambda path. This option is ignored if the lb.stars flag is FALSE.
#' @param ncores number of cores to use for subsampling. See \code{batch.pulsar} for more paralellization options.
#'
#' @return an S3 object of class \code{pulsar} with a named member for each stability metric run. Within each of these are:
#' \itemize{
#'    \item summary: the summary statistic over \code{rep.num} graphs at each value of lambda
#'    \item criterion: the stability criterion used
#'    \item merge: the raw statistic over the \code{rep.num} graphs, prior to summarization
#'    \item opt.ind: optimal index of lambda selected by the criterion at the desired threshold. Will return \emph{0} if no optimum is found or if selection for the criterion is not implemented.
#'   }
#' If \code{stars} is included as a criterion then additional arguments include
#' \itemize{
#'    \item lb.index the lambda index of the lower bound at N=2 samples if \code{lb.stars} flag is set to TRUE
#'    \item ub.index the lambda index of the upper bound at N=2 samples if \code{ub.stars} flag is set to TRUE
#'}
#' @return call: the original function call
#' @examples
#'
#' ## Generate the data with huge:
#' library(huge)
#' set.seed(10010)
#' p <- 40 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(.2, .01, len=40)
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
#' @useDynLib pulsar
#' @importFrom Matrix mean triu
#' @export
pulsar <- function(data, fun=huge::huge, fargs=list(), criterion=c("stars"), thresh = 0.1,
                    subsample.ratio = NULL, rep.num = 20, seed=NULL,
                     lb.stars=FALSE, ub.stars=FALSE, ncores = 1)  {
    gcinfo(FALSE)
    n <- nrow(data)
    p <- ncol(data)
    # min requirements for function args
    if (is.null(fargs$lambda)) {
        stop(paste('Error: missing members in fargs:', 
             paste(c('lambda')[c(is.null(fargs$lambda))])))
    } else {
        if (!all(fargs$lambda == cummin(fargs$lambda)))
            warning("Are you sure you don't want the lambda path to be monotonically decreasing")
        if (length(fargs$lambda) < 2)
            warning("Only 1 value of lambda is given. Are you sure you want to do model selection?")
    }
    nlams <- length(fargs$lambda)

    knowncrits <- c("stars", "diss", "estrada", "gcd", "nc", "sufficiency") #"vgraphlet", "egraphlet",
    if (!all(criterion %in% knowncrits))
       stop(paste('Error: unknown criterion', 
            paste(criterion[!(criterion %in% knowncrits)], collapse=", "), sep=": "))

    if (is.null(subsample.ratio)) {
        if (n > 144)
            subsample.ratio = 10 * sqrt(n)/n
        if (n <= 144)
            subsample.ratio = 0.8
    }

    if (!is.null(seed)) set.seed(seed)
    ind.sample <- replicate(rep.num, sample(c(1:n), floor(n * subsample.ratio), 
                    replace = FALSE), simplify=FALSE)
    if (!is.null(seed)) set.seed(NULL)
    estFun <- function(ind.sample, fargs) {
        tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
        if (is.null(tmp$path)) stop('Error: expected data stucture with \'path\' member') 
        return(tmp$path)
    }

    if (lb.stars) {
        if (!("stars" %in% criterion)) # || length(criterion)>1)
            stop('Lower/Upper bound method must be used with StARS')
        minN <- 2 # Hard code for now
    } else {
        minN <- rep.num
    }
    isamp <- ind.sample[1:minN]
    premerge <- parallel::mclapply(isamp, estFun, fargs=fargs, mc.cores=ncores)


    if (lb.stars) {
        lb.premerge       <- premerge
        lb.premerge.reord <- lapply(1:nlams, function(i) lapply(1:minN, 
                              function(j) lb.premerge[[j]][[i]]))
        lb.est            <- stars.stability(lb.premerge.reord, thresh, minN, p)

        if (lb.est$opt.index == 1) 
            warning("Accurate lower bound could not be determine with N=2 subsamples")
        if (ub.stars) {
            # upper bound is determined by equivilent of MaxEnt of Poisson Binomial
            pmean <- pmin(sapply(lb.est$merge, function(x) {
                             mean(triu(x, k=1))
                            }), 1)
            ub.summary <- cummax(4*pmean*(1-pmean))
            tmpub      <- which.max(ub.summary >= thresh)[1] - 2
            if (any(ub.summary == 0))  ## adjust upper bound to exclude empty graphs
                ub.index <- max(tmpub, max(which(ub.summary == 0))+1)
            else
                ub.index <- max(tmpub, 1)
        } else ub.index <- 1
        # reselect lambda between bounds
        fargs$lambda <- fargs$lambda[ub.index:lb.est$opt.index]
        nlams <- length(fargs$lambda)
        lb.premerge  <- lapply(lb.premerge, function(ppm) ppm[ub.index:lb.est$opt.index])
        isamp <- ind.sample[-(1:minN)]
        tmp   <- parallel::mclapply(isamp, estFun, fargs=fargs, mc.cores=ncores)
        premerge     <- c(lb.premerge, tmp)
    }

    premerge.reord <- lapply(1:nlams, function(i) lapply(1:rep.num, 
                              function(j) premerge[[j]][[i]]))
    rm(premerge) ; gc()
    est <- list()
    
    for (i in 1:length(criterion)) {
      crit <- criterion[i]

      if (crit == "stars")
        est$stars <- stars.stability(premerge.reord, thresh, rep.num, p)

      else if (crit == "diss")
        est$diss <-  diss.stability(premerge.reord, thresh, rep.num, p, nlams)

#      else if (crit == "estrada")
#        est$estrada <- estrada.stability(premerge.reord, thresh, rep.num, p, nlams)

      else if (crit == "estrada") {
        if (!("stars" %in% criterion)) 
            warning('Need StaRS for computing Estrada classes... not run')
        else  est$estrada <- estrada.stability(est$stars$merge, thresh, rep.num, p, nlams)
      }

      else if (crit == "sufficiency") {
        if (!("stars" %in% criterion)) warning('Need StaRS for computing sufficiency... not run')
        else  est$sufficiency <- sufficiency(est$stars$merge, rep.num, p, nlams)
      }

#      else if (crit == "egraphlet")
#        est$egraphlet <- egraphlet.stability(premerge.reord, thresh, rep.num, p, nlams)

#      else if (crit == "vgraphlet")
#        est$vgraphlet <- vgraphlet.stability(premerge.reord, thresh, rep.num, p, nlams)

      else if (crit == "gcd")
        est$gcd <- gcd.stability(premerge.reord, thresh, rep.num, p, nlams)

      else if (crit == "nc")
        est$nc <- nc.stability(premerge.reord, thresh, rep.num, p, nlams)

    }

    if (lb.stars) {
      find <- 1:length(lb.est$summary)
      pind <- ub.index:lb.est$opt.index
      p2ind <- setdiff(find, pind)
      tmpsumm <- vector('numeric', length(lb.est$summary))
      tmpsumm[p2ind] <- lb.est$summary[p2ind]
      tmpsumm[pind]  <- est$stars$summary
      est$stars$summary <- tmpsumm

      tmpmerg <- vector('list', length(lb.est$summary))
      tmpmerg[p2ind] <- lb.est$merge[p2ind]
      tmpmerg[pind]  <- est$stars$merge
      est$stars$merge <- tmpmerg

      est$stars$lb.index <- lb.est$opt.index
      est$stars$ub.index <- ub.index
      est$stars$opt.index <- est$stars$opt.index + ub.index - 1
    }
    if ("stars" %in% criterion) {
        if (est$stars$opt.index == 1)
            warning("Optimal lambda may not be included within the supplied path")
    }
    est$call  <- match.call()
    est$envir <- parent.frame()
    return(structure(est, class="pulsar"))
}


#' @keywords internal
stars.stability <- function(premerge, stars.thresh, rep.num, p, merge=NULL) {
    if (is.null(stars.thresh)) stars.thresh <- 0.1
    est <- list()

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
    est$summary <- cummax(est$summary)
    est$opt.index    <- max(which.max(est$summary >= stars.thresh)[1] - 1, 1)
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
    if (is.null(merge)) {
        est$merge <- lapply(premerge, function(pm) dist(t(sapply(pm, gcvec))))
    } else {
        est$merge <- merge
    }

    est$summary <- vector('numeric', nlams)
    for (i in 1:nlams) {
        est$summary[i] <- mean(est$merge[[i]])
    }

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
