#' pulsar: batch mode
#'
#' Run pulsar using stability selection, or another criteria, to select an undirected graphical model over a lambda-path.
#'
#' @param wkdir set the working directory if different than \code{\link{getwd}}
#' @param regdir directory to store intermediate batch job files. Default will be a tempory directory
#' @param init text string appended to basename of the regdir path to store the batch jobs for the initial StARS variability estimate (ignored if `regdir` is NA)
#' @param conffile path to or string that identifies a \code{\link[batchtools:batchtools-package]{batchtools}} configuration file. This argument is passed directly to the \code{name} argument of the \code{\link[pulsar]{findConfFile}} function. See that help for detailed explanation.
#' @param job.res named list of resources needed for each job (e.g. for PBS submission script). The format and members depends on configuration and template. See examples section for a Torque example
#' @param cleanup Flag for removing batchtools registry files. Recommended FALSE unless you're sure intermediate data shouldn't be saved.
#' @param refit Boolean flag to refit on the full dataset after pulsar is run. (see also \code{\link{refit}})
#' @return an S3 object of class \code{\link{batch.pulsar}} with a named member for each stability criterion/metric. Within each of these are:
#' \itemize{
#'    \item summary: the summary criterion over \code{rep.num} graphs at each value of lambda
#'    \item criterion: the stability metric
#'    \item merge: the raw criterion merged over the \code{rep.num} graphs (constructed from \code{rep.num} subsamples), prior to summarization
#'    \item opt.ind: index (along the path) of optimal lambda selected by the criterion at the desired threshold. Will return \eqn{0} if no optimum is found or \code{NULL} if selection for the criterion is not implemented.
#'   }
#' If \code{stars} is included as a criterion then additional arguments include
#' \itemize{
#'    \item lb.index: the lambda index of the lower bound at \eqn{N=2} samples if \code{lb.stars} flag is set to TRUE
#'    \item ub.index: the lambda index of the upper bound at \eqn{N=2} samples if \code{ub.stars} flag is set to TRUE
#'}
#' @return reg: Registry object. See \code{batchtools::makeRegistry}
#' @return id: Identifier for mapping graph estimation function. See \code{batchtools::batchMap}
#' @return call: the original function call
#' @examples
#' \dontrun{
#' ## Generate the data with huge:
#' library(huge)
#' set.seed(10010)
#' p <- 400 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(.2, .01, len=40)
#' hugeargs  <- list(lambda=lams, verbose=FALSE)
#'
#' ## Run batch.pulsar using snow on 5 cores, and show progress.
#' options(mc.cores=5)
#' options(batchtools.progress=TRUE, batchtools.verbose=FALSE)
#' out <- batch.pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                  rep.num=20, criterion='stars', conffile='snow')

#' ## Run batch.pulsar on a Torque cluster
#' ## Give each job 1gb of memory and a limit of 30 minutes
#' resources <- list(mem="1GB", nodes="1", walltime="00:30:00")
#' out.p <- batch.pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                  rep.num=100, criterion=c('stars', 'gcd'), conffile='torque'
#'                  job.res=resources, regdir=file.path(getwd(), "testtorq"))
#' plot(out.p)

#' ## take a look at the default torque config and template files we just used
#' file.show(findConfFile('torque'))
#' file.show(findTemplateFile('simpletorque'))
#' }
#' @references Müller, C. L., Bonneau, R., & Kurtz, Z. (2016). Generalized Stability Approach for Regularized Graphical Models. arXiv https://arxiv.org/abs/1605.07072
#' @references Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. Proceedings of the Twenty-Third Annual Conference on Neural Information Processing Systems (NIPS).
#' @references Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012). The huge Package for High-dimensional Undirected Graph Estimation in R. The Journal of Machine Learning Research, 13, 1059–1062.
#' @references Michel Lang, Bernd Bischl, Dirk Surmann (2017). batchtools: Tools for R to work on batch systems. The Journal of Open Source Software, 2(10). URL https://doi.org/10.21105/joss.00135.
#' @inheritParams pulsar
#' @importFrom Matrix mean triu
#' @seealso \code{\link{pulsar}} \code{\link{refit}}
#' @export
batch.pulsar <- function(data, fun=huge::huge, fargs=list(),
                    criterion=c("stars"), thresh = 0.1, subsample.ratio = NULL,
                    lb.stars=FALSE, ub.stars=FALSE, rep.num = 20, seed=NULL,
                    wkdir=getwd(), regdir=NA, init="init", conffile='',
                    job.res=list(), cleanup=FALSE, refit=TRUE) {

    if (!requireNamespace('batchtools', quietly=TRUE)) {
      stop("'batchtools' package required to run 'batch.pulsar'")
    }
    gcinfo(FALSE)
    if (!is.na(regdir))
      if (file.exists(regdir)) stop('Registry directory already exists')

    n <- nrow(data)
    p <- ncol(data)
    # min requirements for function args
    knowncrits <- c("stars", "gcd", "estrada", "sufficiency")
    .lamcheck(fargs$lambda)
    .critcheck0(criterion, knowncrits)
    subsample.ratio <- .ratcheck(subsample.ratio, n)
    nlams    <- length(fargs$lambda)
    conffile <- findConfFile(conffile)

    if (!is.null(seed)) set.seed(seed)
    ind.sample <- replicate(rep.num, sample(c(1:n),
                    floor(n*subsample.ratio), replace=FALSE), simplify=FALSE)
    if (refit) {
      tmp <- 1L:n
      attr(tmp, 'full') <- TRUE
      ind.sample <- c(list(tmp), ind.sample)
    }
    if (!is.null(seed)) set.seed(NULL)

    ## build the estimator function that takes the randomized sample index
    ## and the full data
    estFun <- function(ind.sample, fargs, data, fun) {
      tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
      if (!('path' %in% names(tmp)))
        stop('Error: expected data stucture with \'path\' member')

      if (isTRUE(attr(ind.sample, 'full')))
        return(tmp)
      else
        return(tmp$path)
    }


    est <- list()
    reduceargs <- list() ; reduceGCDargs <- list()
    if (lb.stars) {
      if (!("stars" %in% criterion)) # || length(criterion)>1)
        stop('Lower/Upper bound method must be used with StARS')
      minN   <- 2 + refit
      if (!is.na(regdir)) regdir <- paste(regdir, init, sep="_")
    } else
        minN <- rep.num + refit

    isamp <- ind.sample[1:minN]
    out <- batchply(data, estFun, fun, fargs, isamp,
                    wkdir, regdir, conffile, job.res)
    reg <- out$reg ; id  <- out$id
    ## jobs w/ no errors
    doneRun <- batchtools::waitForJobs(reg=reg, id)
    jdone   <- batchtools::findDone(reg=reg, id)
    pulsar.jobs <- intersect((1+refit):minN, jdone$job.id)

    if (refit) {
      fullmodel <- batchtools::loadResult(id=1, reg=reg)
      minN <- minN - 1L
    } else {
      fullmodel <- NULL
    }
    ## Use this for stars crits
    starsaggfun <- function(res, aggr)
      lapply(1:length(aggr), function(i) aggr[[i]] + res[[i]])

    if (lb.stars) {
      est$init.reg <- reg ; est$init.id  <- id

      if (!doneRun)
        stop('Errors in batch jobs for computing initial stability')

      # collect initial results
      lb.starsmerge <- batchtools::reduceResults(reg=reg, ids=pulsar.jobs, fun=starsaggfun)
      lb.est <- stars.stability(NULL, thresh, minN, p, lb.starsmerge)
      gc(FALSE)
      # compute initial gcd if selected
      if ('gcd' %in% criterion) {
        aggfun <- function(job, res) lapply(res, gcvec)
        lb.gcdpremerge <- do.call(batchtools::reduceResultsList,
                           c(list(reg=reg, ids=pulsar.jobs, fun=aggfun), reduceGCDargs))
      }
      if (cleanup) unlink(reg$file.dir, recursive=TRUE)
      if (lb.est$opt.index == 1)
          warning("Accurate lower bound could not be determined with the first 2 subsamples")
      if (ub.stars) {
        # upper bound is determined by equivilent of MaxEnt of Poisson Binomial
        pmean <- sapply(lb.est$merge, function(x) sum(x)/(p*(p-1)))
        ub.summary <- cummax(4*pmean*(1-pmean))
        tmpub      <- .starsind(ub.summary, thresh, 1)
        if (any(ub.summary == 0)) ## adjust upper bound to exclude empty graphs
          ub.index <- max(tmpub, max(which(ub.summary == 0))+1)
        else
          ub.index <- max(tmpub, 1)

      } else ub.index <- 1
      ## select middle of the lambda path
      fargs$lambda <- fargs$lambda[ub.index:lb.est$opt.index]
      nlams        <- length(fargs$lambda)
      reduceargs   <- list(init=lb.starsmerge[ub.index:lb.est$opt.index])

      ## process initial estimate if gcd is selected
      if ('gcd' %in% criterion) {
        reduceGCDargs <- list(init=lapply(lb.gcdpremerge,
                            function(gcdpm) gcdpm[ub.index:lb.est$opt.index]))
      }
      regdir <- gsub(paste("_", init, sep=""), "", regdir)
      ## Run graph estimation on the rest of the subsamples
      isamp <- ind.sample[-(1:minN)]
      out   <- batchply(data, estFun, fun, fargs, isamp,
                        wkdir, regdir, conffile, job.res)
      reg <- out$reg ; id <- out$id
      doneRun <- batchtools::waitForJobs(reg=reg, id)
      jdone   <- batchtools::findDone(reg=reg, id)
      pulsar.jobs <- intersect((1+refit):rep.num, jdone$job.id)
    }
    rep.num <- length(pulsar.jobs) # adjust denominator to successfull jobs
    if (lb.stars) rep.num <- rep.num + minN
    if (!doneRun)
      warning(paste("Only", jdone, "jobs completed... proceeding anyway"))

    for (i in 1:length(criterion)) {
      crit <- criterion[i]
      if (crit == "stars") {
        ## Reduce results, include init estimate from N=2 if lb/ub is used ##
        starsmerge <- do.call(batchtools::reduceResults,
                         c(list(reg=reg, ids=pulsar.jobs, fun=starsaggfun), reduceargs))
        est$stars  <- stars.stability(NULL, thresh, rep.num, p, starsmerge)
      }

      if (crit == "gcd") {
        gcdaggfun   <- function(res) lapply(res, gcvec)
        gcdpremerge <- c(reduceGCDargs$init,
                        batchtools::reduceResultsList(reg=reg, ids=pulsar.jobs, fun=gcdaggfun))
        gcdmerge    <- lapply(1:nlams, function(i) dist(t(sapply(1:rep.num, function(j) gcdpremerge[[j]][[i]]))))
        est$gcd <- gcd.stability(NULL, thresh, rep.num, p, nlams, gcdmerge)
      }

      else if (crit == "estrada") {
        if (!("stars" %in% criterion))
          warning('Need StaRS for computing Estrada classes... not run')
        else
          est$estrada <- estrada.stability(est$stars$merge, thresh, rep.num, p, nlams)
      }

      else if (crit == "sufficiency") {
        if (!("stars" %in% criterion)) warning('Need StaRS for computing sufficiency... not run')
        else  est$sufficiency <- sufficiency(est$stars$merge, rep.num, p, nlams)
      }

    }

    if (lb.stars) {
      ## split indices of init and full stars estimate
      pind <- ub.index:lb.est$opt.index
      pinv <- setdiff(1:length(lb.est$summary), pind)
      ## stitch back together init and full stars summaries
      tmpsumm       <- vector('numeric', length(lb.est$summary))
      tmpsumm[pinv] <- lb.est$summary[pinv]
      tmpsumm[pind] <- est$stars$summary
      est$stars$summary <- tmpsumm
      ## stitch back together init and full stars merges
      tmpmerg <- vector('list', length(lb.est$summary))
      tmpmerg[pinv]   <- lb.est$merge[pinv]
      tmpmerg[pind]   <- est$stars$merge
      est$stars$merge <- tmpmerg
      ## record stars-related indices
      est$stars$lb.index  <- lb.est$opt.index
      est$stars$ub.index  <- ub.index
      est$stars$opt.index <- est$stars$opt.index + ub.index - 1
    }

    if (cleanup) unlink(reg$file.dir, recursive=TRUE)
    est$id  <- id
    est$reg <- reg

    if ("stars" %in% criterion) {
      if (est$stars$opt.index == 1) {
        direction <- if (any(est$stars$summary >= .1)) "larger" else "smaller"
        warning(paste("Optimal lambda may be", direction,
                      "than the supplied values"))
      }
    }
    est$call  <- match.call()
    est$est   <- fullmodel
    est$envir <- parent.frame()
    return(structure(est, class = c("batch.pulsar", "pulsar")))
}

#' @keywords internal
.get_batchtools_conffile <- function(conffile) {
  ## try to eval batchtools' default for makeRegistry
  ## otherwise use pulsar's findConfFile
  defconffile <- batchtools::findConfFile()
  if (length(defconffile)==0 || is.null(defconffile))
    defconffile <- findConfFile(conffile)
  defconffile
}

#' @keywords internal
batchply <- function(data, estFun, fun, fargs, ind.sample, wkdir, regdir,
                     conffile, job.res) {
  reg <- batchtools::makeRegistry(file.dir=regdir, work.dir=wkdir,
                                  conf.file=findConfFile(conffile))
  args <- list(fargs=fargs, data=data, fun=fun)
  id   <- batchtools::batchMap(estFun, ind.sample, more.args=args, reg=reg)
  doneSub <- batchtools::submitJobs(reg=reg, resources=job.res)
  return(list(reg=reg, id=id))
}
