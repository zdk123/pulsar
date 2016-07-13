#' Generate a string for a temporary directory
#'
#' Generate a string to create a random temporary directory, in a platform indepdendent manner. By default, this directory will live under the subdirectory of the per-session temporary directory given by \code{tempdir} from base R. 
#'
#' @param base the base path for the temporary directory.
#' @param len the number of letters to randomly generate the directory name
#' @param fsep the path separator to use (platform dependent)

#' @details This function creates a random path intended for temporary directories. E.g. for testing pulsar's batch mode. This function is useful if you need a safe place to store (and delete) files without endangering important directories or R's per session tmp directory, given by \code{tempdir}, which may be needed for other uses.

#' @return a character vector representing a file path for a randomly generated directory.
#' @export
getTempDir <- function(base=tempdir(), len=6, fsep=.Platform$file.sep) {
    file.path(base, 
        paste("Rtmp", toupper(paste(sample(letters, len, replace=TRUE), collapse="")), sep=""),
        fsep=fsep)
}


#' pulsar: batch mode
#'
#' Run pulsar stability selection, or other criteria, in batch mode to select the sparsity of an undirected gaussian
#' graphical model. 
#'
#' @param criterion a character vector of selection statistics. Multiple criteria can be supplied. By default, StARS is run. Currently, there are no selection criterion available for summary statistics, aside for stars, so the entire path and summary is returned. The options are:
#' \itemize{
#'    \item stars (Stability approach to regularization selection)
#'    \item gcd   (Graphet correlation distance, requires the \pkg{orca} package)
#' }
#' @param regdir directory to store intermediate batch job files
#' @param regid text string representing a unique registry ID for your project
#' @param init text string appended to basename of regdir to store the batch jobs for the initial StARS variability estimate
#' @param conffile path to BatchJobs configuration file
#' @param job.res named list of resources needed for each job (e.g. for intermediate PBS script). The format and members depends on configuation and template.
#' @param progressbars Flag to show various BatchJobs progress bars. Set FALSE for less output.
#' @param cleanup Flag for removing BatchJob registry files. Recommended FALSE unless you're sure indetermediate data shouldn't be saved.
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
#' @return reg Registry object. See \code{BatchJobs::makeRegistry}
#' @return id Identifier for mapping graph estimation function. See \code{BatchJobs::batchMap}
#' @return call: the original function call
#' @inheritParams pulsar
#' @importFrom Matrix mean triu
#' @export
batch.pulsar <- function(data, fun=huge::huge, fargs=list(), criterion=c("stars"),
                            thresh = 0.1, subsample.ratio = NULL, lb.stars=FALSE, ub.stars=FALSE,
                            rep.num = 20, seed=NULL, regdir=getTempDir(), regid=basename(regdir), 
                            init="subtwo", conffile = ".BatchJobs.R", job.res=list(), 
                            progressbars=TRUE, cleanup=FALSE)  {
    gcinfo(FALSE)
    if (file.exists(regdir)) stop('Registry directory already exists')
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

    knowncrits <- c("stars", "gcd")
    if (!all(criterion %in% knowncrits)) {
       stop(paste('Error: unknown criterion', paste(criterion[!(criterion %in% knowncrits)], 
                   collapse=", "), sep=": "))
    }
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

    estFun <- function(ind.sample, fargs, data, fun) {
        tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
        if (is.null(tmp$path)) stop('Error: expected data stucture with \'path\' member') 
        return(tmp$path)
    }

    est <- list()
    reduceargs <- list() ; reduceGCDargs <- list()
    if (lb.stars) {
       if (!("stars" %in% criterion)) # || length(criterion)>1)
       stop('Lower/Upper bound method must be used with StARS')
       minN <- 2
       regid  <- paste(regid, init, sep="_")
       regdir <- paste(regdir, init, sep="_")
    } else 
        minN <- rep.num

    isamp <- ind.sample[1:minN]
    out <- batchply(data, estFun, fun, fargs, isamp, regid, regdir, conffile, job.res, 
                    progressbar=progressbars)
    reg <- out$reg ; id  <- out$id


    if (lb.stars) {
        est$init.reg <- reg ; est$init.id  <- id

        doneRun <- BatchJobs::waitForJobs(reg, id, progressbar=progressbars)
        if (!doneRun) {
            stop('Errors in batch jobs for computing initial stability')
        }
        # collect initial results
        lb.starsmerge <- BatchJobs::reduceResults(reg, progressbar=progressbars, 
                         fun=function(job, res, aggr) 
                             lapply(1:length(aggr), function(i) aggr[[i]] + res[[i]]))
        lb.est <- stars.stability(NULL, thresh, minN, p, lb.starsmerge)

        # compute initial gcd if selected
        if ('gcd' %in% criterion) {
            aggfun <- function(job, res) lapply(res, gcvec)
            lb.gcdpremerge <- do.call(BatchJobs::reduceResultsList, 
                               c(list(reg=reg, fun=aggfun, progressbar=progressbars),
                                      reduceGCDargs))
        }
        if (cleanup) unlink(regdir, recursive=TRUE)
        if (lb.est$opt.index == 1)
            warning("Accurate lower bound could not be determine with N=2 subsamples")

        if (ub.stars) {
            # upper bound is determined by equivilent of MaxEnt of Poisson Binomial
            pmean <- pmin(sapply(lb.est$merge, function(x) {
                             mean(triu(x, k=1))
                            }), 1)

            ub.summary   <- cummax(4*pmean*(1-pmean))
            ## adjust upper bound to exclude empty graphs
            tmpub <- which.max(ub.summary >= thresh)[1] - 2
            if (any(ub.summary == 0))  ## adjust upper bound to exclude empty graphs
                ub.index <- max(tmpub, max(which(ub.summary == 0))+1)
            else
                ub.index <- max(tmpub, 1)

        } else ub.index <- 1
        # select middle of the lambda path
        fargs$lambda <- fargs$lambda[ub.index:lb.est$opt.index]
        nlams <- length(fargs$lambda)
        reduceargs <- list(init=lb.starsmerge[ub.index:lb.est$opt.index])
        
        # process initial estimate if gcd is selected
        if ('gcd' %in% criterion) {
          reduceGCDargs <- list(init=lapply(lb.gcdpremerge, 
                                     function(gcdpm) gcdpm[ub.index:lb.est$opt.index]))
        }
        regid  <- gsub(paste("_", init, sep=""), "", regid)
        regdir <- gsub(paste("_", init, sep=""), "", regdir)
        # Run graph estimation on the rest of the subsamples
        isamp <- ind.sample[-(1:minN)]
        out <- batchply(data, estFun, fun, fargs, isamp, regid, regdir, conffile, job.res, 
                        progressbar=progressbars)
        reg <- out$reg ; id  <- out$id
    }

    # jobs w/ no errors
    doneRun <- BatchJobs::waitForJobs(reg, id, progressbar=progressbars)
    jdone   <- BatchJobs::findDone(reg, id)
    rep.num <- length(jdone)
    if (lb.stars) rep.num <- rep.num + minN
    if (!doneRun)
        warning(paste("Only", jdone, "jobs completed... proceeding anyway"))

    for (i in 1:length(criterion)) {
      crit <- criterion[i]
      if (crit == "stars") {
         # Reduce results, include init estimate from N=2 if lb/ub is used
         starsaggfun <- function(job, res, aggr) lapply(1:length(aggr), function(i) aggr[[i]] + res[[i]])
         starsmerge <- do.call(BatchJobs::reduceResults, 
                         c(list(reg=reg, fun=starsaggfun, progressbar=progressbars), reduceargs))
         est$stars <- stars.stability(NULL, thresh, rep.num, p, starsmerge)
     }
     
     if (crit == "gcd") {
        gcdaggfun <- function(job, res) lapply(res, gcvec)
        gcdpremerge <- c(reduceGCDargs$init, 
                        BatchJobs::reduceResultsList(reg=reg, fun=gcdaggfun, progressbar=progressbars))
        gcdmerge    <- lapply(1:nlams, function(i) dist(t(sapply(1:rep.num, function(j) 
                                    gcdpremerge[[j]][[i]]))))
        est$gcd <- gcd.stability(NULL, thresh, rep.num, p, nlams, gcdmerge)
     }

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

    if (cleanup) unlink(regdir, recursive=TRUE)
    est$id  <- id
    est$reg <- reg

    if ("stars" %in% criterion) {
        if (est$stars$opt.index == 1)
            warning("Optimal lambda may not be included within the supplied path")
    }
    est$call  <- match.call()
    est$envcl <- parent.frame()
    return(structure(est, class = "batch.pulsar"))
}

#' @keywords internal
batchply <- function(data, estFun, fun, fargs, ind.sample, regid, regdir, 
                     conffile, job.res, progressbar) {
    BatchJobs::loadConfig(conffile)
    reg <- BatchJobs::makeRegistry(id=regid, file.dir=regdir)
    id  <- BatchJobs::batchMap(reg, estFun, ind.sample, 
                more.args = list(fargs=fargs, data=data, fun=fun))
    doneSub <- BatchJobs::submitJobs(reg=reg, resources=job.res, progressbar=progressbar)
    return(list(reg=reg, id=id))
}
