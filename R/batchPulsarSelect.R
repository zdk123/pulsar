#' Generate a string for a temporary directory
#'
#' Generate a string to create a random temporary directory, in a platform indepdendent manner. By default, this directory will live under the subdirectory of the per-session temporary directory given by \code{tempdir} from base R. 
#'
#' @param base the base path for the temporary directory.
#' @param len the number of letters to randomly generate the directory name
#' @param fsep the path separator to use (platform dependent)

#' @details This function creates a random path intended for temporary directories. E.g. for testing pulsar's batch mode. This function is useful if you need a safe place to store (and delete) files without endangering important directories or R's per session tmp directory, given by \code{tempdir}, which may be needed for other uses.

#' @return a character vector representing a file path for a randomly generated directory.
#' @seealso batch.pulsar
#' @export
getTempDir <- function(base=tempdir(), len=6, fsep=.Platform$file.sep) {
    file.path(base, 
        paste("Rtmp", toupper(paste(sample(letters, len, replace=TRUE), collapse="")), 
              sep=""), fsep=fsep)
}


#' pulsar: batch mode
#'
#' Run pulsar using stability selection, or another criteria, to select an undirected graphical model over a lambda-path.
#'
#' @param regdir directory to store intermediate batch job files
#' @param regid text string representing a unique registry ID for your project
#' @param init text string appended to basename of regdir to store the batch jobs for the initial StARS variability estimate
#' @param conffile path to BatchJobs configuration file
#' @param job.res named list of resources needed for each job (e.g. for intermediate PBS script). The format and members depends on configuation and template. See examples section for a Torque example
#' @param progressbars Flag to show various BatchJobs progress bars. Set FALSE for less output.
#' @param cleanup Flag for removing BatchJob registry files. Recommended FALSE unless you're sure indetermediate data shouldn't be saved.
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
#' @return reg: Registry object. See \code{BatchJobs::makeRegistry}
#' @return id: Identifier for mapping graph estimation function. See \code{BatchJobs::batchMap}
#' @return call: the original function call
#' @details 
#' The options for \code{criterion} statistics are:
#' \itemize{
#'    \item stars (Stability approach to regularization selection)
#'    \item gcd   (Graphet correlation distance, requires the \pkg{orca} package)
#'    \item estrada (estrada class) see \code{\link{estrada.class}}
#'    \item sufficiency (Tandon & Ravikumar's sufficiency statistic)
#' }
#' @examples
#' \dontrun{
#' ## Generate the data with huge:
#' library(huge)
#' set.seed(10010)
#' p <- 40 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(.2, .01, len=40)
#'
#' ## Run batch.pulsar on a Torque cluster
#'
#' ## Get example template and config files
#' url <- "https://raw.githubusercontent.com/zdk123/pulsar/master/inst/extdata"
#' download.file(paste(url, "BatchJobsTorque.R", sep="/"), destfile=".BatchJobs.R")
#' download.file(paste(url, "simpletorque.tml", sep="/"),  destfile="simpletorque.tml")
#'
#' ## Give each job 1gb of memory and a limit of 30 minutes
#' resources <- list(memory="1GB", nodes="1", walltime="00:30:00")
#' hugeargs  <- list(lambda=lams, verbose=FALSE)
#' out.p <- batch.pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                       rep.num=100, criterion=c('stars', 'gcd'),
#'                       job.res=resources, regdir=file.path(getwd(), "testtorq"))
#' plot(out.p)
#' }
#' @references Müller, C. L., Bonneau, R., & Kurtz, Z. (2016). Generalized Stability Approach for Regularized Graphical Models. arXiv http://arxiv.org/abs/1605.07072
#' @references Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. Proceedings of the Twenty-Third Annual Conference on Neural Information Processing Systems (NIPS).
#' @references Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012). The huge Package for High-dimensional Undirected Graph Estimation in R. The Journal of Machine Learning Research, 13, 1059–1062.
#' @references Bischl, B., Lang, M., Mersmann, O., Rahnenführer, J., & Weihs, C. (2015). BatchJobs and BatchExperiments : Abstraction Mechanisms for Using R in Batch Environments. Journal of Statistical Software, 64(11), 1–25. doi:10.18637/jss.v064.i11
#' @inheritParams pulsar
#' @importFrom Matrix mean triu
#' @seealso \code{\link{pulsar}}
#' @export
batch.pulsar <- function(data, fun=huge::huge, fargs=list(), criterion=c("stars"), thresh = 0.1,
                         subsample.ratio = NULL, lb.stars=FALSE, ub.stars=FALSE, rep.num = 20,
                         seed=NULL, regdir=getTempDir(), regid=basename(regdir), init="subtwo",
                         conffile = ".BatchJobs.R", job.res=list(), progressbars=TRUE, cleanup=FALSE) {
    gcinfo(FALSE)
    if (file.exists(regdir)) stop('Registry directory already exists')
    n <- nrow(data)
    p <- ncol(data)
    # min requirements for function args
    knowncrits <- c("stars", "gcd", "estrada", "sufficiency")
    .lamcheck(fargs$lambda)
    .critcheck0(criterion, knowncrits)
    subsample.ratio <- .ratcheck(subsample.ratio, n)
    nlams <- length(fargs$lambda)

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
        if (!doneRun)
            stop('Errors in batch jobs for computing initial stability')

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
            warning("Accurate lower bound could not be determined with the first 2 subsamples")
        if (ub.stars) {
            # upper bound is determined by equivilent of MaxEnt of Poisson Binomial
            pmean <- sapply(lb.est$merge, function(x) { sum(x)/(p*(p-1)) })
            ub.summary   <- cummax(4*pmean*(1-pmean))
            tmpub      <- .starsind(ub.summary, thresh, 1)
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
        out   <- batchply(data, estFun, fun, fargs, isamp, regid, regdir, conffile, job.res, 
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

      else if (crit == "estrada") {
        if (!("stars" %in% criterion))
            warning('Need StaRS for computing Estrada classes... not run')
        else  est$estrada <- estrada.stability(est$stars$merge, thresh, rep.num, p, nlams)
      }

      else if (crit == "sufficiency") {
        if (!("stars" %in% criterion)) warning('Need StaRS for computing sufficiency... not run')
        else  est$sufficiency <- sufficiency(est$stars$merge, rep.num, p, nlams)
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
      tmpmerg[p2ind]  <- lb.est$merge[p2ind]
      tmpmerg[pind]   <- est$stars$merge
      est$stars$merge <- tmpmerg

      est$stars$lb.index  <- lb.est$opt.index
      est$stars$ub.index  <- ub.index
      est$stars$opt.index <- est$stars$opt.index + ub.index - 1
    }

    if (cleanup) unlink(regdir, recursive=TRUE)
    est$id  <- id
    est$reg <- reg

    if ("stars" %in% criterion) {
        if (est$stars$opt.index == 1) {
            direction <- if (any(est$stars$summary >= .1)) "larger" else "smaller"
            warning(paste("Optimal lambda may be", direction, "than the  supplied values"))
        }
    }
    est$call  <- match.call()
    est$envir <- parent.frame()
    return(structure(est, class = c("batch.pulsar", "pulsar")))
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
