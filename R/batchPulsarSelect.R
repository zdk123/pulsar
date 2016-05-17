## /TODO: reimplement
batch.pulsar <- function(data, fun=huge::huge, fargs=list(), criterion=c("stars"),
                            thresh = 0.1, subsample.ratio = NULL, lb.stars=FALSE, ub.stars=FALSE,
                            rep.num = 20, regid = "batchtest", regdir="./", init="subtwo",
                            conffile = ".BatchJobs.R", job.res=list())  {
    gcinfo(FALSE)
    n <- nrow(data)
    p <- ncol(data)
    # min requirements for function args
    if (is.null(fargs$lambda)) {
        stop(paste('Error: missing members in fargs:', 
             paste(c('lambda')[c(is.null(fargs$lambda))])))
    }
    nlams <- length(fargs$lambda)

    knowncrits <- c("stars", "gcd")
    if (!all(criterion %in% knowncrits))
       stop(paste('Error: unknown criterion', paste(criterion[!(criterion %in% knowncrits)], collapse=", "), sep=": "))

    if (is.null(subsample.ratio)) {
        if (n > 144)
            subsample.ratio = 10 * sqrt(n)/n
        if (n <= 144)
            subsample.ratio = 0.8
    }
    ind.sample <- replicate(rep.num, sample(c(1:n), floor(n * subsample.ratio),
                    replace = FALSE), simplify=FALSE)


    estFun <- function(ind.sample, fargs, data, fun) {
        tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
        if (is.null(tmp$path)) stop('Error: expected data stucture with \'path\' member') 
        return(tmp$path)
    }


    if (lb.stars) {
       if (!("stars" %in% criterion)) # || length(criterion)>1)
       stop('Lower/Upper bound method must be used with StARS')
       minN <- 2
       regid  <- paste(regid, init, sep="_")
       regdir <- paste(regdir, init, sep="_")
    } else 
        minN <- rep.num

    isamp <- ind.sample[1:minN]
    out <- batchply(data, estFun, fun, fargs, isamp, regid, regdir, conffile, job.res)
    reg <- out$reg ; id  <- out$id


    if (lb.stars) {
        if (is.null(thresh)) {warning("no threshold provided, using th=0.1") ; thresh <- .1}
        doneRun <- waitForJobs(reg, id)
        if (!doneRun) {
            stop('Errors in batch jobs for computing initial stability')
        }

 
        lb.starsmerge <- reduceResults(reg, fun=function(job, res, aggr) 
                         lapply(1:length(aggr), function(i) aggr[[i]] + res[[i]]))
        print(rep.num)
        lb.est <- stars.stability(NULL, thresh, minN, p, lb.starsmerge)

        if (ub.stars) {
            # upper bound is determined by equivilent of MaxEnt of Poisson Binomial
            pmean <- unlist(lapply(lb.est$merge, function(x) 2*sum(Matrix::summary(x)[,3]) / (p*(p-1))))
            ub.summary   <- cummax(4*pmean*(1-pmean))
            ub.index <- max(which.max(ub.summary >= thresh)[1] - 2, 1)
        } else ub.index <- 1

        fargs$lambda <- fargs$lambda[ub.index:lb.est$opt.index]
        nlams <- length(fargs$lambda)

        regid  <- gsub(paste("_", init, sep=""), "", regid)
        regdir <- gsub(paste("_", init, sep=""), "", regdir)
        isamp <- ind.sample[-(1:minN)]
        out <- batchply(data, estFun, fun, fargs, isamp, regid, regdir, conffile, job.res)
        reg <- out$reg ; id  <- out$id
    }

    # jobs w/ no errors
    doneRun <- waitForJobs(reg, id)
    jdone   <- findDone(reg, id)
    rep.num <- length(jdone)
    if (!doneRun) {
        warning('Not all batch jobs completed... proceeding anyway')
    }

    est <- list()
    for (i in 1:length(criterion)) {
      crit <- criterion[i]
      if (crit == "stars") {
         starsmerge <- reduceResults(reg, fun=function(job, res, aggr) 
            lapply(1:length(aggr), function(i) aggr[[i]] + res[[i]]))
         est$stars <- stars.stability(NULL, thresh, rep.num, p, starsmerge)
     }
     
     if (crit == "gcd") {
        gcdpremerge <- reduceResultsList(reg, fun=function(job, res) lapply(res, gcdvec))
        gcdmerge <- lapply(1:nlams, function(i) t(dist(sapply(1:rep.num, 
                           function(j) gcdpremerge[[j]][[i]]))))
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

    class(est) <- "batch.pulsar"
    return(c(est, list(id=id, reg=reg)))
}



batchply <- function(data, estFun, fun, fargs, ind.sample, regid, regdir, conffile, job.res) {
    loadConfig(conffile)
    reg <- makeRegistry(id=regid, file.dir=regdir)
    id  <- batchMap(reg, estFun, ind.sample, more.args = list(fargs=fargs, data=data, fun=fun))
    doneSub <- submitJobs(reg, resources=job.res)
    return(list(reg=reg, id=id))
}
