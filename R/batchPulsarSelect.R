## /TODO: reimplement
batch.pulsar <- function(data, fun=huge::huge, fargs=list(), criterion=c("stars"),
                            stars.thresh = 0.1, stars.subsample.ratio = NULL, 
                            rep.num = 20, verbose = TRUE, regid = "batchtest", regdir="./",
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

    knowncrits <- c("stars")
    if (!all(criterion %in% knowncrits))
       stop(paste('Error: unknown criterion', paste(criterion[!(criterion %in% knowncrits)], collapse=", "), sep=": "))

    if (is.null(stars.subsample.ratio)) {
        if (n > 144)
            stars.subsample.ratio = 10 * sqrt(n)/n
        if (n <= 144)
            stars.subsample.ratio = 0.8
    }
    ind.sample <- replicate(rep.num, sample(c(1:n), floor(n * stars.subsample.ratio),
                    replace = FALSE), simplify=FALSE)

    estFun <- function(ind.sample, fargs) {
        tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
        if (is.null(tmp$path)) stop('Error: expected data stucture with \'path\' member') 
        return(tmp$path)
    }


    if (lb.stars) warning('lower bound rule not yet implemented for batch mode')
         
#        if (!("stars" %in% criterion) || length(criterion)>1) stop('Lower bound method only available for StARS (and currently, only StARS)')
#        minN <- 2 # Hard code for now
#        lb.premerge  <- parallel::mclapply(ind.sample[1:minN], estFun, 
#                      fargs=fargs, mc.cores=ncores)
#        lb.premerge.reord <- lapply(1:nlams, function(i) lapply(1:minN, 
#                              function(j) lb.premerge[[j]][[i]]))
#        lb.est       <- stars.stability(lb.premerge.reord, thresh, minN, p)
#        fargs$lambda <- fargs$lambda[1:lb.est$opt.index]
#        nlams <- length(fargs$lambda)
#        lb.premerge  <- lapply(lb.premerge, function(ppm) ppm[1:lb.est$opt.index])
#        tmp <- parallel::mclapply(ind.sample[-(1:minN)], 
#                      estFun, fargs=fargs, mc.cores=ncores)
#        premerge <- c(lb.premerge, tmp)

#    } else {
        loadConfig(conffile)

        reg <- makeRegistry(id=regid, file.dir=regdir)
        id  <- batchMap(reg, estFun, ind.sample, more.args = list(fargs))
        doneSub <- submitJobs(reg, resources=job.res)
        doneRun <- waitForJobs(reg, id)
                
        # jobs w/ no errors
        if (!doneRun) {
            warning('Error in batch jobs')
        }

#    }


    est <- list()
    est$merge <- reduceResults(reg, fun=function(job, res, aggr) 
                    lapply(1:length(aggr), function(i) aggr[[i]]+res[[i]])) 

    est <- list()

    est$summary <- rep(0, length(est$merge))
    for (i in 1:length(est$merge)) {
        est$merge[[i]] <- est$merge[[i]]/rep.num
        # TODO: add normalization constants for non glasso problems
        est$summary[i] <- 4 * sum(est$merge[[i]] * (1 - est$merge[[i]]))/(p * (p - 1))
    }
    est$opt.index    <- max(which.max(est$summary >= stars.thresh)[1] - 1, 1)
#    est$refit        <- est$path[[est$opt.index]]
#    est$opt.lambda   <- est$lambda[est$opt.index]
#    est$opt.sparsity <- est$sparsity[est$opt.index]
    est$criterion <- "stars"
    class(est)    <- "batch.select"
    return(c(est, list(id=id, reg=reg)))
}
