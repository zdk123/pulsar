# run stars selection serially for testing purposes
stars <- function(data, fun=huge::huge, fargs=list(),
                            stars.thresh = 0.05, stars.subsample.ratio = NULL, 
                            rep.num = 20, verbose = TRUE, ncores = 1)  {
    gcinfo(FALSE)
    n <- nrow(data)
    p <- ncol(data)
    # min requirements for function args
    if (is.null(fargs$lambda)) {
        stop(paste('Error: missing members in fargs:', 
             paste(c('lambda')[c(is.null(fargs$lambda))])))
    }
    
    if (is.null(stars.subsample.ratio)) {
        if (n > 144)
            stars.subsample.ratio = 10 * sqrt(n)/n
        if (n <= 144)
            stars.subsample.ratio = 0.8
    }
    ind.sample <- replicate(rep.num, sample(c(1:n), floor(n * stars.subsample.ratio), 
                    replace = FALSE), simplify=FALSE)

    estFun <- function(ind.sample) {
        tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
        if (is.null(tmp$path)) stop('Error: expected data stucture with \'path\' member') 
        return(tmp$path)
    }

    premerge <- parallel::mclapply(ind.sample, estFun, mc.cores=ncores)
    est$merge <- Reduce(function(l1, l2) lapply(1:length(l1), 
                    function(i) l1[[i]] + l2[[i]]), premerge, accumulate=FALSE)
    est$merge <- reduceResults(reg, fun=function(job, res, aggr) 
                    lapply(1:length(aggr), function(i) aggr[[i]]+res[[i]])) 
    gc() # flush

    est$variability <- rep(0, length(est$merge))
    for (i in 1:length(est$merge)) {
        est$merge[[i]] <- est$merge[[i]]/rep.num
        # TODO: add normalization constants for non glasso problems
        est$variability[i] <- 4 * sum(est$merge[[i]] * (1 - est$merge[[i]]))/(p * (p - 1))
    }
    est$opt.index    <- max(which.max(est$variability >= stars.thresh)[1] - 1, 1)
#    est$refit        <- est$path[[est$opt.index]]
#    est$opt.lambda   <- est$lambda[est$opt.index]
#    est$opt.sparsity <- est$sparsity[est$opt.index]
    est$criterion <- "stars"
    class(est)    <- "batch.select"
    return(c(est, list(id=id, reg=reg)))
}
