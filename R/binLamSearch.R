## DEPRICATED as of 02/2016 for lbStARS in the pulsar function

##bin.lam.search <- function(data, lambda, fun=huge::huge, fargs=list(), sel.fun=batch.stars, 
##                            stars.thresh = 0.05, regid="batchtest", regdir="./batchtest", rep0=10, offset0=1,...) {

##    fargs$lambda <- NULL
##    lambda <- sort(lambda, decreasing=TRUE)

##    est <- list()
##    if (as.character(substitute(sel.fun)) != "batch.stars") 
##        stop('Only batch.stars currently supported')
##    
##    args0 <- list(...)
##    args0$rep.num <- rep0
##    est0 <- do.call(batch.stars, c(list(data, fun, fargs=c(fargs, list(lambda=lambda)), 
##                           regid=regid, regdir=paste(regdir, regid, sep="/"), stars.thresh=stars.thresh),
##                           args0))
##    m <- min(est0$opt.index+offset0, length(lambda))
##    i <- 1
##    j <- round(m/2)
##    nit <- 1
##    while (i <= j) {
##        rid <- paste("binsearch", nit, sep="")
##        
##        esttmp <- batch.stars(data, fun, fargs=c(fargs, list(lambda=lambda[i:j])), 
##                           regid=rid, regdir=paste(regdir, regid, rid, sep="/"), stars.thresh=stars.thresh,
##                            ...)
##        esttmp$lambda <- lambda[i:j]
##        print(lambda[i:j])
##        est <- .combine.stars(est, esttmp)
##        # recalculate opt index
##        est$opt.index <- max(which.max(est$variability >= stars.thresh)[1] - 1, 1)
##        if (est$opt.index > 1) break
##        else {
##            i <- j+1
##            j <- round((j+m-1)/2)
##        }
##        
##        nit <- nit+1
##    }
##    est <- .combine.init0(est0, est)
##    if (est$opt.index == 1) warning("Optimal lambda not found")
##    return(est)

##}


##.combine.stars <- function(sel1, sel2) {
##    if (length(sel1) == 0) return(sel2)

##    else {
##        sel1$merge <- c(sel1$merge, sel2$merge)
##        sel1$variability <- c(sel1$variability, sel2$variability)
##        sel1$reg <- c(list(sel1$reg), list(sel2$reg))
##        sel1$lambda <- c(sel1$lambda, sel2$lambda)
##        return(sel1)
##    }
##}


##.combine.init0 <- function(sel0, sel) {
##    n <- length(sel0$variability)
##    m <- length(sel$variability)
##    sel0$variability <- c(sel$variability, tail(sel0$variability, n-m))
##    sel0$merge <- c(sel$merge, tail(sel0$merge, n-m))
##    sel0$reg  <- c(sel$reg, tail(sel0$reg, n-m))
##    sel0$lambda  <- c(sel$lambda, tail(sel0$lambda, n-m))
##    sel0$opt.index <- sel$opt.index
##    sel0
##}
