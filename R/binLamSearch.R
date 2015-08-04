bin.lam.search <- function(data, lambda, fun=huge::huge, fargs=list(), sel.fun=batch.stars, 
                            stars.thresh = 0.05, regid="batchtest", regdir=".", ...) {

    fargs$lambda <- NULL
    lambda <- sort(lambda, decreasing=TRUE)
    m <- length(lambda)
    i <- 1
    j <- round(m/2)
    est <- list()
    if (as.character(substitute(sel.fun)) != "batch.stars") 
        stop('Only batch.stars currently supported')
    
    nit <- 1
    while (i <= j) {
        rid <- paste("binsearch", nit, sep="")
        esttmp <- batch.stars(data, fun, fargs=c(fargs, list(lambda=lambda[i:j])), 
                           regid=rid, regdir=paste(regdir, regid, rid, sep="/"), ...)
        est <- .combine.stars(est, esttmp)
        # recalculate opt index
        est$opt.index <- max(which.max(est$variability >= stars.thresh)[1] - 1, 1)
        if (est$opt.index > 1) return(est)
        else {
            i <- j+1
            j <- round((j+m)/2)
        }
    nit <- nit+1
    }
    if (est$opt.index == 1) warning("Optimal lambda not found")
    return(est)

}


.combine.stars <- function(sel1, sel2) {
    if (length(sel1) == 0) return(sel2)

    else {
        sel1$merge <- c(sel1$merge, sel2$merge)
        sel1$variability <- c(sel1$variability, sel2$variability)
        sel1$reg <- c(list(sel1$reg), list(sel2$reg))
        return(sel1)
    }
}
