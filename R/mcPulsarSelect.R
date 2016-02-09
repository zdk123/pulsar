pulsar <- function(data, fun=huge::huge, fargs=list(), criterion=c("stars"), 
                    thresh = NULL, subsample.ratio = NULL, 
                    rep.num = 20, lb.stars=FALSE, ncores = 1, ...)  {
    gcinfo(FALSE)
    n <- nrow(data)
    p <- ncol(data)
    # min requirements for function args
    if (is.null(fargs$lambda)) {
        stop(paste('Error: missing members in fargs:', 
             paste(c('lambda')[c(is.null(fargs$lambda))])))
    }
    nlams <- length(fargs$lambda)

    knowncrits <- c("stars", "diss", "estrada", "graphlet", "nc", "sufficiency")
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

    estFun <- function(ind.sample) {
        tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
        if (is.null(tmp$path)) stop('Error: expected data stucture with \'path\' member') 
        return(tmp$path)
    }

    if (lb.stars) {
        if (!("stars" %in% criterion)) stop('Lower bound method only available for StARS')
        minN <- 2 # Hard code for now
        lb.premerge <- parallel::mclapply(ind.sample[1:minN], estFun, mc.cores=ncores)
        lb.est <- stability(lb.premerge, thresh, minN, p)
        fargs$lambda <- fargs$lambda[1:lb.est$opt.index]
        lb.premerge <- lapply(lb.premerge, function(ppm) ppm[1:lb.est$opt.index])
        premerge <- c(lb.premerge, parallel::mclapply(ind.sample[-(1:minN)], estFun, mc.cores=ncores))

    } else {
        premerge <- parallel::mclapply(ind.sample, estFun, mc.cores=ncores)
    }

    premerge.reord <- lapply(1:nlams, function(i) lapply(1:rep.num, function(j) premerge[[j]][[i]]))
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
        if (!("stars" %in% criterion)) warning('Need StaRS for computing Estrada classes... not run')
        else  est$estrada <- estrada.stability(est$stars$merge, thresh, rep.num, p, nlams)
      }

      else if (crit == "sufficiency") {
        if (!("stars" %in% criterion)) warning('Need StaRS for computing sufficiency... not run')
        else  est$sufficiency <- sufficiency(est$stars$merge, rep.num, p, nlams)
      }

      else if (crit == "graphlet")
        est$graphlet <- graphlet.stability(premerge.reord, thresh, rep.num, p, nlams)

      else if (crit == "nc")
        est$nc <- nc.stability(premerge.reord, thresh, rep.num, p, nlams)

    }

    if (lb.stars) {
      est$stars$summary <- c(est$stars$summary, lb.est$summary[-(1:lb.est$opt.index)])
      est$stars$merge   <- c(est$stars$merge, lb.est$merge[-(1:lb.est$opt.index)])
      est$stars$lb.opt.index <- lb.est$opt.index
    }

    class(est) <- "pulsar"
    return(est)
}


stars.stability <- function(premerge, stars.thresh, rep.num, p) {
    if (is.null(stars.thresh)) stars.thresh <- 0.05
    est <- list()
    est$merge <- lapply(premerge, function(x) Reduce("+", x))
    gc() # flush
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


sufficiency <- function(merge, rep.num, p) {
## Merge solution from StARS
  est <- list()
  est$merge <- sapply(merge, function(x) apply(x*(1-x), 2, max))
  est$summary <- colMeans(est$merge)
  est$criterion <- 'sufficiency'
  return(est)
}


.sumsq <- function(x,y) x + y^2

diss.stability <- function(premerge, diss.thresh, rep.num, p, nlams) {
    est <- list()
    disslist  <- lapply(premerge, function(pm) lapply(pm, graph.diss))
    est$merge <- lapply(disslist, function(dissmat) Reduce("+", dissmat)/rep.num)
    mergesq   <- lapply(disslist, function(dissmat) Reduce(.sumsq, dissmat)/rep.num)

    gc() # flush
    est$summary <- rep(0, length(est$merge))
    for (i in 1:length(est$merge)) {
        est$merge[[i]] <- mergesq[[i]] - est$merge[[i]]^2
#        est$summary[i] <- 4 * sum(est$merge[[i]] * (1 - est$merge[[i]])) / (p * (p - 1))
        est$summary[i] <- sum(triu(mergesq[[i]])) / (p * (p - 1))
    }
    ## monotonize variability
###    est$summary <- cummax(est$summary)
    
    if (!is.null(diss.thresh))
      est$opt.index <- max(which.max(est$summary >= diss.thresh)[1] - 1, 1)
    else
      est$opt.index <- 0

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


estrada.stability <- function(merge, thresh, rep.num, p, nlams) {
    est <- list()
    est$summary <- unlist(lapply(merge, function(x) estrada.class(x >= .05)))
    ## monotonize variability
#    est$summary <- cummax(est$summary)
    if (!is.null(thresh))
      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
    else
      est$opt.index <- 0

    est$criterion <- "estrada.stability"
    return(est)
}


nc.stability <- function(premerge, thresh, rep.num, p, nlams) {
    est <- list()
    est$merge <- sapply(premerge, function(x) sapply(x, natural.connectivity))
    est$summary <- colMeans(est$merge)
#    return(est$summary)
#    gc() # flush
#    est$summary <- rep(0, length(est$merge))
#    for (i in 1:length(est$merge)) {
#        est$summary[i] <- mean(est$merge[[i]])/rep.num
#    }
    ## monotonize variability
#    est$summary <- cummax(est$summary)
    if (!is.null(thresh))
      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
    else
      est$opt.index <- 0

    est$criterion <- "nc.stability"
    return(est)

}

graphlet.stability <- function(premerge, thresh, rep.num, p, nlams) {

    est <- list()
#    estrlist    <- lapply(premerge.reord, function(pm) lapply(pm, estrada))
#    est$merge <- lapply(estrlist, function(estrvec) Reduce("+", estrvec)/rep.num)
   est$summary <- lapply(premerge, function(x) (((sapply(x, graphletvec)))))
  return(est)
    if (!is.null(thresh))
      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
    else
      est$opt.index <- 0

    est$criterion <- "graphlet.stability"
    return(est)
}
