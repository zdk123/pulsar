## Supporting and generic functions

#' @keywords internal
getArgs <- function(call, envir=parent.frame()) {
    fin    <- lapply(call, eval, envir=envir)
    forms  <- formals(fin[[1]])
    iscall <- sapply(forms, class) == 'call'
    iscall <- iscall & !(names(forms) %in% c('regdir', 'regid'))
    forms[iscall] <- lapply(forms[iscall], eval)
    c(forms[!(names(forms) %in% names(fin))], fin)
}

#' @keywords internal
.pcheck <- function(obj) {
    if (!inherits(obj, 'pulsar'))
        stop("obj must be pulsar output")
}

#' @keywords internal
.critcheck <- function(obj, criterion=NULL) {
    if (!(criterion %in% names(obj)))
        warning('desired criterion was not used in the pulsar run')
}


#' Get or evaluate an optimal index
#'
#' If the optimal index for the lambda path is not already assigned, then use a validated method to
#' select the optimal index of the lambda path for alternate criteria  (i.e. other than StARS).
#'
#' @param obj the pulsar/batch.pulsar object to evaluate
#' @param criterion a character argument for the desired summary criterion
#' @param ... Ignored
#' @details Automated optimal index selection is [currently] only implemented for \code{gcd} (graphlet stability).
#'
#' Criterion:
#' \itemize{
#'  \item gcd: Select the minimum gcd summary score within the lower and upper StARS bounds.
#' }
#' @return index of the lambda path
#' @seealso \code{\link{opt.index}}
#' @export
get.opt.index <- function(obj, criterion="gcd", ...) {
    optind <- opt.index(obj, criterion)
    if (!is.null(optind)) optind
    if (criterion == 'gcd') {
        if (is.null(obj$stars$lb.index) || !obj$stars$lb.index)
            stop('Lower bound needed for gcd metric (run with lb.stars=TRUE)')
        gcdind <- which.min(getElement(obj, criterion)$summary)
        gcdind <- gcdind + obj$stars$ub.index - 1
        names(gcdind) <- criterion
        return(gcdind)
    } else {
        stop("Currently, gcd is the only supported criterion")
    }
}

#' Optimal index
#'
#' Get or set the optimal index of the lambda path, as determined by a given criterion. \code{value} must be a numeric/int.
#'
#' @param obj a pulsar or batch.pulsar object
#' @param criterion a summary statistic criterion for lambda selection. If value is not named, default to gcd.
#' @seealso \code{\link{get.opt.index}}
#' @export
opt.index <- function(obj, criterion='gcd') {
    .pcheck(obj)
    .critcheck(obj, criterion)
    getElement(obj, criterion)$opt.index
}

#' @param value Integer index for optimal lambda by criterion
#' @rdname opt.index
#' @export
"opt.index<-" <- function(obj, criterion=names(value), value) {
    .pcheck(obj)
    fin <- getArgs(obj$call, obj$envir)
    .critcheck(obj, criterion)
    if (length(criterion) > 1) stop("Select one criterion")
    if (is.null(criterion)) criterion <- 'gcd'
    if (!is.null(value)) {
      if (!is.numeric(value) || value < 1 || value >= length(fin$fargs$lambda))
        stop('Index value must be positive int within range length of lambda path')
    }
    obj[[ criterion ]]$opt.index <- value
    obj
}

#' Lambda path
#'
#' Generate a lambda path sequence in descending order, equally or log-spaced.
#'
#' @param max numeric, maximum lambda value
#' @param min numeric, minimum lambda value
#' @param len numeric/int, length of lambda path
#' @param log logical, should the lambda path be log-spaced
#' @return numeric vector of lambdas
#' @examples
#' ## Generate the data with huge:
#' library(huge)
#' set.seed(10010)
#' p <- 40 ; n <- 100
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#'
#' ## Theoretical lamda max is the maximum abs value of the empirical covariance matrix
#' maxCov <- getMaxCov(dat$data)
#' lams   <- getLamPath(maxCov, 5e-2*maxCov, len=40)
#'
#' @seealso \code{\link{getMaxCov}}
#' @export
getLamPath <- function(max, min, len, log=FALSE) {
    if (max < min) stop('Did you flip min and max?')
    if (log) { min <- log(min) ; max <- log(max) }
    lams  <- seq(max, min, length.out=len)
    if (log) exp(lams)
    else lams
}

#' Max value of cov
#'
#' Get the maximum [absolute] value of a covariance matrix.
#'
#' @param x A matrix/Matrix of data or covariance
#' @param cov Flag if \code{x} is a covariance matrix, Set False is \code{x} is an nxp data matrix. By default, if \code{x} is symmetric, assume it is a covariance matrix.
#' @param abs Flag to get max absolute value
#' @param diag Flag to include diagonal entries in the max
#' @details This function is useful to determine the theoretical value for lambda_max - for Gaussian data, but may be a useful starting point in the general case as well.
#' @seealso \code{\link{getLamPath}}
#' @export
getMaxCov <- function(x, cov=isSymmetric(x), abs=TRUE, diag=FALSE) {
    if (!cov) x <- cov(x)
    tmp <- Matrix::triu(x, k=if (diag) 0 else 1)
    tmp <- if (abs) abs(tmp) else tmp
    max(tmp)
}
