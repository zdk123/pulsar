#' Print a \code{pulsar.refit} S3 object
#'
#' Print information about the model, path length, graph dimension, criterion and optimal indices and graph sparsity.
#'
#' @param x a \code{pulsar.refit}. output from \code{refit}
#' @param ... ignored
#' @importFrom utils capture.output
#' @export
print.pulsar.refit <- function(x, ...) {
    cat("Pulsar-selected refit of", capture.output(print(x$fun)), "\n")
    cat("Path length:", length(x$est$path), "\n")
    cat("Graph dim:  ", ncol(x$est$path[[1]]), "\n")
    crits <- names(x$refit)
    if (length(crits) > 0) {
        critext  <- ifelse(length(crits) > 1, "Criteria:", "Criterion:")
        critext2 <- lapply(crits, function(cr) {
            sp <- sum(x$refit[[cr]]) / ncol(x$refit[[cr]])^2
            optext <- paste(cr, "... sparsity ", signif(sp, 3), sep="")
            paste("  ", optext, sep="")
        })
        cat(critext, "\n", paste(critext2, collapse="\n"), "\n", sep="")
    }
}


#' Print a \code{pulsar} and \code{batch.pulsar} S3 object
#'
#' Print information about the model, path length, graph dimension, criterion and optimal indices, if defined.
#'
#' @param x a fitted \code{pulsar} or \code{batch.pulsar} object
#' @param ... ignored
#' @export
print.pulsar <- function(x, ...) {
    fin <- getArgs(getCall(x), getEnvir(x))
    mode <- ifelse(fin$ncores > 1, "parallel", "serial")
    cat("Mode:", mode)
    .print.pulsar(x, fin)
}

#' @rdname print.pulsar
#' @export
print.batch.pulsar <- function(x, ...) {
    fin <- getArgs(getCall(x), getEnvir(x))
    cat("Mode: batch")
    .print.pulsar(x, fin)
}

#' @keywords internal
.print.pulsar <- function(x, fin) {
    if (fin$lb.stars) {
        cat("... bound index: lower ", x$stars$lb.index,
            ", upper ", x$stars$ub.index, "\n", sep="")
    } else cat("\n")
    cat("Path length:", length(fin$fargs$lambda), "\n")
    cat("Subsamples: ", fin$rep.num, "\n")
    cat("Graph dim:  ", ncol(fin$data), "\n")
    critext  <- ifelse(length(fin$criterion) > 1, "Criteria:", "Criterion:")
    critext2 <- lapply(fin$criterion, function(cr) {
        opt.ind <- x[[cr]]$opt.ind
        optext  <- ifelse(is.null(opt.ind), "",
          paste("... opt: index ", opt.ind, ", lambda ",
          signif(fin$fargs$lambda[opt.ind], 3), sep=""))
        paste("  ", cr, optext, sep="")
        })
    cat(critext, "\n", paste(critext2, collapse="\n"), "\n", sep="")
}

#' Plot a \code{pulsar} S3 object
#'
#' @param x a \code{pulsar} or \code{batch.pulsar} object
#' @param scale Flag to scale non-StARS criterion to max StARS value (or 1)
#' @param invlam Flag to plot 1/lambda
#' @param loglam Flag to plot log[lambda]
#' @param legends Flag to plot legends
#' @param ... ignored
#'
#' @details If both invlam and loglam are given, log[1/lambda] is plotted
#' @export
plot.pulsar <- function(x, scale=TRUE, invlam=FALSE, loglam=FALSE, legends=TRUE, ...) {
    .plot.pulsar(x, scale, invlam, loglam, legends)
}

#' @importFrom graphics plot points legend
#' @keywords internal
.plot.pulsar <- function(x, scale=TRUE, invlam=FALSE, loglam=FALSE, legends=TRUE) {
    fin  <- getArgs(getCall(x), getEnvir(x))
    lams <- fin$fargs$lambda
    xlab <- "lambda"
    if (invlam) {lams <- 1/lams ; xlab <- paste("1/", xlab, sep="")}
    if (loglam) {lams <- log(lams) ; xlab <- paste("log[ ", xlab, " ]", sep="")}

    nlam  <- length(lams)
    crits <- fin$criterion
    n     <- length(crits)
    if (scale) {
        ylab <- "summary (scaled)"
        if ("stars" %in% crits)
            ymax <- max(x$stars$summary)
        else ymax <- 1
    } else {
        ylab <- "summary"
        ymax <- max(unlist(lapply(crits, function(c) x[[ c ]]$summary)))
    }

    yrange <- c(0, ymax)
    plot(lams, seq(yrange[1], yrange[2], length.out=nlam),
         xlab=xlab, ylab=ylab, type='n')
    if (!is.null(x$stars$lb.index)) {
        ilams <- 1:length(lams)
        range1 <- ilams < x$stars$ub.index
        range2 <- ilams > x$stars$lb.index
        range  <- !(range1 | range2)
        ccol   <- vector('numeric', n+1)
        ltys   <- vector('numeric', n+1)
        legs   <- vector('numeric', n+1)
    } else {
        range1 <- rep(FALSE, nlam) ; range2 <- range1
        range  <- !range1
        ccol   <- vector('numeric', n)
        ltys   <- vector('numeric', n)
        legs   <- vector('numeric', n)
    }

    i <- 1 ; lcol <- 1
    optcrits <- c() ; optcols <- c()
    for (cr in crits) {
        summs <- x[[ cr ]]$summary
        optind <- opt.index(x, cr)
        if (scale && cr != "stars") summs <- ymax*summs/max(summs)
        if (length(summs) == nlam) {
            points(lams[range],  summs[range],  type='b', col=lcol)
            points(lams[range1], summs[range1], type='b', col=lcol, lty=2)
            points(lams[range2], summs[range2], type='b', col=lcol, lty=2)
            optind2 <- optind

            if (any(range1 | range2)) {
                ccol[i:(i+1)] <- c(lcol,lcol)
                ltys[i:(i+1)] <- c(2,1)
                legs[i:(i+1)] <- c(paste("b-", cr, sep=""), cr)
                i <- i+1
            } else {
                ccol[i] <- lcol
                ltys[i] <- 1
                legs[i] <- cr
            }
        } else {
            points(lams[range], summs, type='b', col=lcol)
            optind2 <- optind-which(range)[1]+1
            ccol[i] <- lcol
            ltys[i] <- 1
            legs[i] <- cr
        }

        if (!is.null(optind)) {
            points(lams[optind], summs[optind2], type='p', cex=1.5, pch=16, col=lcol)
            optcrits <- c(optcrits, cr)
            optcols  <- c(optcols , lcol)
        }
        lcol <- lcol + 1 ; i <- i + 1
    }

    if (legends) {
        legend('bottomleft', legs, col=ccol, pch=1, lty=ltys, cex=1.4)
        if (length(optcrits) > 0)
          legend('topright', optcrits, pch=16, col=optcols, cex=1.5, title="opt lambda")
    }
}

#' Update a pulsar call
#'
#' Update a pulsar model with new or altered arguments. It does this by extracting the call stored in the object, updating the call and (by default) evaluating it in the environment of the original \code{pulsar} call.
#'
#' @param object a n existing pulsar or batch.pulsar object
#' @param ... arguments to \code{pulsar} to update
#' @param evaluate Flag to evaluate the function. If \code{FALSE}, the updated call is returned without evaluation
#' @details The \code{update} call is evaluated in the environment specified by the \code{pulsar} or \code{batch.pulsar} object, so if any variables were used for arguments to the original call, unless they are purposefully updated, should not be altered. For example, if the variable for the original data is reassigned, the output of \code{update} will not be on the original dataset.
#' @return If \code{evaluate = TRUE}, the fitted object - the same output as \code{pulsar} or \code{batch.pulsar}. Otherwise, the updated call.
#' @examples
#' \dontrun{p <- 40 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(getMaxCov(dat$data), .01, len=20)
#'
#' ## Run pulsar with huge
#' hugeargs <- list(lambda=lams, verbose=FALSE)
#' out.p <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                 rep.num=20, criterion='stars')
#'
#' ## update call, adding bounds
#' out.b <- update(out.p, lb.stars=TRUE, ub.stars=TRUE)
#' }
#' @importFrom stats update
#' @seealso \code{eval}, \code{\link{update}}, \code{\link{pulsar}}, \code{\link{batch.pulsar}}
#' @export
update.pulsar <- function(object, ..., evaluate=TRUE) {
    extras <- match.call(expand.dots=FALSE)$...
    .update.pulsar(object, extras, evaluate)
}

#' @importFrom stats getCall
#' @keywords internal
.update.pulsar <- function(object, extras, evaluate) {
    call <- getCall(object)
    if (is.null(getEnvir(object))) object$envir <- parent.frame()
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
            if (any(!existing)) {
               call <- c(as.list(call), extras[!existing])
               call <- as.call(call)
            }
    }
    if (evaluate)
        eval(call, getEnvir(object))
    else call
}

#' Get calling environment
#'
#' Generic S3 method for extracting an environment from an S3 object. A getter for an explicitly stored environment from an S3 object or list... probably the environment where the original function that created the object was called from. The default method is a wrapper for \code{x$envir}.
#'
#' @param x S3 object to extract the environment
#' @seealso \code{getCall}, \code{environment}, \code{parent.env}, \code{eval}
#' @export
getEnvir <- function(x) {
    UseMethod("getEnvir")
}

#' @rdname getEnvir
#' @export
getEnvir.default <- function(x) {
    getElement(x, "envir")
}
