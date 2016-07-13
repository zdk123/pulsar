#' Refit pulsar model
#' 
#' Run the supplied graphical model function on the whole dataset and refit with the selected lambda(s)
#'
#' @param obj a \code{pulsar} or \code{batch.pulsar} object
#' @param criterion a vector of criteria for refitting on full data. An optimal index must be defined for each or a warning will result.
#' @return an S3 \code{pulsar.refit} object with members:
#' \itemize{
#'   \item est: the raw output from the graphical model function (what was used to generate \code{obj}) but applied to the full dataset.
#'   \item refit: a named list of adjacency matrices, for each optimal criterion in \code{obj} or specified in the \code{criterion} argument
#'   \item fun: the original function for estimating the graphical model
#'}
#' @export
refit <- function(obj, criterion) {
    UseMethod("refit")
}

#' @export
refit.pulsar <- function(obj, criterion) {
    .refit.pulsar(obj, criterion)
}

#' @export
refit.batch.pulsar <- function(obj, criterion) {
    .refit.pulsar(obj, criterion)
}

#' @keywords internal
.refit.pulsar <- function(obj, criterion) {
    est <- vector('list', 2)
    names(est) <- c('est', 'refit')
    fin <- getArgs(obj$call, obj$envcl)
    ## call est function on original dataset
    est$est <- do.call(eval(fin$fun), c(fin$fargs, list(fin$data)))
    if (missing(criterion)) criterion <- eval(fin$criterion)
    est$refit <- vector('list', length(criterion))
    names(est$refit) <- criterion
    for (crit in criterion) {
        optind <- obj[[crit]]$opt.index
        if (!is.null(optind)) {
            est$refit[[crit]] <- est$est$path[[optind]]
        } else {
            est$refit[[crit]] <- NULL
            if (crit %in% names(obj)) {
                warning(paste('No optimal index selected for', 
                                crit, 'criterion', sep=" "), call.=FALSE)
            } else
                warning(paste('Unknown criterion', crit, sep=" "), call.=FALSE)
        }
    }
    est$fun <- obj$call$fun
    structure(est, class='pulsar.refit')
}
