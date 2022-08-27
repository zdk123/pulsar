#' The pulsar package
#'
#' Graphical model selection with the pulsar package
#'
#' @details
#' This package provides methods to select a sparse, undirected graphical model by choosing a penalty parameter (lambda or \eqn{\lambda}) among a list of ordered values of lambda. We use an implementation of the Stability Approach to Regularization Selection (StARS, see references) inspired by the \pkg{huge} package.
#'
#' However, \pkg{pulsar} includes some major differences from other R packages for graphical model estimation and selection (\pkg{glasso}, \pkg{huge}, \pkg{QUIC}, \pkg{XMRF}, \pkg{clime}, \pkg{flare}, etc). The underlying graphical model is computed by passing a function as an argument to \code{\link{pulsar}}. Thus, any algorithm for penalized graphical models can be used in this framework (see \code{\link{pulsar-function}} for more details), including those from the above packages. \pkg{pulsar} brings computational experiments under one roof by separating subsampling and calculation of summary criteria from the user-specified core model. The typical workflow in \pkg{pulsar} is to perform subsampling first (via the \code{\link{pulsar}}) and then refit the model on the full dataset using \code{\link{refit}}.
#'
#' Previous StARS implementations can be inefficient for large graphs or when many subsamples are required. \code{\link{pulsar}} can compute upper and lower bounds on the regularization path for the StARS criterion after only \eqn{2} subsamples which makes it possible to neglect lambda values that are far from the desired StARS regularization parameter, reducing computation time for the rest of the \eqn{N-2} subsamples (Bounded StARS (B-StARS)).
#'
#' We also implement additional subsampling-based graph summary criteria which can be used for more informed model selection. For example, we have shown that induced subgraph (graphlet) stability (G-StARS) improves empirical performance over StARS but other criteria are also offered.
#'
#' Subsampling amounts to running the specified core model for \eqn{N} independent computations. Using the \pkg{batchtools} framework, we provide a simple wrapper, \code{batch.pulsar}, for running \code{\link{pulsar}} in embarrassingly parallel mode in an hpc environment. Summary criteria are computed using a Map/Reduce strategy, which lowers memory footprint for large models.
#' @name pulsar-package
#' @seealso \code{\link{pulsar-function}}, \code{\link{pulsar}}, \code{\link{batch.pulsar}}
#' @docType package
#' @references Müller, C. L., Bonneau, R. A., & Kurtz, Z. D. (2016).Generalized Stability Approach for Regularized Graphical Models.arXiv: https://arxiv.org/abs/1605.07072.
NULL

#' Graphical model functions for pulsar
#'
#' Correctly specify a function for graphical model estimation that is compatible with the pulsar package.
#'
#' @details
#' It is easy to construct your own function for penalized model estimation that can be used with this package. The R function must have correctly specified inputs and outputs and is passed into the \code{fun} argument to \code{\link{pulsar}} or \code{\link{batch.pulsar}}. Any function that does not follow these rules will fail to give the desired output and may trigger an error.
#'
#' These packages on CRAN have functions that work out of the box, so you won't need to construct a wrapper:
#'
#' \tabular{ll}{
#'   ~function~ \tab ~package~\cr
#'    huge     \tab   huge   \cr
#'    sugm     \tab   flare
#' }
#'
#'
#' Inputs:
#'
#' The function may take arbitrary, named arguments but the first argument must be the data \eqn{n*p} data matrix with the \eqn{n} samples in rows and \eqn{p} features in the columns.
#' At least one argument must be named "lambda", which is expected to be a decreasing numeric vector of penalties. The non-data arguments should be passed into \code{\link{pulsar}} or \code{\link{batch.pulsar}} as a named list (the names must match function arguments exactly) to the \code{fargs} argument.
#'
#' Outputs:
#'
#' The output from the function must be a list or another S3 object inherited from a list. At least one member must be named \code{path}. This \code{path} object itself must be a list of \eqn{p*p} adjacency matrices, one for each value of lambda. Each cell in the adjacency matrix contains a 1 or TRUE if there is an edge between two nodes or 0/FALSE otherwise. It is highly recommended (though not enforced by \pkg{pulsar}) that each adjacency matrix be a column-oriented, compressed, sparse matrix from the \pkg{Matrix} package. For example, \code{dgCMatrix}/\code{dsCMatrix} (general/symmetric numeric Matrix) or the 1-bit \code{lgCMatrix}/\code{lsCMatrix} classes.
#' The function may return other named outputs, but these will be ignored.
#'
#' @examples
#' ## Generate a hub example
#'  dat <- huge::huge.generator(100, 40, 'hub', verbose=FALSE)
#'
#' ## Simple correlation thresholding
#' corrthresh <- function(data, lambda) {
#'   S <- cor(data)
#'   path <- lapply(lambda, function(lam) {
#'     tmp <- abs(S) > lam
#'     diag(tmp) <- FALSE
#'     as(tmp, 'lMatrix')
#'   })
#'   list(path=path)
#' }
#'
#' ## Inspect output
#' lam <- getLamPath(getMaxCov(dat$sigmahat), 1e-4, 10)
#' out.cor  <- pulsar(dat$data, corrthresh, fargs=list(lambda=lam))
#' out.cor
#'
#' \dontrun{
#' ## Additional examples
#' ## quic
#' library(QUIC)
#' quicr <- function(data, lambda, ...) {
#'     S    <- cov(data)
#'     est  <- QUIC(S, rho=1, path=lambda, msg=0, tol=1e-2, ...)
#'     est$path <-  lapply(seq(length(lambda)), function(i) {
#'                    ## convert precision array to adj list
#'                    tmp <- est$X[,,i]; diag(tmp) <- 0
#'                  as(tmp!=0, "lMatrix")
#'     })
#'     est
#' }
#' ## clime
#' library(clime)
#' climer <- function(data, lambda, tol=1e-5, ...) {
#'      est <- clime(data, lambda, ...)
#'      est$path <- lapply(est$Omegalist, function(x) {
#'                      diag(x) <- 0
#'                      as(abs(x) > tol, "lMatrix")
#'                  })
#'      est
#' }
#'
#' ## inverse cov shrinkage Schafer and Strimmer, 2005
#' library(corpcor)
#' icovshrink <- function(data, lambda, tol=1e-3, ...) {
#'      path <- lapply(lambda, function(lam) {
#'                      tmp <- invcov.shrink(data, lam, verbose=FALSE)
#'                      diag(tmp) <- 0
#'                      as(abs(tmp) > tol, "lMatrix")
#'                  })
#'      list(path=path)
#' }
#'
#' ## Penalized linear model, only
#' library(glmnet)
#' lasso <- function(data, lambda, respind=1, family="gaussian", ...) {
#'          n <- length(lambda)
#'          tmp <- glmnet(data[,-respind], data[,respind],
#'                                    family=family, lambda=lambda, ...)
#'          path <-lapply(1:n, function(i) as(tmp$beta[,i,drop=FALSE], "lMatrix"))
#'          list(path=path)
#' }
#'
#' ## alternative stability selection (DIFFERENT from hdi package)
#' out <- pulsar(dat$data, lasso, fargs=list(lambda=lam))
#' mergmat <- do.call('cbind', tmp$stars$merge)
#' image(mergmat)
#'}
#' @name pulsar-function
#' @aliases pulsar-function
#' @seealso \code{\link{pulsar}}, \code{\link{batch.pulsar}}, \pkg{huge}, \pkg{Matrix}
#' @references Müller, C. L., Bonneau, R. A., & Kurtz, Z. D. (2016). Generalized Stability Approach for Regularized Graphical Models. arXiv: https://arxiv.org/abs/1605.07072.
NULL
