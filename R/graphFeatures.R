### Supporting functions to compute features of a graph ##
##   represented as an adjacency matrix (can be sparse) ##
##


#' Graph dissimilarity
#'
#' Dissimilarity matrix of a graph is here defined as the number of neighbors shared by any two nodes.
#'
#' @param G a \eqn{p*p} adjacency matrix (dense or sparse) of a graph.
#' @param sim Flag to return Graph similarity instead (1-dissimilarity)
#' @param loops Flag to consider self loops
#'
#' @return a \eqn{p*p} dissimilarity matrix
#' @references Bochkina, N. (2015). Selection of the Regularization Parameter in Graphical Models using a Priori Knowledge of Network Structure, arXiv: 1509.05326.
#' @export
graph.diss <- function(G, sim=FALSE, loops=FALSE) {
  dmat <- GraphDiss2(G)
  dmat[is.na(dmat)] <- 1
  if (!loops) diag(dmat) <- 0
  if (sim) dmat <- 1-dmat
  dmat
}

#' @keywords internal
GraphDiss2 <- function(G) {
    Gprod <- G %*% G
    Gdiag   <- Matrix::diag(Gprod)
    degProd <- Gdiag %*% t(Gdiag)
    1 - (Gprod / sqrt(degProd))
}


#' Natural Connectivity
#'
#' Compute the natural connectivity of a graph
#'
#' @param G a \eqn{p*p} adjacency matrix (dense or sparse) of a graph. Ignored if \code{eig} is given
#' @param eig precomputed list of eigen vals/vectors (output from \code{eigen}). If NULL, compute for \code{G}.
#' @param norm should the natural connectivity score be normalized
#'
#' @details The natural connectivity of a graph is a useful robustness measure of complex networks, corresponding to the average eigenvalue of the adjacency matrix.
#' @return numeric natural connectivity score
#' @references Jun, W., Barahona, M., Yue-Jin, T., & Hong-Zhong, D. (2010). Natural Connectivity of Complex Networks. Chinese Physics Letters, 27(7), 78902. doi:10.1088/0256-307X/27/7/078902
#' @export
natural.connectivity <- function(G, eig=NULL, norm=TRUE) {
  if (is.null(eig)) {
      eig <- eigen(G)
  }
  estrind <- exp(eig$values)
  nc <- log(mean(estrind))
  if (norm) {
    n <- length(estrind)
    nc <- nc / (n - log(n))
  }
  return(nc)
}


#' @keywords internal
.adj2elist <- function(G) {
    if (inherits(G, "sparseMatrix")) {
        G <- Matrix::triu(G, k=1)
        return(Matrix::summary(G)[,-3])
    } else {
        p <- ncol(G)
        return(arrayInd(which(as.logical(triu(G))), c(p,p)))
    }
}

#' Graphlet correlation vector
#'
#' Compute graphlet correlations over the desired orbits (default is 11 non-redundant orbits of graphlets of size <=4) for a single graph \code{G}
#'
#' @param G a \eqn{p*p} adjacency matrix (dense or sparse) of a graph.
#' @param orbind index vector for which orbits to use for computing pairwise graphlet correlations. Default is from Yaveroğlu et al, 2014 (see References), but 1 offset needed for R-style indexing.
#'
#' @references Hočevar, T., & Demšar, J. (2014). A combinatorial approach to graphlet counting. Bioinformatics (Oxford, England), 30(4), 559–65. doi:10.1093/bioinformatics/btt717
#' @references Yaveroğlu, Ö. N., Malod-Dognin, N., Davis, D., Levnajic, Z., Janjic, V., Karapandza, R., … Pržulj, N. (2014). Revealing the hidden language of complex networks. Scientific Reports, 4, 4547. doi:10.1038/srep04547
#' @importFrom stats cor
#' @importFrom methods as
#' @export
gcvec <- function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1) {
  if (length(orbind) < 2) stop("Only one orbit selected, need at least two to calculate graphlet correlations")
  if (any(orbind > 15))   stop("Only 15 orbits, from 4-node graphlets, can be selected")
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
      return(rep(0, n*(n-1)/2))
  }

  p <- ncol(G)
  gcount <- orca::count4(Elist)
  ## expand missing nodes
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
# deprecate direct call to count4 for CRAN submission
#  gcount <- .C("count4", Elist, dim(Elist),
#      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  gcor <- suppressWarnings(cor(rbind(gcount[,orbind],1), method='spearman'))
  gcor[upper.tri(gcor)]
}

#' @keywords internal
subgraph.centrality <- function(Graph, eigs=NULL, rmdiag=FALSE) {
## Code from estrada, for undirected graph represented by adjacency matrix M
## also return odd/even contributions and the first eigen vector
  if (rmdiag) diag(Graph) <- 0
# Calculate the subgraph centrality
#  if (inherits(Graph, 'Matrix'))
#    eigs <- eigs_sym(Graph, ncol(Graph)-1)
#  else
  if (is.null(eigs))
    eigs <- eigen(Graph)
  l <- eigs$value
  v <- eigs$vector
  v2 <- v^2
  dl  <- l
  edl <- exp(dl)
  fb  <- v2 %*% edl #A vector of sugraph centralities

# Partition into odd and even contributions of the subgraph centrality
  sinhl <- sinh(dl)
  fbodd <- v2 %*% sinhl

  coshl <- cosh(dl)
  feven <- v2 %*% coshl
  out <- list(central=fb, odd=fbodd, even=feven, evec=v, evals=l)
  class(out) <- 'subgraph.centrality'
  return(out)
}

#' @keywords internal
.SMA <- function(x) ((mean(abs(x))))


#' Estrada class
#'
#' Estrada proposes that graphs can be classified into four different classes. We call this the Estrada class.
#' These are:
#'   I. Expander-like
#'  II. Cluster
#' III. Core-Periphery
#' IV.  Mixed.
#' @param G a \eqn{p*p} adjacency matrix of a Graph
#' @param evthresh tolerance for a zero eigenvalue
#' @return Estrada class (\eqn{1-4})
#' @references Estrada, E. (2007). Topological structural classes of complex networks. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 75(1), 1-12. doi:10.1103/PhysRevE.75.016103
#' @export
estrada.class <- function(G, evthresh=1e-3) {
  if (!inherits(G, "subgraph.centrality"))
      G <- subgraph.centrality(G)

  ev1  <- G$evec[,1]
  eval <- G$evals[1]

  if (length(unique(sign(ev1))) == 1 && G$evals[2] > 0) { ## try the second eigen vector
      ev1  <- G$evec[,2]
      eval <- G$evals[2]
  }
  subgodd <- G$odd
  Evratio <- pmax(ev1^2 * sinh(eval) / subgodd, evthresh)
  Evratio[is.nan(Evratio)] <- 0
  if (sum(Evratio==evthresh) > (2/3)*length(Evratio)) return(0)
  delLogEv1 <- log10(sqrt(Evratio))
  delSplit <- split(delLogEv1, sign(delLogEv1))
  Devs  <- lapply(delSplit, .SMA)
  if (is.null(Devs$`-1`) && is.null(Devs$`1`)) return(0)
  if (is.null(Devs$`-1`))  Devs$`-1` <- 0
  else if (is.null(Devs$`1`))   Devs$`1` <-  0

  if (length(Devs) != 2) return(0) ## will we ever get here?
  if (log10(Devs$`1`+ 1e-3)  > -2.1) eclass <- c(3,4)    else eclass <- c(1,2)
  if (log10(Devs$`-1`+1e-3) > -2.1)  eclass <- eclass[2] else eclass <- eclass[1]
  return(eclass)
}
