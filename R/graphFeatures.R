### Supporting functions to compute features of a graph ## 
##   represented as an adjacency matrix (can be sparse) ##
## 

graph.diss <- function(G, sim=FALSE, loops=FALSE) {
  dmat <- GraphDiss(G)
  dmat[is.nan(dmat)] <- 1
  if (!loops) diag(dmat) <- 0
  if (sim) dmat <- 1-dmat
  dmat
}


eigs_sym.dsCMatrix <- function(A, ...) {
    eigs_sym(as(A, 'dgCMatrix'), ...)
}


natural.connectivity <- function(G, eig=NULL, norm=TRUE) {
## eigs -> precomputed eigendecomp

  if (is.null(eig)) {
  
#    arploaded <- tryCatch(library(rARPACK), error=function(e) FALSE)
#    if (!arploaded) {
#      warning('install rARPACK for fast, sparse eigen decomp, proceeding with eigen')
      eig <- eigen(G)
#    } else {
#      eig <- eigs_sym(G, k=ncol(G))
#    }
  }

  estrind <- exp(eig$values)
  nc <- log(mean(estrind))
  
  if (norm) {
    n <- length(estrind)
    nc <- nc / (n - log(n))
  }
  return(nc)
}


graphcorvec <- function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1) {
  ## assume G1, G2 are sparse Matrix objects
  Elist  <-  orca:::convert.graph(Matrix::summary(G)[,-3])
  n <- length(orbind)
  if (ncol(Elist) == 0) return(rep(0, n*(n-1)/2))

  gcount <- .C("count4", Elist, dim(Elist), 
      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits

  gcor <- cor(rbind(gcount[,orbind],1), method='spearman')
  gcor[upper.tri(gcor)]
}


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
  fbodd <- v2 %*% sinhl;

  coshl <- cosh(dl)
  feven <- v2 %*% coshl
  out <- list(central=fb, odd=fbodd, even=feven, evec=v, evals=l)
  class(out) <- 'subgraph.centrality'
  return(out)
}

.LSMA <- function(x) log10(sqrt(mean(abs(x)))+1e-3)

estrada.class <- function(Graph, evthresh=1e-3, oddthresh=1e-3) {

  if (class(Graph) != "subgraph.centrality")
      Graph <- subgraph.centrality(Graph)

  ev1  <- Graph$evec[,1]
  eval <- Graph$evals[1]

#  if (length(unique(sign(ev1))) == 1 && Graph$evals[2] > 0) { ## try the second eigen vector
#      ev1  <- Graph$evec[,2]
#      eval <- Graph$evals[2]
#  }
  subgodd <- Graph$odd
##  keepind <- which(ev1 > evthresh & subgodd > oddthresh)


  ev1[abs(ev1) <= evthresh] <- evthresh
  subgodd[subgodd <= oddthresh] <- oddthresh
#  if ((sum(ev1==evthresh) > 3) | sum(subgodd == oddthresh) > 3) return(0)
  delLogEv1 <- log10(sqrt(ev1^2 * sinh(eval) / subgodd))

  if (any(is.nan(delLogEv1))) print(subgodd)

  delSplit <- split(delLogEv1, sign(delLogEv1))
  logDevs  <- lapply(delSplit, .LSMA)
  if (length(logDevs) != 2) return(0)
  if (logDevs$`1`  > -2.1) eclass <- c(3,4)    else eclass <- c(1,2)
  if (logDevs$`-1` > -2.1) eclass <- eclass[2] else eclass <- eclass[1]
  return(eclass)
}


.nonzero <- function(x, tol=.Machine$double.eps) 
    pmax(x, .Machine$double.eps) > .Machine$double.eps*2


largestCC <- function(G, ncc=1, tol=.Machine$double.eps) {
  n  <- ncol(G)
  if (ncc > n) stop('more components than nodes')
  diag(G) <- 1
  Mn <- matPow(G, n)
  Mnsvd <- svd(G, nu=0, nv=ncc)
  
  keepncc <- which(.nonzero(Mnsvd$d[1:ncc]))
  if (length(keepncc) < ncc) warning("Fewer Connected components than desired")
  
  vects <- Mnsvd$v[,keepncc,drop=FALSE]
  apply(abs(vects), 2, function(x) which(.nonzero(x)))
}



