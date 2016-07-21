#p <- 8
### generate synthetic data
#set.seed(10010)
#dat <- huge::huge.generator(p*10, p, "random", prob=.25, verbose=FALSE, v=.1, u=.5)

#source('pulsarfuns.R')


#context("pulsar: huge, serial mode")
#runtests(pulsar, "pulsar", dat, fun=huge::huge, fargs=list(verbose=FALSE, scr=TRUE), seed=10010)


#context("pulsar: quic, serial mode")
#quicr <- function(data, lambda) {
#    S    <- cov(data)
#    est  <- QUIC::QUIC(S, rho=1, path=lambda, msg=0, tol=1e-2)
#    path <-  lapply(seq(length(lambda)), function(i) {
#                tmp <- est$X[,,i]; diag(tmp) <- 0
#                as(tmp!=0, "lgCMatrix")
#    })
#    est$path <- path
#    est
#}
#runtests(pulsar, "pulsar", dat, fun=quicr, fargs=list(), seed=10010)
