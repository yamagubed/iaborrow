
#' Generating multivariate correlated data
#'
#' Multivariate correlated data are generated. Gaussian copula is used to
#' specify the correlation between variables. Any probability distributions
#' available in R STAT is applicable.
#' @usage
#' datagen(margdist,corvec,nsim=1000)
#' @param margdist List of distributions to be used for the data generation.
#' @param corvec Vector of Gaussian copula correlation parameter.
#' @param nsim Number of simulation.
#' @return Data having the specified distribution.
#' @import copula
#' @export

datagen <- function(margdist,corvec,nsim)
{
  varnum <- length(margdist)
  cormat <- copula::normalCopula(param=corvec,dim=varnum,dispstr="un")

  dist <- NULL
  parm <- NULL
  for(i in 1:varnum){
    dist <- c(dist,(margdist[[i]])$dist)
    parm <- append(parm,list((margdist[[i]])$parm))
  }

  mycop <- copula::mvdc(copula=cormat,margins=dist,paramMargins=parm)
  return(copula::rMvdc(nsim,mycop))
}
