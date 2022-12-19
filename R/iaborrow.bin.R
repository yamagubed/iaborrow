
#' Conducting simulation study of Bayesian adaptive design incorporating data
#' from historical controls: binary outcome
#'
#' Bayesian adaptive design proposed by Psioda (2018) is implemented, where
#' historical control data is incorporated at interim decision. The binary
#' outcome is applicable.
#' @usage
#' iaborrow.bin(
#'   n.CT, n.CC, n.EC, n.CT1, n.CC1,
#'   out.prob.CT, out.prob.CC, driftOR,
#'   cov.CT, cov.CC, cov.EC, cormat,
#'   chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   a0, alternative="greater", sig.level=0.025, nsim=10)
#' @param n.CT Number of patients in concurrent treatment at final analysis
#' opportunity.
#' @param n.CC Number of patients in concurrent control at final analysis
#' opportunity.
#' @param n.EC Number of patients in external control.
#' @param n.CT1 Number of patients in concurrent treatment at interim analysis
#' opportunity.
#' @param n.CC1 Number of patients in concurrent control at interim analysis
#' opportunity.
#' @param out.prob.CT True probability of outcome in concurrent treatment.
#' @param out.prob.CC True probability of outcome in concurrent control.
#' @param driftOR Odds ratio between external control and concurrent control
#' for which the bias should be plotted.
#' @param cov.CT List of covariates in concurrent treatment. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.CC List of covariates in concurrent control. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.EC List of covariates in historical control. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cormat Matrix of correlation coefficients for outcome and covariates,
#' specified as Gaussian copula parameter.
#' @param a0 Parameter of power prior.
#' @param chains Number of Markov chains in MCMC sampling. The default value is
#' \code{chains=2}.
#' @param iter Number of iterations for each chain (including warmup) in MCMC
#' sampling. The default value is \code{iter=4000}.
#' @param warmup Number of warmup (aka burnin) iterations per chain in MCMC
#' sampling. The default value is \code{warmup=floor(iter/2)}.
#' @param thin Period for saving samples in MCMC sampling. The default value
#' is \code{thin=1}.
#' @param alternative Alternative hypothesis to be tested ("greater" or "less").
#' The default value is \code{alternative="greater"}.
#' @param sig.level Significance level. The default value is
#' \code{sig.level=0.025}.
#' @param nsim Number of simulated trials. The default value is \code{nsim=10}.
#' @return
#' \item{w}{Likelihood ratio statistics.}
#' \item{p1}{Probability meeting that the treatment effect is above 0
#' (if \code{alternative="greater"}) or below 0 (if \code{alternative="less"})
#' at interim analysis opportunity.}
#' \item{p2}{Probability meeting that the treatment effect is above 0
#' (if \code{alternative="greater"}) or below 0 (if \code{alternative="less"})
#' at final analysis opportunity.}
#' \item{r1}{\code{TRUE} when significant at interim analysis opportunity;
#' otherwise \code{FALSE}.}
#' \item{r2}{\code{TRUE} when significant at final analysis opportunity;
#' otherwise \code{FALSE}.}
#' @references
#' Psioda MA, Soukup M, Ibrahim JG. A practical Bayesian adaptive design
#' incorporating data from historical controls *Statistics in Medicine*
#' 2018; 37:4054-4070.
#' @examples
#' n.CT  <- 100
#' n.CC  <- 100
#' n.EC  <- 100
#' n.CT1 <- 50
#' n.CC1 <- 50
#'
#' out.prob.CT <- c(0.23,0.35)
#' out.prob.CC <- 0.23
#' driftOR     <- c(0.6,1.0,1.4)
#'
#' cov.CT <- list(list(dist="norm", mean=0,sd=1),
#'                list(dist="binom",prob=0.4))
#'
#' cov.CC <- list(list(dist="norm", mean=0,sd=1),
#'                list(dist="binom",prob=0.4))
#'
#' cov.EC <- list(list(dist="norm", mean=0,sd=1),
#'                list(dist="binom",prob=0.4))
#'
#' cormat <- rbind(c(  1,0.1,0.1),
#'                 c(0.1,  1,0.1),
#'                 c(0.1,0.1,  1))
#'
#' a0 <- 0.5
#'
#' iaborrow.bin(
#'   n.CT=n.CT, n.CC=n.CC, n.EC=n.EC, n.CT1=n.CT1, n.CC1=n.CC1,
#'   out.prob.CT=out.prob.CT, out.prob.CC=out.prob.CC, driftOR=driftOR,
#'   cov.CT=cov.CT, cov.CC=cov.CC, cov.EC=cov.EC, cormat=cormat,
#'   a0=a0)
#' @import rstan boot
#' @export

iaborrow.bin <- function(
  n.CT, n.CC, n.EC, n.CT1, n.CC1,
  out.prob.CT, out.prob.CC, driftOR,
  cov.CT, cov.CC, cov.EC, cormat,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  a0, alternative="greater", sig.level=0.025, nsim=10)
{
  n.CT2 <- n.CT-n.CT1
  n.CC2 <- n.CC-n.CC1
  ncov  <- length(cov.CT)

  out.prob.EC <- (driftOR*out.prob.CC/(1-out.prob.CC))/(1+driftOR*out.prob.CC/(1-out.prob.CC))

  ls1 <- length(out.prob.EC)
  ls2 <- length(out.prob.CT)

  w  <- array(0,dim=c(nsim,ls1,ls2))
  p1 <- array(0,dim=c(nsim,ls1,ls2))
  p2 <- array(0,dim=c(nsim,ls1,ls2))

  for(ss in 1:nsim){
    for(s1 in 1:ls1){
      for(s2 in 1:ls2){

        marg.CT <- list(list(dist="binom",parm=list(size=1,prob=out.prob.CT[s2])))
        marg.CC <- list(list(dist="binom",parm=list(size=1,prob=out.prob.CC)))
        marg.EC <- list(list(dist="binom",parm=list(size=1,prob=out.prob.EC[s1])))

        for(i in 1:ncov){
          if(cov.CT[[i]]$dist=="norm"){
            marg.CT <- append(marg.CT,list(list(dist=cov.CT[[i]]$dist,parm=list(mean=cov.CT[[i]]$mean,sd=cov.CT[[i]]$sd))))
            marg.CC <- append(marg.CC,list(list(dist=cov.CC[[i]]$dist,parm=list(mean=cov.CC[[i]]$mean,sd=cov.CC[[i]]$sd))))
            marg.EC <- append(marg.EC,list(list(dist=cov.EC[[i]]$dist,parm=list(mean=cov.EC[[i]]$mean,sd=cov.EC[[i]]$sd))))
          }else if(cov.CT[[i]]$dist=="binom"){
            marg.CT <- append(marg.CT,list(list(dist=cov.CT[[i]]$dist,parm=list(size=1,prob=cov.CT[[i]]$prob))))
            marg.CC <- append(marg.CC,list(list(dist=cov.CC[[i]]$dist,parm=list(size=1,prob=cov.CC[[i]]$prob))))
            marg.EC <- append(marg.EC,list(list(dist=cov.EC[[i]]$dist,parm=list(size=1,prob=cov.EC[[i]]$prob))))
          }
        }

        cvec <- cormat[lower.tri(cormat)]
        data.CT1 <- datagen(margdist=marg.CT,corvec=cvec,nsim=n.CT1)
        data.CC1 <- datagen(margdist=marg.CC,corvec=cvec,nsim=n.CC1)
        data.CT2 <- datagen(margdist=marg.CT,corvec=cvec,nsim=n.CT2)
        data.CC2 <- datagen(margdist=marg.CC,corvec=cvec,nsim=n.CC2)

        data.CT <- rbind(data.CT1,data.CT2)
        data.CC <- rbind(data.CC1,data.CC2)
        data.EC <- datagen(margdist=marg.EC,corvec=cvec,nsim=n.EC)

        dat1 <- list(
          nCT = n.CT1,
          nCC = n.CC1,
          p   = ncov,
          yCT = data.CT1[,1],
          yCC = data.CC1[,1],
          xCT = data.CT1[,-1],
          xCC = data.CC1[,-1])

        dat2 <- list(
          nEC = n.EC,
          p   = ncov,
          yEC = data.EC[,1],
          xEC = data.EC[,-1],
          a0  = a0)

        dat3 <- list(
          nCT = n.CT1,
          nCC = n.CC1,
          nEC = n.EC,
          p   = ncov,
          yCT = data.CT1[,1],
          yCC = data.CC1[,1],
          yEC = data.EC[,1],
          xCT = data.CT1[,-1],
          xCC = data.CC1[,-1],
          xEC = data.EC[,-1],
          a0  = a0)

        dat4 <- list(
          nCT = n.CT,
          nCC = n.CC,
          p   = ncov,
          yCT = data.CT[,1],
          yCC = data.CC[,1],
          xCT = data.CT[,-1],
          xCC = data.CC[,-1])

        mcmc1 <- rstan::sampling(stanmodels$BinConcurrent,
                                 data          = dat1,
                                 chains        = chains,
                                 iter          = iter,
                                 warmup        = warmup,
                                 thin          = thin,
                                 show_messages = FALSE,
                                 cores         = 1,
                                 refresh       = 0)
        mcmc1.sample <- rstan::extract(mcmc1)

        mcmc2 <- rstan::sampling(stanmodels$BinExternal,
                                 data          = dat2,
                                 chains        = chains,
                                 iter          = iter,
                                 warmup        = warmup,
                                 thin          = thin,
                                 show_messages = FALSE,
                                 cores         = 1,
                                 refresh       = 0)
        mcmc2.sample <- rstan::extract(mcmc2)

        mcmc3 <- rstan::sampling(stanmodels$BinHybrid,
                                 data          = dat3,
                                 chains        = chains,
                                 iter          = iter,
                                 warmup        = warmup,
                                 thin          = thin,
                                 show_messages = FALSE,
                                 cores         = 1,
                                 refresh       = 0)
        mcmc3.sample <- rstan::extract(mcmc3)

        mcmc4 <- rstan::sampling(stanmodels$BinConcurrent,
                                 data          = dat4,
                                 chains        = chains,
                                 iter          = iter,
                                 warmup        = warmup,
                                 thin          = thin,
                                 show_messages = FALSE,
                                 cores         = 1,
                                 refresh       = 0)
        mcmc4.sample <- rstan::extract(mcmc4)

        hat.gammaCC0 <- median(mcmc2.sample$gammaCC)
        hat.beta0    <- apply(mcmc2.sample$beta,2,median)

        hat.gammaCC1 <- median(mcmc1.sample$gammaCC)
        hat.theta1   <- median(mcmc1.sample$theta)
        hat.beta1    <- apply(mcmc1.sample$beta,2,median)

        hat.gammaCC <- median(mcmc3.sample$gammaCC)
        hat.theta   <- median(mcmc3.sample$theta)
        hat.beta    <- apply(mcmc3.sample$beta,2,median)

        m1.CT <- hat.gammaCC1+hat.theta1+(dat1$xCT)%*%hat.beta1
        m1.CC <- hat.gammaCC1           +(dat1$xCC)%*%hat.beta1
        m1    <- rbind(m1.CT,m1.CC)

        m2    <- hat.gammaCC0           +(dat2$xEC)%*%hat.beta0

        m3.CT <- hat.gammaCC +hat.theta +(dat1$xCT)%*%hat.beta
        m3.CC <- hat.gammaCC            +(dat1$xCC)%*%hat.beta
        m3    <- rbind(m3.CT,m3.CC)

        m4    <- hat.gammaCC            +(dat2$xEC)%*%hat.beta

        w[ss,s1,s2] <- (   sum(log(dbinom(c(dat1$yCT,dat1$yCC),1,boot::inv.logit(m1))))
                       +a0*sum(log(dbinom(dat2$yEC,            1,boot::inv.logit(m2))))
                       -   sum(log(dbinom(c(dat1$yCT,dat1$yCC),1,boot::inv.logit(m3))))
                       -a0*sum(log(dbinom(dat2$yEC,            1,boot::inv.logit(m4)))))

        if(alternative=="greater"){
          p1[ss,s1,s2] <- mean(mcmc3.sample$theta>0)
          p2[ss,s1,s2] <- mean(mcmc4.sample$theta>0)
        }else if(alternative=="less"){
          p1[ss,s1,s2] <- mean(mcmc3.sample$theta<0)
          p2[ss,s1,s2] <- mean(mcmc4.sample$theta<0)
        }
  }}}

  r1 <- (p1>(1-sig.level))
  r2 <- (p2>(1-sig.level))

  return(list(w=w,p1=p1,p2=p2,r1=r1,r2=r2))
}
