
#' Simulation study of Bayesian adaptive design incorporating data from
#' historical controls: univariate continuous outcome
#'
#' Bayesian adaptive design proposed by Psioda (2018) is implemented, where
#' historical control data is incorporated at interim decision. The univariate
#' continuous outcome is only applicable.
#' @usage
#' iaborrow.cont(
#'   n.HC, n.C1, n.T1, n.C2, n.T2,
#'   out.mean.HC, out.sd.HC, out.mean.C, out.sd.C, out.mean.T, out.sd.T,
#'   cov.mean.HC, cov.sd.HC, cov.cor.HC,
#'   cov.mean.C, cov.sd.C, cov.cor.C, cov.mean.T, cov.sd.T, cov.cor.T,
#'   a0, chains=4, iter=2000, warmup=floor(iter/2), thin=1,
#'   alternative="greater", sig.level=0.05,
#'   w0, accept.t1e=0.1, accept.pow=0.7, nsim=100)
#' @param n.HC Number of patients in external control.
#' @param n.C1 Number of patients in concurrent control at interim analysis
#' occasion.
#' @param n.T1 Number of patients in treatment at interim analysis
#' occasion.
#' @param n.C2 Number of patients in concurrent control at final analysis
#' occasion.
#' @param n.T2 Number of patients in treatment at final analysis
#' occasion.
#' @param out.mean.HC True mean of outcome in historical control.
#' @param out.sd.HC True sd of outcome in historical control.
#' @param out.mean.C True mean of outcome in concurrent control.
#' @param out.sd.C True sd of outcome in concurrent control.
#' @param out.mean.T True mean of outcome in treatment.
#' @param out.sd.T True sd of outcome in treatment.
#' @param cov.mean.HC Vector of true mean of covariates in historical control.
#' @param cov.sd.HC Vector of true sd of covariates in historical control.
#' @param cov.cor.HC True correlation matrix between covariates in historical
#' control.
#' @param cov.mean.C Vector of true mean of covariates in concurrent control.
#' @param cov.sd.C Vector of true sd of covariates in historical control.
#' @param cov.cor.C True correlation matrix between covariates in concurrent
#' control.
#' @param cov.mean.T Vector of true mean of covariates in treatment.
#' @param cov.sd.T Vector of true sd of covariates in treatment.
#' @param cov.cor.T True correlation matrix between covariates in treatment.
#' @param a0 Parameter of power prior.
#' @param chains Number of Markov chains in MCMC sampling. The default value is
#' \code{chains=4}.
#' @param iter Number of iterations for each chain (including warmup) in MCMC
#' sampling. The default value is \code{iter=2000}.
#' @param warmup Number of warmup (aka burnin) iterations per chain in MCMC
#' sampling. The default value is \code{warmup=floor(iter/2)}.
#' @param thin Period for saving samples in MCMC sampling. The default value
#' is \code{thin=1}.
#' @param alternative Alternative hypothesis to be tested ("greater" or "less").
#' The default value is \code{alternative="greater"}.
#' @param sig.level Significance level. The default value is
#' \code{sig.level=0.05}.
#' @param w0 Vector of range of stopping criteria.
#' @param accept.t1e Acceptable maximum type I error rate. The default value is
#' \code{accept.t1e=0.1}.
#' @param accept.pow Acceptable minimum power. The default value is
#' \code{accept.pow=0.7}.
#' @param nsim Number of simulated trials. The default value is \code{nsim=100}.
#' @details The \code{iaborrow.cont} is a function which generates the operating
#' characteristics of Bayesian adaptive design incorporating data from
#' historical controls. This fixed-borrowing adaptive design incorporates
#' subject-level control data from a previously completed clinical trial. The
#' interim analysis evaluates the extent of prior-data conflict based on a
#' similarity measure for new and historical trial data, and then (i) if the
#' conflict is too great, the new trial is continued and the data are analyzed
#' with non-informative prior, (ii) Otherwise, the new trial is stopped early
#' and the data are analyzed with informative prior. The univariate
#' continuous outcome is only applicable.
#' @return
#' \item{w}{Likelihood ratio statistics.}
#' \item{p1}{Probability meeting that the treatment effect is above 0
#' (if \code{alternative="greater"}) or below 0 (if \code{alternative="less"})
#' at interim analysis occasion.}
#' \item{p2}{Probability meeting that the treatment effect is above 0
#' (if \code{alternative="greater"}) or below 0 (if \code{alternative="less"})
#' at final analysis occasion.}
#' \item{sig.rate}{Reject rate.}
#' \item{max.t1e}{Type I error at every \code{w0}}
#' \item{min.pow}{Power at every \code{w0}}
#' \item{stop.rate}{Stoppage probability at every \code{w0}.}
#' \item{w0.opt}{Optimal value of \code{w0}.}
#' \item{t1e}{Type I error rate at optimal value of \code{w0}.}
#' \item{stop.null}{Stoppage probability at optimal value of \code{w0} under
#' null hypothesis.}
#' \item{exp.size.null}{Expected total sample size under null hypothesis.}
#' \item{pow}{Power at optimal value of \code{w0}.}
#' \item{stop.alt}{Stoppage probability at optimal value of \code{w0} under
#' alternative hypothesis.}
#' \item{exp.size.alt}{Expected total sample size under alternative hypothesis.}
#' @references
#' Psioda MA, Soukup M, Ibrahim JG. A practical Bayesian adaptive design
#' incorporating data from historical controls *Statistics in Medicine*
#' 2018; 37:4054-4070.
#' @examples
#' n.HC        <- 200
#' n.C1        <- 50
#' n.T1        <- 50
#' n.C2        <- 100
#' n.T2        <- 100
#' out.mean.HC <- c(-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6)
#' out.sd.HC   <- 1
#' out.mean.C  <- 0
#' out.sd.C    <- 1
#' out.mean.T  <- 0.46
#' out.sd.T    <- 1
#' cov.mean.HC <- c(0,0)
#' cov.sd.HC   <- c(1,1)
#' cov.cor.HC  <- cbind(c(1.0,0.5),c(0.5,1.0))
#' cov.mean.C  <- c(0,0)
#' cov.sd.C    <- c(1,1)
#' cov.cor.C   <- cbind(c(1.0,0.5),c(0.5,1.0))
#' cov.mean.T  <- c(0,0)
#' cov.sd.T    <- c(1,1)
#' cov.cor.T   <- cbind(c(1.0,0.5),c(0.5,1.0))
#' a0          <- 0.5
#' w0          <- seq(0,5,0.1)
#'
#' iaborrow.cont(
#'   n.HC=n.HC, n.C1=n.C1, n.T1=n.T1, n.C2=n.C2, n.T2=n.T2,
#'   out.mean.HC=out.mean.HC, out.sd.HC=out.sd.HC,
#'   out.mean.C=out.mean.C, out.sd.C=out.sd.C,
#'   out.mean.T=out.mean.T, out.sd.T=out.sd.T,
#'   cov.mean.HC=cov.mean.HC, cov.sd.HC=cov.sd.HC, cov.cor.HC=cov.cor.HC,
#'   cov.mean.C=cov.mean.C, cov.sd.C=cov.sd.C, cov.cor.C=cov.cor.C,
#'   cov.mean.T=cov.mean.T, cov.sd.T=cov.sd.T, cov.cor.T=cov.cor.T,
#'   a0=a0, w0=w0)
#' @import rstan MBESS MASS
#' @export

iaborrow.cont <- function(
  n.HC, n.C1, n.T1, n.C2, n.T2,
  out.mean.HC, out.sd.HC, out.mean.C, out.sd.C, out.mean.T, out.sd.T,
  cov.mean.HC, cov.sd.HC, cov.cor.HC,
  cov.mean.C, cov.sd.C, cov.cor.C, cov.mean.T, cov.sd.T, cov.cor.T,
  a0, chains=4, iter=2000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.05,
  w0, accept.t1e=0.1, accept.pow=0.7, nsim=100)
{
  norm.reg1 <- norm.reg.conc()
  norm.reg2 <- norm.reg.ext()
  norm.reg3 <- norm.reg.hyb()

  out.mean.T <- c(out.mean.C,out.mean.T)

  ls1 <- length(out.mean.HC)
  ls2 <- length(out.mean.T)

  cov.vcov.HC <- MBESS::cor2cov(cor.mat=cov.cor.HC,sd=cov.sd.HC)
  cov.vcov.C  <- MBESS::cor2cov(cor.mat=cov.cor.C,sd=cov.sd.C)
  cov.vcov.T  <- MBESS::cor2cov(cor.mat=cov.cor.T,sd=cov.sd.T)

  w  <- array(0,dim=c(nsim,ls1,ls2))
  p1 <- array(0,dim=c(nsim,ls1,ls2))
  p2 <- array(0,dim=c(nsim,ls1,ls2))

  for(ss in 1:nsim){
    for(s1 in 1:ls1){
      for(s2 in 1:ls2){

        out.HC <- rnorm(n.HC,out.mean.HC[s1],out.sd.HC)
        out.C1 <- rnorm(n.C1,out.mean.C,out.sd.C)
        out.T1 <- rnorm(n.T1,out.mean.T[s2],out.sd.T)
        out.C2 <- c(out.C1,rnorm(n.C2-n.C1,out.mean.C,out.sd.C))
        out.T2 <- c(out.T1,rnorm(n.T2-n.T1,out.mean.T[s2],out.sd.T))

        cov.HC <- MASS::mvrnorm(n.HC,cov.mean.HC,cov.vcov.HC)
        cov.C1 <- MASS::mvrnorm(n.C1,cov.mean.C,cov.vcov.C)
        cov.T1 <- MASS::mvrnorm(n.T1,cov.mean.T,cov.vcov.T)
        cov.C2 <- rbind(cov.C1,MASS::mvrnorm(n.C2-n.C1,cov.mean.C,cov.vcov.C))
        cov.T2 <- rbind(cov.T1,MASS::mvrnorm(n.T2-n.T1,cov.mean.T,cov.vcov.T))

        dat1 <- list(
          n = n.C1+n.T1,
          p = ncol(cov.C1),
          x = rbind(cov.C1,cov.T1),
          z = c(rep(0,n.C1),rep(1,n.T1)),
          y = c(out.C1,out.T1))

        dat2 <- list(
          nHC = n.HC,
          pHC = ncol(cov.HC),
          xHC = cov.HC,
          yHC = out.HC,
          a0  = a0)

        dat3 <- list(
          n   = n.C1+n.T1,
          p   = ncol(cov.C1),
          x   = rbind(cov.C1,cov.T1),
          z   = c(rep(0,n.C1),rep(1,n.T1)),
          y   = c(out.C1,out.T1),
          nHC = n.HC,
          xHC = cov.HC,
          yHC = out.HC,
          a0  = a0)

        dat4 <- list(
          n = n.C2+n.T2,
          p = ncol(cov.C2),
          x = rbind(cov.C2,cov.T2),
          z = c(rep(0,n.C2),rep(1,n.T2)),
          y = c(out.C2,out.T2))

        mcmc1 <- rstan::sampling(norm.reg1,dat1,chains=chains,iter=iter,warmup=warmup,thin=thin,show_messages=FALSE)
        mcmc1.sample <- rstan::extract(mcmc1)

        mcmc2 <- rstan::sampling(norm.reg2,dat2,chains=chains,iter=iter,warmup=warmup,thin=thin,show_messages=FALSE)
        mcmc2.sample <- rstan::extract(mcmc2)

        mcmc3 <- rstan::sampling(norm.reg3,dat3,chains=chains,iter=iter,warmup=warmup,thin=thin,show_messages=FALSE)
        mcmc3.sample <- rstan::extract(mcmc3)

        mcmc4 <- rstan::sampling(norm.reg1,dat4,chains=chains,iter=iter,warmup=warmup,thin=thin,show_messages=FALSE)
        mcmc4.sample <- rstan::extract(mcmc4)

        hat.alpha0 <- median(mcmc2.sample$alpha)
        hat.beta0  <- apply(mcmc2.sample$beta,2,median)

        hat.alpha1 <- median(mcmc1.sample$alpha)
        hat.gamma1 <- median(mcmc1.sample$gamma)
        hat.beta1  <- apply(mcmc1.sample$beta,2,median)

        hat.alpha <- median(mcmc3.sample$alpha)
        hat.gamma <- median(mcmc3.sample$gamma)
        hat.beta  <- apply(mcmc3.sample$beta,2,median)

        m1 <- hat.alpha1+dat1$z*hat.gamma1+(dat1$x)%*%hat.beta1
        m2 <- hat.alpha0                  +(dat2$x)%*%hat.beta0
        m3 <- hat.alpha +dat1$z*hat.gamma +(dat1$x)%*%hat.beta
        m4 <- hat.alpha                   +(dat2$x)%*%hat.beta

        sig1 <- median(mcmc1.sample$sigma)
        sig2 <- median(mcmc2.sample$sigma)
        sig3 <- median(mcmc3.sample$sigma)
        sig4 <- median(mcmc3.sample$sigma)

        w[ss,s1,s2] <- (sum(log(dnorm(dat1$y,m1,sig1)))+a0*sum(log(dnorm(dat2$y,m2,sig2)))
                       -sum(log(dnorm(dat1$y,m3,sig3)))-a0*sum(log(dnorm(dat2$y,m4,sig4))))

        if(alternative=="greater"){
          p1[ss,s1,s2] <- mean(mcmc3.sample$gamma>0)
          p2[ss,s1,s2] <- mean(mcmc4.sample$gamma>0)
        }else if(alternative=="less"){
          p1[ss,s1,s2] <- mean(mcmc3.sample$gamma<0)
          p2[ss,s1,s2] <- mean(mcmc4.sample$gamma<0)
        }
  }}}

  r1 <- (p1>(1-sig.level))
  r2 <- (p2>(1-sig.level))

  lw        <- length(w0)
  sig.rate  <- array(0,dim=c(2,lw,ls1))
  max.t1e   <- numeric(lw)
  min.pow   <- numeric(lw)
  stop.rate <- array(0,dim=c(2,lw,ls1))

  for(i in 1:lw){
    sig.rate[1,i,]  <- apply(r1[,,1]*(w[,,1]<=w0[i])+r2[,,1]*(w[,,1]>w0[i]),2,mean)
    sig.rate[2,i,]  <- apply(r1[,,2]*(w[,,2]<=w0[i])+r2[,,2]*(w[,,2]>w0[i]),2,mean)
    max.t1e[i]      <- max(sig.rate[1,i,])
    min.pow[i]      <- min(sig.rate[2,i,])
    stop.rate[1,i,] <- apply(w[,,1]<=w0[i],2,mean)
    stop.rate[2,i,] <- apply(w[,,2]<=w0[i],2,mean)
  }

  w0.opt <- max(w0[(max.t1e<=accept.t1e)&(min.pow>=accept.pow)])
  w0.flg <- which(w0==w0.opt)

  t1e           <- sig.rate[1,w0.flg,]
  stop.null     <- apply(w[,,1]<=w0.opt,2,mean)
  exp.size.null <- stop.null*(n.C1+n.T1)+(1-stop.null)*(n.C2+n.T2)

  pow          <- sig.rate[2,w0.flg,]
  stop.alt     <- apply(w[,,2]<=w0.opt,2,mean)
  exp.size.alt <- stop.alt*(n.C1+n.T1)+(1-stop.alt)*(n.C2+n.T2)

  return(list(w=w,p1=p1,p2=p2,
              sig.rate=sig.rate,max.t1e=max.t1e,min.pow=min.pow,
              stop.rate=stop.rate,w0.opt=w0.opt,t1e=t1e,stop.null=stop.null,
              exp.size.null=exp.size.null,pow=pow,stop.alt=stop.alt,
              exp.size.alt=exp.size.alt))
}



