
#' Conducting simulation study of Bayesian adaptive design incorporating data
#' from historical controls: time-to-event outcome
#'
#' Bayesian adaptive design proposed by Psioda (2018) is implemented, where
#' historical control data is incorporated at interim decision. The univariate
#' continuous outcome is only applicable.
#' @usage
#' iaborrow.t2e <- function(
#'   n.CT, n.CC, nevent.C, n.EC, nevent.E, nevent.C1, accrual,
#'   out.mevent.CT, out.mevent.CC, driftHR,
#'   cov.CT, cov.CC, cov.EC, cormat,
#'   chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   a0, alternative="greater", sig.level=0.05, nsim=10)
#' @param n.CT Number of patients in concurrent treatment at final analysis
#' occasion.
#' @param n.CC Number of patients in concurrent control at final analysis
#' occasion.
#' @param nevent.C Number of events in concurrent data at final analysis
#' occasion.
#' @param n.EC Number of patients in external control.
#' @param nevent.E Number of events in external control.
#' @param nevent.C1 Number of events in concurrent data at interim analysis
#' occasion.
#' @param accrual Accrual rate.
#' @param out.mevent.CT True median time to event in concurrent treatment.
#' @param out.mevent.CC True median time to event in concurrent control.
#' @param driftHR HR between external control and concurrent control
#' for which the bias should be plotted.
#' @param cov.CT List of covariates in concurrent treatment. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.CC List of covariates in concurrent control. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.EC List of covariates in historical control. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cormat Matrix of correlation between outcome and covariates,
#' specified as Gaussian copula parameter.
#' @param chains Number of Markov chains in MCMC sampling. The default value is
#' \code{chains=4}.
#' @param iter Number of iterations for each chain (including warmup) in MCMC
#' sampling. The default value is \code{iter=2000}.
#' @param warmup Number of warmup (aka burnin) iterations per chain in MCMC
#' sampling. The default value is \code{warmup=floor(iter/2)}.
#' @param thin Period for saving samples in MCMC sampling. The default value
#' is \code{thin=1}.
#' @param a0 Parameter of power prior.
#' @param alternative Alternative hypothesis to be tested ("greater" or "less").
#' The default value is \code{alternative="greater"}.
#' @param sig.level Significance level. The default value is
#' \code{sig.level=0.05}.
#' @param nsim Number of simulated trials. The default value is \code{nsim=100}.
#' @return
#' \item{w}{Likelihood ratio statistics.}
#' \item{p1}{Probability meeting that the treatment effect is above 0
#' (if \code{alternative="greater"}) or below 0 (if \code{alternative="less"})
#' at interim analysis occasion.}
#' \item{p2}{Probability meeting that the treatment effect is above 0
#' (if \code{alternative="greater"}) or below 0 (if \code{alternative="less"})
#' at final analysis occasion.}
#' \item{r1}{Reject the null hypothesis at interim analysis occasion.}
#' \item{r2}{Reject the null hypothesis at final analysis occasion.}
#' @references
#' Psioda MA, Soukup M, Ibrahim JG. A practical Bayesian adaptive design
#' incorporating data from historical controls *Statistics in Medicine*
#' 2018; 37:4054-4070.
#' @import rstan
#' @export

iaborrow.t2e <- function(
  n.CT, n.CC, nevent.C, n.EC, nevent.E, nevent.C1, accrual,
  out.mevent.CT, out.mevent.CC, driftHR,
  cov.CT, cov.CC, cov.EC, cormat,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  a0, alternative="greater", sig.level=0.05, nsim=10)
{
  ncov  <- length(cov.CT)

  out.lambda.CT <- log(2)/out.mevent.CT
  out.lambda.CC <- log(2)/out.mevent.CC
  out.lambda.EC <- out.lambda.CC/driftHR

  ls1 <- length(out.lambda.EC)
  ls2 <- length(out.lambda.CT)

  w  <- array(0,dim=c(nsim,ls1,ls2))
  p1 <- array(0,dim=c(nsim,ls1,ls2))
  p2 <- array(0,dim=c(nsim,ls1,ls2))

  for(ss in 1:nsim){
    for(s1 in 1:ls1){
      for(s2 in 1:ls2){

        marg.CT <- list(list(dist="exp",parm=list(rate=out.lambda.CT[s2])))
        marg.CC <- list(list(dist="exp",parm=list(rate=out.lambda.CC)))
        marg.EC <- list(list(dist="exp",parm=list(rate=out.lambda.EC[s1])))

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
        data.CT <- datagen(margdist=marg.CT,corvec=cvec,nsim=n.CT)
        data.CC <- datagen(margdist=marg.CC,corvec=cvec,nsim=n.CC)
        data.EC <- datagen(margdist=marg.EC,corvec=cvec,nsim=n.EC)

        enroll.period.C <- floor((n.CT+n.CC)/accrual)

        enroll.time.CT <- NULL
        enroll.time.CC <- NULL
        for(i in 1:enroll.period.C){
          tmp <- runif(1)
          if(tmp<0.5){
            tmp.n.CT <- floor(accrual*(n.CT/(n.CT+n.CC)))
          }else{
            tmp.n.CT <- ceiling(accrual*(n.CT/(n.CT+n.CC)))
          }
          tmp.n.CC <- accrual-tmp.n.CT

          e.st <- i-1
          e.en <- i
          enroll.time.CT <- c(enroll.time.CT,runif(tmp.n.CT,e.st,e.en))
          enroll.time.CC <- c(enroll.time.CC,runif(tmp.n.CC,e.st,e.en))
        }
        enroll.time.CT <- c(enroll.time.CT,runif(n.CT-length(enroll.time.CT),enroll.period.C,enroll.period.C+1))
        enroll.time.CC <- c(enroll.time.CC,runif(n.CC-length(enroll.time.CC),enroll.period.C,enroll.period.C+1))

        enroll.period.E <- floor(n.EC/accrual)

        enroll.time.EC <- NULL
        for(i in 1:enroll.period.E){
          e.st <- i-1
          e.en <- i
          enroll.time.EC <- c(enroll.time.EC,runif(accrual,e.st,e.en))
        }
        enroll.time.EC <- c(enroll.time.EC,runif(n.EC-length(enroll.time.EC),enroll.period.E,enroll.period.E+1))

        obs.time.CT <- data.CT[,1]+enroll.time.CT
        obs.time.CC <- data.CC[,1]+enroll.time.CC
        obs.time.EC <- data.EC[,1]+enroll.time.EC

        last.sub.C  <- sort(c(obs.time.CT,obs.time.CC))[nevent.C]
        last.sub.C1 <- sort(c(obs.time.CT,obs.time.CC))[nevent.C1]
        last.sub.E  <- sort(obs.time.EC)[nevent.E]

        censor.CT <- as.numeric(obs.time.CT>last.sub.C)
        censor.CC <- as.numeric(obs.time.CC>last.sub.C)
        censor.EC <- as.numeric(obs.time.EC>last.sub.E)

        data.CT1   <- data.CT[enroll.time.CT<=last.sub.C1,]
        data.CC1   <- data.CC[enroll.time.CC<=last.sub.C1,]
        censor.CT1 <- as.numeric(obs.time.CT>last.sub.C1)[enroll.time.CT<=last.sub.C1]
        censor.CC1 <- as.numeric(obs.time.CC>last.sub.C1)[enroll.time.CC<=last.sub.C1]

        cenf <- (obs.time.CT>last.sub.C)
        data.CT[cenf,1] <- data.CT[cenf,1]-(obs.time.CT[cenf]-last.sub.C)

        cenf <- (obs.time.CC>last.sub.C)
        data.CC[cenf,1] <- data.CC[cenf,1]-(obs.time.CC[cenf]-last.sub.C)

        cenf <- (obs.time.EC>last.sub.E)
        data.EC[cenf,1] <- data.EC[cenf,1]-(obs.time.EC[cenf]-last.sub.E)

        cenf1 <- (enroll.time.CT<=last.sub.C1)
        cenf2 <- (obs.time.CT[cenf1]>last.sub.C1)
        cenf3 <- cenf1&(obs.time.CT>last.sub.C1)
        data.CT1[cenf2,1] <- data.CT1[cenf2,1]-(obs.time.CT[cenf3]-last.sub.C1)

        cenf1 <- (enroll.time.CC<=last.sub.C1)
        cenf2 <- (obs.time.CC[cenf1]>last.sub.C1)
        cenf3 <- cenf1&(obs.time.CC>last.sub.C1)
        data.CC1[cenf2,1] <- data.CC1[cenf2,1]-(obs.time.CC[cenf3]-last.sub.C1)

        dat1 <- list(
          nCT_o = sum(censor.CT1==0),
          nCT_c = sum(censor.CT1==1),
          nCC_o = sum(censor.CC1==0),
          nCC_c = sum(censor.CC1==1),
          p     = ncov,
          yCT_o = data.CT1[censor.CT1==0,1],
          yCT_c = data.CT1[censor.CT1==1,1],
          yCC_o = data.CC1[censor.CC1==0,1],
          yCC_c = data.CC1[censor.CC1==1,1],
          xCT_o = data.CT1[censor.CT1==0,-1],
          xCT_c = data.CT1[censor.CT1==1,-1],
          xCC_o = data.CC1[censor.CC1==0,-1],
          xCC_c = data.CC1[censor.CC1==1,-1])

        dat2 <- list(
          nEC_o = sum(censor.EC==0),
          nEC_c = sum(censor.EC==1),
          p     = ncov,
          yEC_o = data.EC[censor.EC==0,1],
          yEC_c = data.EC[censor.EC==1,1],
          xEC_o = data.EC[censor.EC==0,-1],
          xEC_c = data.EC[censor.EC==1,-1],
          a0    = a0)

        dat3 <- list(
          nCT_o = sum(censor.CT1==0),
          nCT_c = sum(censor.CT1==1),
          nCC_o = sum(censor.CC1==0),
          nCC_c = sum(censor.CC1==1),
          nEC_o = sum(censor.EC==0),
          nEC_c = sum(censor.EC==1),
          p     = ncov,
          yCT_o = data.CT1[censor.CT1==0,1],
          yCT_c = data.CT1[censor.CT1==1,1],
          yCC_o = data.CC1[censor.CC1==0,1],
          yCC_c = data.CC1[censor.CC1==1,1],
          yEC_o = data.EC[censor.EC==0,1],
          yEC_c = data.EC[censor.EC==1,1],
          xCT_o = data.CT1[censor.CT1==0,-1],
          xCT_c = data.CT1[censor.CT1==1,-1],
          xCC_o = data.CC1[censor.CC1==0,-1],
          xCC_c = data.CC1[censor.CC1==1,-1],
          xEC_o = data.EC[censor.EC==0,-1],
          xEC_c = data.EC[censor.EC==1,-1],
          a0    = a0)

        dat4 <- list(
          nCT_o = sum(censor.CT==0),
          nCT_c = sum(censor.CT==1),
          nCC_o = sum(censor.CC==0),
          nCC_c = sum(censor.CC==1),
          p     = ncov,
          yCT_o = data.CT[censor.CT==0,1],
          yCT_c = data.CT[censor.CT==1,1],
          yCC_o = data.CC[censor.CC==0,1],
          yCC_c = data.CC[censor.CC==1,1],
          xCT_o = data.CT[censor.CT==0,-1],
          xCT_c = data.CT[censor.CT==1,-1],
          xCC_o = data.CC[censor.CC==0,-1],
          xCC_c = data.CC[censor.CC==1,-1])

        mcmc1 <- rstan::sampling(stanmodels$T2EConcurrent,
                                 data          = dat1,
                                 chains        = chains,
                                 iter          = iter,
                                 warmup        = warmup,
                                 thin          = thin,
                                 show_messages = FALSE,
                                 cores         = 1,
                                 refresh       = 0)
        mcmc1.sample <- rstan::extract(mcmc1)

        mcmc2 <- rstan::sampling(stanmodels$T2EExternal,
                                 data          = dat2,
                                 chains        = chains,
                                 iter          = iter,
                                 warmup        = warmup,
                                 thin          = thin,
                                 show_messages = FALSE,
                                 cores         = 1,
                                 refresh       = 0)
        mcmc2.sample <- rstan::extract(mcmc2)

        mcmc3 <- rstan::sampling(stanmodels$T2EHybrid,
                                 data          = dat3,
                                 chains        = chains,
                                 iter          = iter,
                                 warmup        = warmup,
                                 thin          = thin,
                                 show_messages = FALSE,
                                 cores         = 1,
                                 refresh       = 0)
        mcmc3.sample <- rstan::extract(mcmc3)

        mcmc4 <- rstan::sampling(stanmodels$T2EConcurrent,
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

        m1.CT_o <- exp(hat.gammaCC1+hat.theta1+(dat1$xCT_o)%*%hat.beta1)
        m1.CC_o <- exp(hat.gammaCC1           +(dat1$xCC_o)%*%hat.beta1)
        m1_o    <- rbind(m1.CT_o,m1.CC_o)

        m2_o    <- exp(hat.gammaCC0           +(dat2$xEC_o)%*%hat.beta0)

        m3.CT_o <- exp(hat.gammaCC +hat.theta +(dat1$xCT_o)%*%hat.beta)
        m3.CC_o <- exp(hat.gammaCC            +(dat1$xCC_o)%*%hat.beta)
        m3_o    <- rbind(m3.CT_o,m3.CC_o)

        m4_o    <- exp(hat.gammaCC            +(dat2$xEC_o)%*%hat.beta)

        m1.CT_c <- exp(hat.gammaCC1+hat.theta1+(dat1$xCT_c)%*%hat.beta1)
        m1.CC_c <- exp(hat.gammaCC1           +(dat1$xCC_c)%*%hat.beta1)
        m1_c    <- rbind(m1.CT_c,m1.CC_c)

        m2_c    <- exp(hat.gammaCC0           +(dat2$xEC_c)%*%hat.beta0)

        m3.CT_c <- exp(hat.gammaCC +hat.theta +(dat1$xCT_c)%*%hat.beta)
        m3.CC_c <- exp(hat.gammaCC            +(dat1$xCC_c)%*%hat.beta)
        m3_c    <- rbind(m3.CT_c,m3.CC_c)

        m4_c    <- exp(hat.gammaCC            +(dat2$xEC_c)%*%hat.beta)

        alpha1 <- median(mcmc1.sample$alpha)
        alpha2 <- median(mcmc2.sample$alpha)
        alpha3 <- median(mcmc3.sample$alpha)
        alpha4 <- median(mcmc3.sample$alpha)

        w[ss,s1,s2] <- (    sum(c(log(dweibull(c(dat1$yCT_o,dat1$yCC_o),alpha1,m1_o)),log(pweibull(c(dat1$yCT_c,dat1$yCC_c),alpha1,m1_c,lower.tail=FALSE))))
                        +a0*sum(c(log(dweibull(c(dat2$yEC_o),           alpha2,m2_o)),log(pweibull(c(dat2$yEC_c),           alpha2,m2_c,lower.tail=FALSE))))
                        -   sum(c(log(dweibull(c(dat1$yCT_o,dat1$yCC_o),alpha3,m3_o)),log(pweibull(c(dat1$yCT_c,dat1$yCC_c),alpha3,m3_c,lower.tail=FALSE))))
                        -a0*sum(c(log(dweibull(c(dat2$yEC_o),           alpha4,m4_o)),log(pweibull(c(dat2$yEC_c),           alpha4,m4_c,lower.tail=FALSE)))))

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
