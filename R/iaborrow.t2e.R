
#' Conducting simulation study of Bayesian adaptive design incorporating data
#' from historical controls: time-to-event outcome
#'
#' Bayesian adaptive design proposed by Psioda (2018) is implemented, where
#' historical control data is incorporated at interim decision.
#' The time-to-event outcome is applicable.
#' @usage
#' iaborrow.t2e(
#'   n.CT, n.CC, nevent.C, n.EC, nevent.E, nevent.C1, accrual,
#'   out.mevent.CT, out.mevent.CC, driftHR,
#'   cov.C, cov.cor.C, cov.effect.C,
#'   cov.E, cov.cor.E, cov.effect.E,
#'   chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   a0, alternative="greater", sig.level=0.025, nsim)
#' @param n.CT Number of patients in concurrent treatment at final analysis
#' opportunity.
#' @param n.CC Number of patients in concurrent control at final analysis
#' opportunity.
#' @param nevent.C Number of events in concurrent data at final analysis
#' opportunity.
#' @param n.EC Number of patients in external control.
#' @param nevent.E Number of events in external control.
#' @param nevent.C1 Number of events in concurrent data at interim analysis
#' opportunity.
#' @param accrual Accrual rate (number of enrolled patients per month).
#' @param out.mevent.CT True median time to event in concurrent treatment.
#' @param out.mevent.CC True median time to event in concurrent control.
#' @param driftHR Hazard ratio between external control and concurrent control
#' for which the bias should be plotted.
#' @param cov.C List of covariates in concurrent data. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.cor.C Matrix of correlation coefficients for covariates in
#' concurrent data, specified as Gaussian copula parameter.
#' @param cov.effect.C Vector of covariate effects in concurrent data,
#' specified as hazard ratio.
#' @param cov.E List of covariates in external data. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.cor.E Matrix of correlation coefficients for covariates in
#' external data, specified as Gaussian copula parameter.
#' @param cov.effect.E Vector of covariate effects in external data,
#' specified as hazard ratio.
#' @param chains Number of Markov chains in MCMC sampling. The default value is
#' \code{chains=2}.
#' @param iter Number of iterations for each chain (including warmup) in MCMC
#' sampling. The default value is \code{iter=4000}.
#' @param warmup Number of warmup (aka burnin) iterations per chain in MCMC
#' sampling. The default value is \code{warmup=floor(iter/2)}.
#' @param thin Period for saving samples in MCMC sampling. The default value
#' is \code{thin=1}.
#' @param a0 Parameter of power prior.
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
#' \item{n.CT}{Number of patients in concurrent treatment at final analysis
#' opportunity.}
#' \item{n.CT1}{Number of patients in concurrent treatment at interim analysis
#' opportunity.}
#' \item{n.CC}{Number of patients in concurrent control at final analysis
#' opportunity.}
#' \item{n.CC1}{Number of patients in concurrent control at interim analysis
#' opportunity.}
#' \item{n.EC}{Number of patients in external control.}
#' @references
#' Psioda MA, Soukup M, Ibrahim JG. A practical Bayesian adaptive design
#' incorporating data from historical controls *Statistics in Medicine*
#' 2018; 37:4054-4070.
#' @examples
#' n.CT      <- 60
#' n.CC      <- 60
#' nevent.C  <- 100
#' n.EC      <- 120
#' nevent.E  <- 100
#' nevent.C1 <- 100
#' accrual   <- 16
#'
#' out.mevent.CT <- c(6,9)
#' out.mevent.CC <- 6
#' driftHR       <- c(0.8,1.0,1.2)
#'
#' cov.C <- list(list(dist="norm", mean=0,sd=1),
#'               list(dist="binom",prob=0.4))
#'
#' cov.cor.C <- rbind(c(  1,0.1),
#'                    c(0.1,  1))
#'
#' cov.effect.C <- c(0.1,0.1)
#'
#' cov.E <- list(list(dist="norm", mean=0,sd=1),
#'               list(dist="binom",prob=0.4))
#'
#' cov.cor.E <- rbind(c(  1,0.1),
#'                    c(0.1,  1))
#'
#' cov.effect.E <- c(0.1,0.1)
#'
#' a0 <- 0.5
#'
#' nsim <- 10
#'
#' iaborrow.t2e(
#'   n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
#'   n.EC=n.EC, nevent.E=nevent.E, nevent.C1=nevent.C1, accrual=accrual,
#'   out.mevent.CT=out.mevent.CT, out.mevent.CC=out.mevent.CC, driftHR=driftHR,
#'   cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
#'   cov.E=cov.E, cov.cor.E=cov.cor.E, cov.effect.E=cov.effect.E,
#'   a0=a0, alternative="less", nsim=nsim)
#' @import rstan
#' @export

iaborrow.t2e <- function(
  n.CT, n.CC, nevent.C, n.EC, nevent.E, nevent.C1, accrual,
  out.mevent.CT, out.mevent.CC, driftHR,
  cov.C, cov.cor.C, cov.effect.C,
  cov.E, cov.cor.E, cov.effect.E,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  a0, alternative="greater", sig.level=0.025, nsim)
{
  ncov          <- length(cov.C)
  out.lambda.CT <- log(2)/out.mevent.CT
  out.lambda.CC <- log(2)/out.mevent.CC
  out.lambda.EC <- out.lambda.CC*driftHR

  ls1 <- length(out.lambda.EC)
  ls2 <- length(out.lambda.CT)

  lcov.effect.C <- (-log(cov.effect.C))
  lcov.effect.E <- (-log(cov.effect.E))

  marg.C <- NULL
  marg.E <- NULL
  mean.C <- NULL
  mean.E <- NULL

  for(i in 1:ncov){
    if(cov.C[[i]]$dist=="norm"){
      marg.C <- append(marg.C,list(list(dist=cov.C[[i]]$dist,parm=list(mean=cov.C[[i]]$mean,sd=cov.C[[i]]$sd))))
      marg.E <- append(marg.E,list(list(dist=cov.E[[i]]$dist,parm=list(mean=cov.E[[i]]$mean,sd=cov.E[[i]]$sd))))

      mean.C <- c(mean.C,cov.C[[i]]$mean)
      mean.E <- c(mean.E,cov.E[[i]]$mean)
    }else if(cov.C[[i]]$dist=="binom"){
      marg.C <- append(marg.C,list(list(dist=cov.C[[i]]$dist,parm=list(size=1,prob=cov.C[[i]]$prob))))
      marg.E <- append(marg.E,list(list(dist=cov.E[[i]]$dist,parm=list(size=1,prob=cov.E[[i]]$prob))))

      mean.C <- c(mean.C,cov.C[[i]]$prob)
      mean.E <- c(mean.E,cov.E[[i]]$prob)
    }
  }

  int.C   <- log(1/out.lambda.CC)-sum(mean.C*lcov.effect.C)
  int.E   <- log(1/out.lambda.EC)-sum(mean.E*lcov.effect.E)
  t.theta <- log(out.lambda.CC/out.lambda.CT)

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.E  <- cov.cor.E[lower.tri(cov.cor.E)]

  w     <- array(0,dim=c(nsim,ls1,ls2))
  p1    <- array(0,dim=c(nsim,ls1,ls2))
  p2    <- array(0,dim=c(nsim,ls1,ls2))
  n.CT1 <- array(0,dim=c(nsim,ls1,ls2))
  n.CC1 <- array(0,dim=c(nsim,ls1,ls2))

  for(ss in 1:nsim){
    for(s1 in 1:ls1){
      for(s2 in 1:ls2){

        data.cov.CT <- datagen(margdist=marg.C,corvec=cvec.C,nsim=n.CT)
        data.cov.CC <- datagen(margdist=marg.C,corvec=cvec.C,nsim=n.CC)
        data.cov.EC <- datagen(margdist=marg.E,corvec=cvec.E,nsim=n.EC)

        sigma.CT <- exp(int.C+t.theta+apply(data.cov.CT,1,function(x){sum(x*lcov.effect.C)}))
        sigma.CC <- exp(int.C        +apply(data.cov.CC,1,function(x){sum(x*lcov.effect.C)}))
        sigma.EC <- exp(int.E        +apply(data.cov.EC,1,function(x){sum(x*lcov.effect.E)}))

        data.CT <- cbind(rweibull(n.CT,shape=1,scale=sigma.CT),data.cov.CT)
        data.CC <- cbind(rweibull(n.CC,shape=1,scale=sigma.CC),data.cov.CC)
        data.EC <- cbind(rweibull(n.EC,shape=1,scale=sigma.EC),data.cov.EC)

        enroll.period.C <- floor((n.CT+n.CC)/accrual)
        mn.CT <- floor(n.CT/enroll.period.C)
        mn.CC <- floor(n.CC/enroll.period.C)

        enroll.time.CT <- NULL
        enroll.time.CC <- NULL
        for(i in 1:enroll.period.C){
          e.st <- i-1
          e.en <- i
          enroll.time.CT <- c(enroll.time.CT,runif(mn.CT,e.st,e.en))
          enroll.time.CC <- c(enroll.time.CC,runif(mn.CC,e.st,e.en))
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

        censor.CT <- censor.CT[data.CT[,1]>0]
        data.CT   <- data.CT[data.CT[,1]>0,]

        censor.CC <- censor.CC[data.CC[,1]>0]
        data.CC <- data.CC[data.CC[,1]>0,]

        censor.EC <- censor.EC[data.EC[,1]>0]
        data.EC <- data.EC[data.EC[,1]>0,]

        dat1 <- list(
          nCT_o = sum(censor.CT1==0),
          nCT_c = sum(censor.CT1==1),
          nCC_o = sum(censor.CC1==0),
          nCC_c = sum(censor.CC1==1),
          p     = ncov,
          yCT_o = as.array(data.CT1[censor.CT1==0,1]),
          yCT_c = as.array(data.CT1[censor.CT1==1,1]),
          yCC_o = as.array(data.CC1[censor.CC1==0,1]),
          yCC_c = as.array(data.CC1[censor.CC1==1,1]),
          xCT_o = matrix(data.CT1[censor.CT1==0,-1],sum(censor.CT1==0),ncov),
          xCT_c = matrix(data.CT1[censor.CT1==1,-1],sum(censor.CT1==1),ncov),
          xCC_o = matrix(data.CC1[censor.CC1==0,-1],sum(censor.CC1==0),ncov),
          xCC_c = matrix(data.CC1[censor.CC1==1,-1],sum(censor.CC1==1),ncov))

        dat2 <- list(
          nEC_o = sum(censor.EC==0),
          nEC_c = sum(censor.EC==1),
          p     = ncov,
          yEC_o = as.array(data.EC[censor.EC==0,1]),
          yEC_c = as.array(data.EC[censor.EC==1,1]),
          xEC_o = matrix(data.EC[censor.EC==0,-1],sum(censor.EC==0),ncov),
          xEC_c = matrix(data.EC[censor.EC==1,-1],sum(censor.EC==1),ncov),
          a0    = a0)

        dat3 <- list(
          nCT_o = sum(censor.CT1==0),
          nCT_c = sum(censor.CT1==1),
          nCC_o = sum(censor.CC1==0),
          nCC_c = sum(censor.CC1==1),
          nEC_o = sum(censor.EC==0),
          nEC_c = sum(censor.EC==1),
          p     = ncov,
          yCT_o = as.array(data.CT1[censor.CT1==0,1]),
          yCT_c = as.array(data.CT1[censor.CT1==1,1]),
          yCC_o = as.array(data.CC1[censor.CC1==0,1]),
          yCC_c = as.array(data.CC1[censor.CC1==1,1]),
          yEC_o = as.array(data.EC[censor.EC==0,1]),
          yEC_c = as.array(data.EC[censor.EC==1,1]),
          xCT_o = matrix(data.CT1[censor.CT1==0,-1],sum(censor.CT1==0),ncov),
          xCT_c = matrix(data.CT1[censor.CT1==1,-1],sum(censor.CT1==1),ncov),
          xCC_o = matrix(data.CC1[censor.CC1==0,-1],sum(censor.CC1==0),ncov),
          xCC_c = matrix(data.CC1[censor.CC1==1,-1],sum(censor.CC1==1),ncov),
          xEC_o = matrix(data.EC[censor.EC==0,-1],sum(censor.EC==0),ncov),
          xEC_c = matrix(data.EC[censor.EC==1,-1],sum(censor.EC==1),ncov),
          a0    = a0)

        if(sum(censor.CC==1)>0){

          dat4 <- list(
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            nCC_c = sum(censor.CC==1),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            yCC_c = as.array(data.CC[censor.CC==1,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov),
            xCC_c = matrix(data.CC[censor.CC==1,-1],sum(censor.CC==1),ncov))

        }else if(sum(censor.CC==1)==0){

          dat4 <- list(
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov))
        }

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

        if(sum(censor.CC==1)>0){

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

        }else if(sum(censor.CC==1)==0){

          mcmc4 <- rstan::sampling(stanmodels$T2EConcurrentC0,
                                   data          = dat4,
                                   chains        = chains,
                                   iter          = iter,
                                   warmup        = warmup,
                                   thin          = thin,
                                   show_messages = FALSE,
                                   cores         = 1,
                                   refresh       = 0)
          mcmc4.sample <- rstan::extract(mcmc4)

        }

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

        loghr3 <- (-mcmc3.sample$alpha*mcmc3.sample$theta)
        loghr4 <- (-mcmc4.sample$alpha*mcmc4.sample$theta)

        if(alternative=="greater"){
          p1[ss,s1,s2] <- mean(loghr3>0)
          p2[ss,s1,s2] <- mean(loghr4>0)
        }else if(alternative=="less"){
          p1[ss,s1,s2] <- mean(loghr3<0)
          p2[ss,s1,s2] <- mean(loghr4<0)
        }

        n.CT[ss,s1,s2]  <- nrow(data.CT)
        n.CT1[ss,s1,s2] <- nrow(data.CT1)

        n.CC[ss,s1,s2]  <- nrow(data.CC)
        n.CC1[ss,s1,s2] <- nrow(data.CC1)

        n.EC[ss,s1,s2]  <- nrow(data.EC)
  }}}

  r1 <- (p1>(1-sig.level))
  r2 <- (p2>(1-sig.level))

  return(list(w=w,p1=p1,p2=p2,r1=r1,r2=r2,
              n.CT=n.CT,n.CT1=n.CT1,
              n.CC=n.CC,n.CC1=n.CC1,n.EC=n.EC))
}
