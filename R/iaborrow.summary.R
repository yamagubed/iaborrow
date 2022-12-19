
#' Summarizing simulation study results
#'
#' Simulation study results are summarized.
#' @usage
#' iaborrow.summary(
#'   n.CT, n.CC, n.CT1, n.CC1, ndrift,
#'   w, p1, p2, r1, r2,
#'   w0, accept.t1e=0.1, accept.pow=0.7)
#' @param n.CT Number of patients in concurrent treatment at final analysis
#' opportunity.
#' @param n.CC Number of patients in concurrent control at final analysis
#' opportunity.
#' @param n.CT1 Number of patients in concurrent treatment at interim analysis
#' opportunity.
#' @param n.CC1 Number of patients in concurrent control at interim analysis
#' opportunity.
#' @param ndrift Number of drift scenarios.
#' @param w Likelihood ratio statistics.
#' @param p1 Probability meeting that the treatment effect is above 0
#' (if \code{alternative="greater"}) or below 0 (if \code{alternative="less"})
#' at interim analysis occasion.
#' @param p2 Probability meeting that the treatment effect is above 0
#' (if \code{alternative="greater"}) or below 0 (if \code{alternative="less"})
#' at final analysis occasion.
#' @param r1 \code{TRUE} when significant at interim analysis opportunity;
#' otherwise \code{FALSE}.
#' @param r2 \code{TRUE} when significant at final analysis opportunity;
#' otherwise \code{FALSE}.
#' @param w0 Candidate values of stoppage criteria.
#' @param accept.t1e Maximum acceptable type I error rate. The default value is
#' \code{accept.t1e=0.1}.
#' @param accept.pow Minimum acceptable power. The default value is
#' \code{accept.pow=0.7}.
#' @export

iaborrow.summary <- function(
  n.CT, n.CC, n.CT1, n.CC1, ndrift,
  w, p1, p2, r1, r2,
  w0, accept.t1e=0.1, accept.pow=0.7)
{
  lw       <- length(w0)
  reject   <- array(0,dim=c(2,lw,ndrift))
  max.t1e  <- numeric(lw)
  min.pow  <- numeric(lw)
  stoppage <- array(0,dim=c(2,lw,ndrift))

  for(i in 1:lw){
    reject[1,i,]   <- apply(r1[,,1]*(w[,,1]<=w0[i])+r2[,,1]*(w[,,1]>w0[i]),2,mean)
    reject[2,i,]   <- apply(r1[,,2]*(w[,,2]<=w0[i])+r2[,,2]*(w[,,2]>w0[i]),2,mean)
    max.t1e[i]     <- max(reject[1,i,])
    min.pow[i]     <- min(reject[2,i,])
    stoppage[1,i,] <- apply(w[,,1]<=w0[i],2,mean)
    stoppage[2,i,] <- apply(w[,,2]<=w0[i],2,mean)
  }

  w0.opt <- max(w0[(max.t1e<=accept.t1e)&(min.pow>=accept.pow)])
  w0.flg <- which(w0==w0.opt)

  t1e           <- reject[1,w0.flg,]
  stop.null     <- apply(w[,,1]<=w0.opt,2,mean)
  exp.size.null <- stop.null*(n.CT1+n.CC1)+(1-stop.null)*(n.CT+n.CC)

  pow           <- reject[2,w0.flg,]
  stop.altr     <- apply(w[,,2]<=w0.opt,2,mean)
  exp.size.altr <- stop.altr*(n.CT1+n.CC1)+(1-stop.altr)*(n.CT+n.CC)

  return(list(reject=reject,max.t1e=max.t1e,min.pow=min.pow,
              stoppage=stoppage,w0.opt=w0.opt,
              t1e=t1e,stop.null=stop.null,exp.size.null=exp.size.null,
              pow=pow,stop.altr=stop.altr,exp.size.altr=exp.size.altr))
}
