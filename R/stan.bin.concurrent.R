
#' Stan function for concurrent data with binary outcome
#'
#' Stan function for logistic regression model using concurrent data is compiled.
#' @return Compiled stan function for concurrent data with binary outcome
#' @import rstan
#' @export

stan.bin.concurrent <- function()
{
  stanmodel <-
  '
  data {
    int<lower=1> nCT;
    int<lower=1> nCC;
    int<lower=1> p;
    int<lower=0,upper=1> yCT[nCT];
    int<lower=0,upper=1> yCC[nCC];
    row_vector[p] xCT[nCT];
    row_vector[p] xCC[nCC];
  }
  parameters {
    real gammaCC;
    vector[p] beta;
    real theta;
  }
  model {
    for (i in 1:nCT)
      target += bernoulli_lpmf(yCT[i]|inv_logit(gammaCC+theta+xCT[i]*beta));
    for (i in 1:nCC)
      target += bernoulli_lpmf(yCC[i]|inv_logit(gammaCC      +xCC[i]*beta));
  }
  '
  return(rstan::stan_model(model_code=stanmodel))
}
