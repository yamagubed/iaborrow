
#' Stan function for external data with binary outcome
#'
#' Stan function for logistic regression model using external data is compiled.
#' @return Compiled stan function for external data with binary outcome
#' @import rstan
#' @export

stan.bin.external <- function()
{
  stanmodel <-
  '
  data {
    int<lower=1> nEC;
    int<lower=1> p;
    int<lower=0,upper=1> yEC[nEC];
    row_vector[p] xEC[nEC];
    real a0;
  }
  parameters {
    real gammaCC;
    vector[p] beta;
  }
  model {
    for (i in 1:nEC)
      target += a0*bernoulli_lpmf(yEC[i]|inv_logit(gammaCC+xEC[i]*beta));
  }
  '
  return(rstan::stan_model(model_code=stanmodel))
}
