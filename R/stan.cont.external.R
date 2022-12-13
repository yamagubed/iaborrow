
#' Stan function for normal regression model using external data
#'
#' Stan function for normal regression model using external data is compiled.
#' @return Compiled stan function for normal regression model using external
#' data
#' @import rstan
#' @export

stan.cont.external <- function()
{
  stanmodel <-
  '
  data {
    int<lower=1> nEC;
    int<lower=1> p;
    vector[nEC] yEC;
    row_vector[p] xEC[nEC];
    real a0;
  }
  parameters {
    real gammaCC;
    vector[p] beta;
    real<lower=0> sigma;
  }
  model {
    for (i in 1:nEC)
      target += a0*normal_lpdf(yEC[i]|gammaCC+xEC[i]*beta,sigma);
  }
  '
  return(rstan::stan_model(model_code=stanmodel))
}



