
#' Stan function for normal regression model using external data
#'
#' Stan function for normal regression model using external data is compiled.
#' @return Compiled stan function for normal regression model using external
#' data
#' @import rstan
#' @export

norm.reg.ext <- function()
{
  stanmodel <-
  '
  data {
    int<lower=0> nHC;
    int<lower=0> pHC;
    matrix[nHC,pHC] xHC;
    vector[nHC] yHC;
    real a0;
  }
  transformed data {
    matrix[nHC,pHC] Q_astHC;
    matrix[pHC,pHC] R_astHC;
    matrix[pHC,pHC] R_astHC_inverse;

    Q_astHC         = qr_Q(xHC)[,1:pHC]*sqrt(nHC-1);
    R_astHC         = qr_R(xHC)[1:pHC,]/sqrt(nHC-1);
    R_astHC_inverse = inverse(R_astHC);
  }
  parameters {
    real alpha;
    vector[pHC] theta;
    real<lower=0> sigma;
  }
  model {
    target += a0*normal_lpdf(yHC|Q_astHC*theta+alpha,sigma);
  }
  generated quantities {
    vector[pHC] beta;
    beta = R_astHC_inverse*theta;
  }
  '
  return(rstan::stan_model(model_code=stanmodel))
}



