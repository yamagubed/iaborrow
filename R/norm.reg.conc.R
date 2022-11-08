
#' Stan function for normal regression model using concurrent data
#'
#' Stan function for normal regression model using concurrent data is compiled.
#' @return Compiled stan function for normal regression model using concurrent
#' data
#' @import rstan
#' @export

norm.reg.conc <- function()
{
  stanmodel <-
  '
  data {
    int<lower=0> n;
    int<lower=0> p;
    matrix[n,p] x;
    vector[n] z;
    vector[n] y;
  }
  transformed data {
    matrix[n,p] Q_ast;
    matrix[p,p] R_ast;
    matrix[p,p] R_ast_inverse;

    Q_ast         = qr_Q(x)[,1:p]*sqrt(n-1);
    R_ast         = qr_R(x)[1:p,]/sqrt(n-1);
    R_ast_inverse = inverse(R_ast);
  }
  parameters {
    real alpha;
    real gamma;
    vector[p] theta;
    real<lower=0> sigma;
  }
  model {
    target += normal_lpdf(y|Q_ast*theta+alpha+z*gamma,sigma);
  }
  generated quantities {
    vector[p] beta;
    beta = R_ast_inverse*theta;
  }
  '
  return(rstan::stan_model(model_code=stanmodel))
}


