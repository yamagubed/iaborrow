
#' Stan function for normal regression model using hybrid data
#'
#' Stan function for normal regression model using hybrid data is compiled.
#' @return Compiled stan function for normal regression model using hybrid
#' data
#' @import rstan
#' @export

norm.reg.hyb <- function()
{
  stanmodel <-
  '
  data {
    int<lower=0> n;
    int<lower=0> p;
    matrix[n,p] x;
    vector[n] z;
    vector[n] y;
    int<lower=0> nHC;
    matrix[nHC,p] xHC;
    vector[nHC] yHC;
    real a0;
  }
  transformed data {
    matrix[n,p] Q_ast;
    matrix[p,p] R_ast;
    matrix[p,p] R_ast_inverse;
    matrix[nHC,p] Q_astHC;
    matrix[p,p] R_astHC;
    matrix[p,p] R_astHC_inverse;

    Q_ast         = qr_Q(x)[,1:p]*sqrt(n-1);
    R_ast         = qr_R(x)[1:p,]/sqrt(n-1);
    R_ast_inverse = inverse(R_ast);
    Q_astHC         = qr_Q(xHC)[,1:p]*sqrt(nHC-1);
    R_astHC         = qr_R(xHC)[1:p,]/sqrt(nHC-1);
    R_astHC_inverse = inverse(R_astHC);
  }
  parameters {
    real alpha;
    real gamma;
    vector[p] theta;
    real<lower=0> sigma;
  }
  model {
    target += normal_lpdf(y|Q_ast*theta+alpha+z*gamma,sigma);
    target += a0*normal_lpdf(yHC|Q_astHC*theta+alpha,sigma);
  }
  generated quantities {
    vector[p] beta;
    beta = R_astHC_inverse*theta;
  }
  '
  return(rstan::stan_model(model_code=stanmodel))
}

