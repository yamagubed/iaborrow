data {
  int<lower=1> nCT_o;
  int<lower=1> nCT_c;
  int<lower=1> nCC_o;
  int<lower=1> nCC_c;
  int<lower=1> nEC_o;
  int<lower=1> nEC_c;
  int<lower=1> p;
  vector[nCT_o] yCT_o;
  vector[nCT_c] yCT_c;
  vector[nCC_o] yCC_o;
  vector[nCC_c] yCC_c;
  vector[nEC_o] yEC_o;
  vector[nEC_c] yEC_c;
  row_vector[p] xCT_o[nCT_o];
  row_vector[p] xCT_c[nCT_c];
  row_vector[p] xCC_o[nCC_o];
  row_vector[p] xCC_c[nCC_c];
  row_vector[p] xEC_o[nEC_o];
  row_vector[p] xEC_c[nEC_c];
  real a0;
}
parameters {
  real gammaCC;
  vector[p] beta;
  real theta;
  real alpha;
}
model {
  alpha ~ exponential(1);

  for (i in 1:nCT_o)
    target +=    weibull_lpdf(yCT_o[i]|alpha,exp(gammaCC+theta+xCT_o[i]*beta));
  for (i in 1:nCC_o)
    target +=    weibull_lpdf(yCC_o[i]|alpha,exp(gammaCC      +xCC_o[i]*beta));
  for (i in 1:nEC_o)
    target += a0*weibull_lpdf(yEC_o[i]|alpha,exp(gammaCC      +xEC_o[i]*beta));

  for (i in 1:nCT_c)
    target +=    weibull_lccdf(yCT_c[i]|alpha,exp(gammaCC+theta+xCT_c[i]*beta));
  for (i in 1:nCC_c)
    target +=    weibull_lccdf(yCC_c[i]|alpha,exp(gammaCC      +xCC_c[i]*beta));
  for (i in 1:nEC_c)
    target += a0*weibull_lccdf(yEC_c[i]|alpha,exp(gammaCC      +xEC_c[i]*beta));
}
