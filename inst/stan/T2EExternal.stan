data {
  int<lower=1> nEC_o;
  int<lower=1> nEC_c;
  int<lower=1> p;
  vector[nEC_o] yEC_o;
  vector[nEC_c] yEC_c;
  row_vector[p] xEC_o[nEC_o];
  row_vector[p] xEC_c[nEC_c];
  real a0;
}
parameters {
  real gammaCC;
  vector[p] beta;
  real<lower=0> alpha;
}
model {
  alpha ~ exponential(1);

  for (i in 1:nEC_o)
    target += a0*weibull_lpdf(yEC_o[i]|alpha,exp(gammaCC+xEC_o[i]*beta));

  for (i in 1:nEC_c)
    target += a0*weibull_lccdf(yEC_c[i]|alpha,exp(gammaCC+xEC_c[i]*beta));
}
