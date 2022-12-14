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
