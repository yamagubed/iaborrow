data {
  int<lower=1> nCT;
  int<lower=1> nCC;
  int<lower=1> p;
  vector[nCT] yCT;
  vector[nCC] yCC;
  row_vector[p] xCT[nCT];
  row_vector[p] xCC[nCC];
}
parameters {
  real gammaCC;
  vector[p] beta;
  real theta;
  real<lower=0> sigma;
}
model {
  for (i in 1:nCT)
    target += normal_lpdf(yCT[i]|gammaCC+theta+xCT[i]*beta,sigma);
  for (i in 1:nCC)
    target += normal_lpdf(yCC[i]|gammaCC      +xCC[i]*beta,sigma);
}
