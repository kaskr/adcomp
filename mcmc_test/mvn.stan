data {
  int<lower=0> Npar;
  matrix[Npar,Npar] covar;
  vector[Npar] x;
}
parameters {
  vector[Npar] X;
}

model {
  x~multi_normal(X, covar);
}
