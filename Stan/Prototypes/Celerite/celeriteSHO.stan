functions {
   //real logLikSHO(vector t, vector y,real S0, real w0, real Q, real eps, vector diag);
   vector dotCholSHO(vector t, vector y,real S0, real w0, real Q, real eps, vector diag);
}

data {
  int<lower=1> N;
  vector[N] t; // time
  vector[N] y; // light curve
  //log uniform priors
  vector[2] S0_prior;
  vector[2] w0_prior;
  vector[2] Q_prior;
  vector[2] err_prior; // precision for the white noise
  // hyperpara
  vector[N] diag;
}

transformed data {
  real eps = 1e-9;
}

parameters{
   vector[N] eta;
   real<lower = S0_prior[1], upper = S0_prior[2]> lS0;
   real<lower = w0_prior[1], upper = w0_prior[2]> lw0;
   real<lower = Q_prior[1], upper = Q_prior[2]> lQ;
   real<lower = 0> err; //error variance
}

transformed parameters{
   real S0;
   real w0;
   real Q;
   S0 = exp(lS0);
   w0 = exp(lw0);
   Q = exp(lQ);
}

model{
   vector[N] trend;
   trend = dotCholSHO(t, eta, S0,w0,Q,eps,diag);
   err ~ inv_gamma(err_prior[1], err_prior[2]);
   y ~ normal(trend, sqrt(err));
   eta ~ normal(0,1);
   //target += logLikRotation(t, y, sigma, period, 
   //                     Q0, dQ, f, eps, diag + err);

}

generated quantities {
   vector[N] trend;

   trend = dotCholSHO(t, eta, S0,w0,Q,eps,diag);
}
