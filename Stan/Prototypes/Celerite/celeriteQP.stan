functions {
   //real logLikSHO(vector t, vector y,real S0, real w0, real Q, real eps, vector diag);
   vector dotCholQuasiPeriod(vector t, vector y,real B, real L, real P, real C, vector diag);
}

data {
  int<lower=1> N;
  vector[N] t; // time
  vector[N] y; // light curve
  //log uniform priors
  vector[2] B_prior;
  vector[2] L_prior;
  vector[2] P_prior;
  vector[2] C_prior;
  vector[2] err_prior; // precision for the white noise
  // hyperpara
  vector[N] diag;
}


parameters{
   vector[N] eta;
   real<lower = B_prior[1], upper = B_prior[2]> lB;
   real<lower = L_prior[1], upper = L_prior[2]> lL;
   real<lower = P_prior[1], upper = P_prior[2]> lP;
   real<lower = C_prior[1], upper = C_prior[2]> lC;
   real<lower = 0> err; //error variance
}

transformed parameters{
   real B;
   real L;
   real P;
   real C;
   B = exp(lB);
   L = exp(lL);
   P = exp(lP);
   C = exp(lC);
}

model{
   vector[N] trend;
   trend = dotCholQuasiPeriod(t, eta, B, L, P, C,diag);
   err ~ inv_gamma(err_prior[1], err_prior[2]);
   y ~ normal(trend, sqrt(err));
   eta ~ normal(0,1);
   //target += logLikRotation(t, y, sigma, period, 
   //                     Q0, dQ, f, eps, diag + err);

}

generated quantities {
   vector[N] trend;

   trend = dotCholQuasiPeriod(t, eta, B,L,P,C,diag);
}
