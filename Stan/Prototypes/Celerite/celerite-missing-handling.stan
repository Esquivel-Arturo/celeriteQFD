functions {
   //real logLikSHO(vector t, vector y,real S0, real w0, real Q, real eps, vector diag);
   vector dotCholRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);
   real logLikRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);


}

data {
  int<lower=1> N;
  vector[N] t; // time
  vector[N] y; // light curve
  int<lower = 0, upper = 1> observed[N];
  //log uniform priors
  vector[2] sigma_prior;
  vector[2] period_prior;
  vector[2] Q0_prior;
  vector[2] dQ_prior;
  vector[2] f_prior;
  vector[2] err_prior; // precision for the white noise
  // hyperpara
  vector[N] diag;
}

transformed data {
  real eps = 1e-9;
}

parameters{
   vector[N] eta;
   real<lower = sigma_prior[1], upper = sigma_prior[2]> lsigma;
   real<lower = period_prior[1], upper = period_prior[2]> lperiod;
   real<lower = Q0_prior[1], upper = Q0_prior[2]> lQ0;
   real<lower = dQ_prior[1], upper = dQ_prior[2]> ldQ;
   real<lower = f_prior[1], upper = f_prior[2]> f;
   real<lower = 0> err; //error variance
}

transformed parameters{
   real sigma;
   real period;
   real Q0;
   real dQ;
   sigma = exp(lsigma);
   period = exp(lperiod);
   Q0 = exp(lQ0);
   //Q0=-0.25;
   dQ = exp(ldQ);
}

model{
   vector[N] trend;
   trend = dotCholRotation(t, eta, sigma, period, 
                          Q0, dQ, f, eps, diag);
   err ~ inv_gamma(err_prior[1], err_prior[2]);
   //y ~ normal(trend, sqrt(err));
   eta ~ normal(0,1);
   //target += logLikRotation(t, y, sigma, period, 
   //                     Q0, dQ, f, eps, diag + err);
   for(tt in 1:N){
      target += observed[tt] * normal_lpdf(y[tt]|trend[tt], sqrt(err));
   }

}

generated quantities {
   real Q1;
   real w1;
   real S1;
   real Q2;
   real w2;
   real S2;

   real sigma1;
   real rho1;
   real tau1;
   
   real sigma2;
   real rho2;
   real tau2;
   
   real amp;
   vector[N] trend;

   amp =  sigma * sigma/ (1 + f);
   Q1 = 0.5 + Q0 + dQ;
   w1 = 4 * 3.1415926 * Q1 / (period * sqrt(4 * Q1 * Q1 - 1));
   S1 = amp / (w1 * Q1);
   Q2 = 0.5 + Q0;
   w2 = 8 * 3.1415926 * Q2 / (period * sqrt(4 * Q2 * Q2 - 1));
   S2 = f * amp / (w2 * Q2);
   
   rho1 = 2*3.1415926/w1;
   rho2 = 2*3.1415926/w2;
   tau1 = 2 * Q1/w1;
   tau2 = 2 * Q2/w2;
   sigma1 = sqrt(S1 * w1 * Q1);
   sigma2 = sqrt(S2 * w2 * Q2);


   
   trend = dotCholRotation(t, eta, sigma, period, 
                          Q0, dQ, f, eps, diag);
}
