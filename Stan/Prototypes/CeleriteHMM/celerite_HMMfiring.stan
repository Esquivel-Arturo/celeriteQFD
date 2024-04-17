functions {
   real logLikRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);
}

data {
    int<lower=1> N;
    int<lower=1> K;// classes of states
    vector[N] t; // time
    vector[N] y; // light curve
    //log uniform priors for celerite
    vector[2] sigma_prior;
    vector[2] period_prior;
    vector[2] Q0_prior;
    vector[2] dQ_prior;
    vector[2] f_prior;
    // prior for random firing
    vector<lower = 0>[K] alpha;
    // prior for usual noise
    real mu0;
    real lambda;
    vector[2] noise_prior; // shape_rate for noise
  // hyperpara
    vector[N] diag;
}

transformed data {
  real eps = 1e-9;
}

parameters{
    // trend model
   vector[N] trend;
   real<lower = sigma_prior[1], upper = sigma_prior[2]> lsigma;
   real<lower = period_prior[1], upper = period_prior[2]> lperiod;
   real<lower = Q0_prior[1], upper = Q0_prior[2]> lQ0;
   real<lower = dQ_prior[1], upper = dQ_prior[2]> ldQ;
   real<lower = f_prior[1], upper = f_prior[2]> f;
    // firing model
   simplex[K] theta[K]; // transition matrix
   ordered[K] mus;
   vector<lower = 0>[K] sigma2_err; // noise

}

transformed parameters{
   real sigma;
   real period;
   real Q0;
   real dQ;
   sigma = exp(lsigma);
   period = exp(lperiod);
   Q0 = exp(lQ0);
   dQ = exp(ldQ);
}

model{
   vector[N] yd;// detrended curve
   // prior settings 
      // No need of trend GP model parameters since coded in parameter section 

   // Firing HMM model parameters
   for(i in 1:K){
      sigma2_err[i] ~ inv_gamma(noise_prior[1], noise_prior[2]);
      theta[i] ~ dirichlet(alpha);
      mu[i] ~ normal(mu0, sigma[i]/lambda);
   }

   // likelihood 
   // GP trend
   target += logLikRotation(t, trend, sigma, period, 
                          Q0, dQ, f, eps, diag);

   // vanilla HMM firing model, need to be changed later
   yd = y - trend; // detrended light curve

   real accu[K];
   real gamma[T,K];// joint likelihood of first t states

   for(i in 1:K){
      gamma[1,i] = normal_lpdf(yd[1]|mu[i] + trend[1], sigma[i]);
   }

   for(t in 2:T){
     for(i in 1:K){
         // time t state is i
         for(j in 1:K){
            // and came from j
            accu[j] = gamma[t-1, j] + log(theta[j,i]) + 
                    normal_lpdf(yd[t] | mu[i] + trend[t], sigma[i]);
         }
         gamma[t,i] = log_sum_exp(accu);
      }
   }
   target += log_sum_exp(gamma[T]);
}
