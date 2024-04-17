functions {
   real logLikRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);
   vector dotCholRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);

   
}

data {
    int<lower=1> N;
    vector[N] t; // time
    vector[N] y; // light curve
    //log uniform priors for celerite
    vector[2] sigma_prior;
    vector[2] period_prior;
    vector[2] Q0_prior;
    vector[2] dQ_prior;
    vector[2] f_prior;
    // prior for random firing
    //  transitioning
    vector<lower = 0>[2] alpha;
    //   prior for "usual" noise
    real mu0;
    real lambda;
    vector[2] noise_prior; // shape_rate for noise
    // prior on the powerlaw
    real<lower = 0> rate;// shape for the exponential of the power law
    
    vector[N] diag;// hyperpara, usually set to 0
}

transformed data {
  real eps = 1e-9;
  real eps_pare = 2e-3;
}

parameters{
    // trend model
   vector[N] eta;
   real<lower = sigma_prior[1], upper = sigma_prior[2]> lsigma;
   real<lower = period_prior[1], upper = period_prior[2]> lperiod;
   real<lower = Q0_prior[1], upper = Q0_prior[2]> lQ0;
   real<lower = dQ_prior[1], upper = dQ_prior[2]> ldQ;
   real<lower = f_prior[1], upper = f_prior[2]> f;
    // firing model
   simplex[2] theta[2]; // transition matrix
   real mu; // mean of usual state
   real<lower = 0> sigma2_usual; // usual noise
   real<lower = 0> k;// power law exponent

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
    vector[N] trend;
   vector[N] yd;// detrended curve
   real accu[2];
   real gamma[N,2];// joint likelihood of first t states
   // get the trend
   trend = dotCholRotation(t, eta, sigma, period, 
                          Q0, dQ, f, eps, diag);
   eta ~ normal(0,1);
   // prior settings 
      // No need of trend GP model parameters since coded in parameter section 

   // Firing HMM model parameters
    sigma2_usual ~ inv_gamma(noise_prior[1], noise_prior[2]);
    theta[1] ~ dirichlet(alpha);
    theta[2] ~ dirichlet(alpha);
    mu ~ normal(mu0, sigma2_usual/lambda);
    k ~ exponential(rate);

   // likelihood 

   // vanilla HMM firing model, need to be changed later
   yd = y - trend - mu; // detrended light curve, also minus the usual mean
    //yd = y-mu;
   

    // usual state
    gamma[1,1] = normal_lpdf(yd[1]|0, sigma2_usual);
    // firing
    gamma[1,2] = double_exponential_lpdf(yd[1] |0 , k);
    for(tt in 2:N){
        // at usual state
        //  came from usual state
        accu[1] = gamma[tt-1,1] + log(theta[1,1]) + 
                  normal_lpdf(yd[tt]|0, sigma2_usual);
        //  came from firing state
        accu[2] = gamma[tt-1,2] + log(theta[2,1]) + 
                  normal_lpdf(yd[tt]|0, sigma2_usual);
        gamma[tt,1] = log_sum_exp(accu);
        // at firing state
        //  came from usual state
        accu[1] = gamma[tt-1,1] + log(theta[1,2]) + 
                  double_exponential_lpdf(yd[tt]|0 , k);
        accu[2] = gamma[tt-1,2] + log(theta[2,2]) + 
                  double_exponential_lpdf(yd[tt]|0 , k); 
        gamma[tt,2] = log_sum_exp(accu);            
    }
   target += log_sum_exp(gamma[N]);
}

// Viterbi
generated quantities {
    int <lower = 1, upper = 2>state[N];
    vector[N] trend;
    vector[N] yd;// detrended curve
    real log_p_state;
    
    {
        
        int back_ptr[N, 2];
        real best_logp[N, 2];
        real best_total_logp;
        real logp;
        trend = dotCholRotation(t, eta, sigma, period, 
                          Q0, dQ, f, eps, diag);
        yd = y - mu - trend ;
        // usual
        best_logp[1,1] = normal_lpdf(yd[1]|0, sigma2_usual);
        // firing
        best_logp[1,2] = double_exponential_lpdf(yd[1]|0 , k);
        
        for(tt in 2:N){
            best_logp[tt, 1] = negative_infinity();
            logp = best_logp[tt-1, 1] + log(theta[1,1]) + 
                    normal_lpdf(yd[tt]|0, sigma2_usual);
            if(logp>best_logp[tt,1]){
                back_ptr[tt,1] = 1;
                best_logp[tt,1] = logp;
            }
            logp = best_logp[tt-1,2] + log(theta[2,1]) + 
                    normal_lpdf(yd[tt]|0, sigma2_usual);
            if(logp>best_logp[tt,1]){
                back_ptr[tt,1] = 2;
                best_logp[tt,1] = logp;
            }

            best_logp[tt, 2] = negative_infinity();

            logp = best_logp[tt-1, 1] + log(theta[1,2]) + 
                    double_exponential_lpdf(yd[tt]|0 , k);
            if(logp>best_logp[tt,2]){
                back_ptr[tt,2] = 1;
                best_logp[tt,2] = logp;
            }
            logp = best_logp[tt-1,2] + log(theta[2,2]) + 
                    double_exponential_lpdf(yd[tt]|0 , k);
            if(logp>best_logp[tt,2]){
                back_ptr[tt,2] = 2;
                best_logp[tt,2] = logp;
            }

        }

        log_p_state = max(best_logp[N]);
        for(i in 1:2){
            if(best_logp[N, i]==log_p_state){
                state[N] = i;
            }
        }
        for(tt in 1:(N-1)){
            state[N-tt] = back_ptr[N-tt+1, state[N-tt+1]];
        }
    }
}
