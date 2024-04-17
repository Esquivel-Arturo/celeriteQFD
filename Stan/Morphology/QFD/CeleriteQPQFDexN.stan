// Celerite with Linear Firing Exponential Decay
// Quiet-Firing-Decay (QFD model)
functions {
    vector dotCholQuasiPeriod(vector t, vector y,real B, real L, real P, real C, vector diag);
    real sigmoid(real x){
        return(1/(exp(-x)+1));
    }
}

data {
    int<lower=1> N;
    vector[N] t; // time
    vector[N] y; // light curve
    //log uniform priors for celerite
    vector[2] B_prior;
    vector[2] L_prior;
    vector[2] P_prior;
    vector[2] C_prior;
    // prior for random firing
    //  transitioning
    // quiet can only goto iself or firing 
    vector<lower = 0>[2] alpha_quiet;
    // firing can only goto itself or decay
    vector<lower = 0>[2] alpha_firing;
    // decay can go anywhere
    vector<lower = 0>[3] alpha_decay;
    //   prior for quite noise
    real mu0_quiet;
    real lambda_quiet;
    vector[2] gamma_noise; // shape_rate for noise
    // prior on linear increasing, a normal random walk, with positive offset and slope 1
    real mu0_rate_firing; // prior on how much increase
    real<lower=0> sigma_rate_firing; // prior on how much increase
    // prior on decreasing, a AR with no offset and slope from 0 to 1
    real mu0_rate_decay;// I will put a beta prior here
    real<lower = 0> sigma_rate_decay;
    
    vector[N] diag;// hyperpara, usually set to 0
}

transformed data {
    real eps = 1e-9;
}

parameters {
    // trend parameter
    vector[N] eta;
    real<lower = B_prior[1], upper = B_prior[2]> lB;
    real<lower = L_prior[1], upper = L_prior[2]> lL;
    real<lower = P_prior[1], upper = P_prior[2]> lP;
    real<lower = C_prior[1], upper = C_prior[2]> lC;
    // quiet parameter
    simplex[2] theta_quiet; // transitioning probability, 1. to quiet, 2. to firing
    real mu_quiet; // work as grand mu
    real<lower = 0> sigma2_noise; // quiet state variance
    // firing parameter
    simplex[2] theta_firing;// 1. to firing, 2. to decay
    real lograte_firing; // must be increasing on average
    // decay parameter
    simplex[3] theta_decay;// 1. to rest, 2. to firing, 3. to decay
    real logitrate_decay; // the exponent of decaying, must be between 0 and 1 (to decrease)
}

transformed parameters{
   real B;
   real L;
   real P;
   real C;
   real rate_firing;
   real rate_decay;
   B = exp(lB);
   L = exp(lL);
   P = exp(lP);
   C = exp(lC);
   rate_firing = exp(lograte_firing);
   rate_decay = sigmoid(logitrate_decay);
}

model{
    vector[N] trend;
    vector[N] yd;// detrended curve
    real accu_quiet[2];// used for quiet
    real accu_firing[3];// used for firing, just to avoid -Inf firing can come from anywhere
    real accu_decay[2]; // used for decay
    real gamma[N-1,3];// joint likelihood of first t states
    // get the trend
    trend = dotCholQuasiPeriod(t, eta, B, L, P, C,diag);
    eta ~ normal(0,1);
    // prior settings 
      // No need of trend GP model parameters since coded in parameter section 
      // Firing HMM model parameters
        // AR of states
    sigma2_noise ~ inv_gamma(gamma_noise[1], gamma_noise[2]);
        // mean parameters
    mu_quiet ~ normal(mu0_quiet, sqrt(sigma2_noise/lambda_quiet));// serve as overall mean 
    //increm_firing ~ exponential(alpha_incre_firing); // the on average increase when firing
    //rate_decay ~ beta(alpha_rate_decay[1], alpha_rate_decay[2]);
    lograte_firing ~ normal(mu0_rate_firing, sigma_rate_firing);
    logitrate_decay ~ normal(mu0_rate_decay, sigma_rate_decay);
        // transition 
    theta_quiet ~ dirichlet(alpha_quiet);
    theta_firing ~ dirichlet(alpha_firing);
    theta_decay ~ dirichlet(alpha_decay);
    
    // likelihood part
    
    // HMM firing with linear increase and exponential decrease AR structure
    yd = y - trend - mu_quiet; // detrended light curve, also minus the usual mean
    

    // The forward algorithm, keep in mind we are going to start with time point 2 
    //   in stead of 1 since we have an AR(1) firing and decaying model

    // quiet state
    gamma[1,1] = normal_lpdf(yd[2]|0, sqrt(sigma2_noise));
    // firing state, a (positive tranded) Gaussian random walk 
    gamma[1,2] = exp_mod_normal_lpdf(yd[2] | yd[1] , sqrt(sigma2_noise) , rate_firing );
    // decay state, exponential decay
    gamma[1,3] = normal_lpdf(yd[2] | rate_decay * yd[1], sqrt(sigma2_noise)); 

    // then the forward algorithm will start from 3
    for(tt in 2:(N-1)){
        // at quiet state
        //  came from quiet state
        accu_quiet[1] = gamma[tt-1,1] + log(theta_quiet[1]) + 
                  normal_lpdf(yd[tt+1]|0, sqrt(sigma2_noise));
        //  came from decay state
        accu_quiet[2] = gamma[tt-1,3] + log(theta_decay[1]) + 
                  normal_lpdf(yd[tt+1]|0, sqrt(sigma2_noise));
        
        gamma[tt,1] = log_sum_exp(accu_quiet);

        // at firing state
        //  came from quiet state
        accu_firing[1] = gamma[tt-1,1] + log(theta_quiet[2]) + 
                  exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );
        //  came from firing state
        accu_firing[2] = gamma[tt-1,2] + log(theta_firing[1]) + 
                  exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );
        //  came from decay state, i.e. compound flaring
        accu_firing[3] = gamma[tt-1,3] + log(theta_decay[2]) + 
                  exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );
        gamma[tt,2] = log_sum_exp(accu_firing);  

        //  at decaying state
        //    came from firing state
        accu_decay[1] = gamma[tt-1, 2] + log(theta_firing[2]) + 
                  normal_lpdf(yd[tt+1] | rate_decay * yd[tt] , sqrt(sigma2_noise));
        //    came from decay state
        accu_decay[2] = gamma[tt-1, 3] + log(theta_decay[3]) + 
                  normal_lpdf(yd[tt+1] | rate_decay * yd[tt], sqrt(sigma2_noise));
        gamma[tt,3] = log_sum_exp(accu_decay);

    }
   target += log_sum_exp(gamma[N-1]);
}

generated quantities {
    int<lower = 1, upper = 3> viterbi[N-1]; 
    vector[N] trend;
    real log_p_state;
    
    {   
        int back_ptr[N-1, 3];
        real best_logp[N-1, 3];
        real best_total_logp;
        real logp;
        vector[N] yd;
        trend = dotCholQuasiPeriod(t, eta, B, L, P, C,diag);
        yd = y-trend-mu_quiet;

        // quiet
        best_logp[1,1] = normal_lpdf(yd[2]|0, sqrt(sigma2_noise));
        // firing
        best_logp[1,2] = exp_mod_normal_lpdf(yd[2] | yd[1] , sqrt(sigma2_noise) , rate_firing );
        // decay
        best_logp[1,3] = normal_lpdf(yd[2] | rate_decay * yd[1], sqrt(sigma2_noise)); 

        for(tt in 2:(N-1)){
            // at quiet
            best_logp[tt, 1] = negative_infinity();
            //  from quiet
            logp = best_logp[tt-1, 1] + log(theta_quiet[1]) + 
                    normal_lpdf(yd[tt+1]|0, sqrt(sigma2_noise));
            if(logp>best_logp[tt,1]){
                back_ptr[tt,1] = 1;
                best_logp[tt,1] = logp;
            }
            //  cannot be from firing
            //  from decay
            logp = best_logp[tt-1,3] + log(theta_decay[1]) + 
                    normal_lpdf(yd[tt+1]|0, sqrt(sigma2_noise));
            if(logp>best_logp[tt,1]){
                back_ptr[tt,1] = 3;
                best_logp[tt,1] = logp;
            }

            // firing state
            best_logp[tt, 2] = negative_infinity();
            // from quiet
            logp = best_logp[tt-1, 1] + log(theta_quiet[2]) + 
                    exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );
            if(logp>best_logp[tt,2]){
                back_ptr[tt,2] = 1;
                best_logp[tt,2] = logp;
            }
            // from firing
            logp = best_logp[tt-1, 2] + log(theta_firing[1]) + 
                    exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );

            if(logp>best_logp[tt,2]){
                back_ptr[tt,2] = 2;
                best_logp[tt,2] = logp;
            }
            // from decay
            logp = best_logp[tt-1, 3] + log(theta_decay[2]) + 
                    exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );

            if(logp>best_logp[tt,2]){
                back_ptr[tt,2] = 3;
                best_logp[tt,2] = logp;
            }

            // at decay
            best_logp[tt, 3] = negative_infinity();
            // cannot from quiet
            // from firing
            logp = best_logp[tt-1, 2] + log(theta_firing[2]) + 
                    normal_lpdf(yd[tt+1] | rate_decay * yd[tt] , sqrt(sigma2_noise));
                    
            if(logp>best_logp[tt,3]){
                back_ptr[tt,3] = 2;
                best_logp[tt,3] = logp;
            }

            // from decay
            logp = best_logp[tt-1, 3] + log(theta_decay[3]) + 
                    normal_lpdf(yd[tt+1] | rate_decay * yd[tt] , sqrt(sigma2_noise));
            if(logp>best_logp[tt,3]){
                back_ptr[tt,3] = 3;
                best_logp[tt,3] = logp;
            }

            
        }
        log_p_state = max(best_logp[N-1]);
        for(i in 1:3){
            if(best_logp[N-1, i]==log_p_state){
                viterbi[N-1] = i;
            }
        }

        for(tt in 2:(N-1)){
            viterbi[N-tt] = back_ptr[N-tt+1, viterbi[N-tt+1]];
        }
    }

}
