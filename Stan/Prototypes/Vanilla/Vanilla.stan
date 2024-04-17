data {
    int<lower = 0> T;
    int<lower = 1> K; // number of components
    vector[T] y;
    vector[T] observed;// is 1 if observed and 0 if missing
    vector<lower = 0>[K] alpha;
    // real<lower = 0, upper = 1> pinit; //prior on initial state
    real mu0;
    real lambda;
    real shape;
    real rate;// normal inverse gamma
}


parameters {
    ordered[K] mu; // mean for two states
    real<lower = 0> sigma[K];// std for two states
    simplex[K] theta[K]; // transition matrix
}

model {
    real accu[K];
    real gamma[T,K];// joint likelihood of first t states
    for(i in 1:K){
        gamma[1,i] = normal_lpdf(y[1]|mu[i], sigma[i]);
    }
    for(t in 2:T){
        for(i in 1:K){
            // time t state is i
            for(j in 1:K){
                // and came from j
                accu[j] = gamma[t-1, j] + log(theta[j,i]) + 
                        observed[t] * normal_lpdf(y[t] | mu[i], sigma[i]);// not observed then don't add likelihood
            }
            gamma[t,i] = log_sum_exp(accu);
        }
    }

    for(i in 1:K){
        theta[i] ~ dirichlet(alpha);
        mu[i] ~ normal(mu0, sigma[i]/lambda);
        sigma[i] ~ inv_gamma(shape,rate);
    }
    target += log_sum_exp(gamma[T]);
}

// Viterbi
generated quantities {
    int <lower = 1, upper = K>state[T];
    real log_p_state;
    {
        int back_ptr[T, K];
        real best_logp[T, K];
        real best_total_logp;
        for(i in 1:K){
            best_logp[1,i] = normal_lpdf(y[1]|mu[i], sigma[i]);
        }
        for(t in 2:T){
            for(i in 1:K){
                best_logp[t, i] = negative_infinity();
                for(j in 1:K){
                    real logp;
                    logp = best_logp[t-1,j]  + log(theta[j,i]) + 
                        observed[t] * normal_lpdf(y[t] | mu[i], sigma[i]);
                    if(logp > best_logp[t,i]){
                        back_ptr[t,i] = j;
                        best_logp[t,i] = logp;
                    }
                }
            }
        }
        log_p_state = max(best_logp[T]);
        for(i in 1:K){
            if(best_logp[T, i]==log_p_state){
                state[T] = i;
            }
        }
        for(t in 1:(T-1)){
            state[T-t] = back_ptr[T-t+1, state[T-t+1]];
        }
    }
}

