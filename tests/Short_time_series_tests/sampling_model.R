simuHMM <- function(T, mu, sigma, theta, p){
    res <- rep(NA, T)
    state <- 1 + (runif(1) > p[1])
    res[1] <- rnorm(1,mu[state],sigma[state])
    for(t in 2:T){
        state <- 1 + (runif(1) > theta[state,1])
        res[t] <- rnorm(1,mu[state],sigma[state])
    }
    return(res)
}

simuQFD <- function(N = 1000, theta_quiet = c(0.99,0.01), 
                        theta_firing = c(0.6, 0.4), 
                        theta_decay = c(0.1,0,0.95), 
                        incrm_firing = 6,
                        rate_decay = 0.85, 
                        sigma_quiet = 2, 
                        sigma_firing = 2,sigma_decay = .3,
                        init = NULL,
                        init_state = 1
                        ){

    if(is.null(init)){
        init <- rnorm(1,0,sigma_quiet)
    }
    state <- init_state
    res <- rnorm(N,0,1)
    states <- 0 * res
    res[1] <- init
    states[1] <- state
    #browser()
    for(i in 2:N){
        if(state == 1){
            state <- sample(c(1,2),1, prob = theta_quiet)# quiet only goto firing and quiet
        }
        else if (state==2) {
            state <- sample(c(2,3),1, prob = theta_firing) # firing only goto firing and decay
        }
        else {
            state <- sample(c(1,2,3),1,prob = theta_decay) # decay can go anywhere
        }
        states[i] <- state
        if(state == 1){
            res[i] <- sigma_quiet * res[i]
        }
        else if (state==2) {
            res[i] <- res[i-1] + incrm_firing + sigma_firing * res[i]
        }
        else {
            res[i] <- res[i-1] * rate_decay + sigma_decay * res[i]
        }
    }
    return(list(timeseries = res, state = states))

}

# QFD with exponentially modified normal increase
simuQFDexN <- function(N = 1000, theta_quiet = c(0.99,0.01), 
                        theta_firing = c(0.6, 0.4), 
                        theta_decay = c(0.1,0,0.95), 
                        rate_firing = 0.1,
                        rate_decay = 0.85, 
                        sigma_noise = 2, 
                        init = NULL,
                        init_state = 1
                        ){

    if(is.null(init)){
        init <- rnorm(1,0,sigma_noise)
    }
    state <- init_state
    res <- rnorm(N,0,1)
    states <- 0 * res
    res[1] <- init
    states[1] <- state
    #browser()
    for(i in 2:N){
        if(state == 1){
            state <- sample(c(1,2),1, prob = theta_quiet)# quiet only goto firing and quiet
        }
        else if (state==2) {
            state <- sample(c(2,3),1, prob = theta_firing) # firing only goto firing and decay
        }
        else {
            state <- sample(c(1,2,3),1,prob = theta_decay) # decay can go anywhere
        }
        states[i] <- state
        if(state == 1){
            res[i] <- sigma_noise * res[i]
        }
        else if (state==2) {
            res[i] <- res[i-1] + rexp(1,rate_firing) + sigma_noise * res[i]
        }
        else {
            res[i] <- res[i-1] * rate_decay + sigma_noise * res[i]
        }
    }
    return(list(timeseries = res, state = states))

}



