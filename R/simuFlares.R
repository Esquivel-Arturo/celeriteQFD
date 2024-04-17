rPareto <- function(n,xm=1,alpha=1,offset=0,upper = Inf){
    n_sample <- 0
    res <- c()
    while(n_sample<n){
        n_remain <- n-n_sample
        U <- runif(n_remain)
        x_temp <- xm/((U)^(1/alpha))+offset
        x_temp <- x_temp[x_temp<=upper]
        res <- c(res,x_temp)
        n_sample <- length(res)
    }
    return(res)
}

rmovexp <- function(n,rate = 1,offset = 0){
    offset+rexp(n,rate)
}

kepler_raising <- function(x){
    1+1.941 * x  - 0.175 * x^2-2.246 * x^3 -1.125 * x^4
}

kepler_decay <- function(x){
    0.6890 * exp(-1.6 * x) + 0.3030 * exp(-0.2783 * x)
}


kepler_flare <- function(tt,t_half,n,flux_dist = rPareto, ...){
    flare <- 0 * tt
    states <- 1 + flare # store states
    flux_all <- flux_dist(n = n, ...)
    peak_time_all <- sample(tt, n)
    peak_time_all <- sort(peak_time_all)
    #browser()
    for(i in 1:n){
        flux_loc <- flux_all[i]
        t_half_loc <- t_half * flux_loc
        peak_time_loc <- peak_time_all[i]
        
        raising_phase <- which(peak_time_loc-tt >=0 & peak_time_loc-tt <= t_half_loc )
        decaying_phase <- which(tt - peak_time_loc >0 & tt - peak_time_loc <= 10 * t_half_loc )

        states[raising_phase] <- 2
        states[decaying_phase] <- 3

        raising_time_loc <- (tt[raising_phase]-peak_time_loc)/t_half_loc
        decaying_time_loc <- (tt[decaying_phase]-peak_time_loc)/t_half_loc

        raising_flux <- flux_loc * kepler_raising(raising_time_loc)
        decaying_flux <- flux_loc * kepler_decay(decaying_time_loc)
        #browser()
        flare[raising_phase] <- raising_flux
        flare[decaying_phase] <- decaying_flux
    }

    return(list(flare = flare, states = states))

}
