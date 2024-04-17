library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

source("./R/sampling_model.R")
N <- 300
set.seed(12345)
res <- simuQFDexN(N = N, rate_firing = .1, 
                theta_quiet = c(0.95,0.05),
                theta_firing = c(0.2, 0.8),
                theta_decay = c(0.2,0,0.8),rate_decay = 0.6, 
                sigma_noise = 1)

tt <- seq(0,(50/300)*N,length.out = N)
f <- 5 * sin(2*tt/(2*pi))
plot(tt,f)
plot(res$timeseries)
simu_signal <- f+res$timeseries
plot(simu_signal)
points(which(res$state==2),simu_signal[res$state==2], col = "red")
points(which(res$state==3),simu_signal[res$state==3], col = "blue")



QFD_data <- list(N=N, t = tt,
                y = simu_signal,
                S0_prior = c(0,5),
                w0_prior = c(-5,5),
                Q_prior = c(3,10),
                alpha_quiet = c(100,1),
                alpha_firing = c(1,1),
                alpha_decay = c(1,1,1),
                mu0_quiet = 0,
                lambda_quiet = 1e-4,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3,
                diag = rep(0,N),
                err_prior = c(0.01, 0.01) # just for celerite, same as gamma_noise
                )

modelSHOQFD <- stan_model(file = './Stan/Morphology/QFD/SHOQFDexN.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))



vbSHOQFD <- vb(modelSHOQFD, data = QFD_data,tol_rel_obj = 0.001)
summ_vbSHOQFD <- summary(vbSHOQFD)
plot(summ_vbQFD[[1]][1:1000 + 1020,1])
summ_vbSHOQFD[[1]][1001:1020,] 
fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=10), iter = 2000)
summQFD <- summary(fitQFD)
summQFD[[1]][1001:1020,] 

fitQFD_array <- as.array(fitQFD)
hist(fitQFD_array[,1,1006])



