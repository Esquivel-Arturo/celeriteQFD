library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

source("./R/sampling_model.R")
N <- 300
set.seed(12345)
res <- simuQFDexN(N = N, rate_firing = .1, 
                theta_quiet = c(0.95,0.05),
                theta_firing = c(0.15, 0.85),
                theta_decay = c(0.2,0,0.8),rate_decay = 0.6, 
                sigma_noise = 1)

tt <- seq(0,10,length.out = N)
f <- 5*sin(3 * tt)
plot(tt,f)
plot(res$timeseries)
simu_signal <- f+res$timeseries
plot(simu_signal)
points(which(res$state==2),simu_signal[res$state==2], col = "red")
points(which(res$state==3),simu_signal[res$state==3], col = "blue")



QFDsimple_data <- list(N=N,
                y = res$timeseries,
                alpha_quiet = c(1,1),
                alpha_firing = c(1,1),
                alpha_decay = c(1,1,1),
                mu0_quiet = 0,
                lambda_quiet = 1e-4,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3
                )


modelQFDsimple <- stan_model(file = './Stan/Morphology/QFD/QFDexN.stan', 
            model_name = "QFTexN")

optQFDsimple <- optimizing(modelQFDsimple, QFDsimple_data)

par(mfrow = c(1,2))
plot(res$timeseries)
points(which(res$state==2),res$timeseries[res$state==2], col = "red")
points(which(res$state==3),res$timeseries[res$state==3], col = "blue")
temp <- res$timeseries[-1]
viterbi <- optQFDsimple$par[1:(N-1) + 13]
plot(temp)
points(which(viterbi==2),temp[viterbi==2], col = "red")
points(which(viterbi==3),temp[viterbi==3], col = "blue")

vbQFDsimple <- vb(modelQFDsimple, data = QFDsimple_data,tol_rel_obj = 0.005)
fitQFDsimple <- sampling(modelQFDsimple, data = QFDsimple_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000)



QFD_data <- list(N=N, t = tt,
                y = simu_signal,
                sigma_prior = c(-2,5),
                #Q0_prior = c(2,4),
                Q0_prior = c(-2,4),# this is key, we need to set quality to be not too small
                dQ_prior = c(-2,4),
                period_prior = c(-3,3),
                #period_prior = c(0,3),
                f_prior = c(1e-6,1-1e-6),
                alpha_quiet = c(1,1), 
                alpha_firing = c(1,1),
                alpha_decay = c(1,1,1),
                mu0_quiet = 0,
                lambda_quiet = .01,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3,
                diag = rep(1e-6,N),
                err_prior = c(0.01, 0.01) # just for celerite, same as gamma_noise
                )

modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celerite.stan', 
            model_name = "celerit", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

vbcelerite <- vb(modelcelerite, data = QFD_data, tol_rel_obj = 0.0001)
fitcelerite <- sampling(modelcelerite, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=10), iter = 2000)
summ_celerite <- summary(fitcelerite)
summ_celerite[[1]][1:23 + N,1]
plot(summ_celerite[[1]][1:N + (N+23),1], type = "l")
lines(f, col = "green")
points(simu_signal)
points(which(res$state==2),simu_signal[res$state==2], col = "red")
points(which(res$state==3),simu_signal[res$state==3], col = "blue")

plot(simu_signal - summ_celerite[[1]][1:N + (N+23),1])

vbQFD <- vb(modelQFD, data = QFD_data,tol_rel_obj = 0.0001, iter = 5e4)
summ_vbQFD <- summary(vbQFD)
# state
par(mfrow = c(3,1))
plot(summ_vbQFD[[1]][1:(N-1) + (N + 22),1])
plot(res$state)
plot(simu_signal[-1])
(summ_vbQFD[[1]][1:22 + N,1])

par(mfrow = c(2,1))
plot(summ_vbQFD[[1]][1:N+2*N+21, 1], type = "l", main = "QFD")
points(simu_signal, col = "blue")
plot(summ_celerite[[1]][1:N + (N+23),1], type = "l")
points(simu_signal,col = "blue")

fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000)
summQFD <- summary(fitQFD)
summQFD[[1]][1:22 + N,1] 

par(mfrow = c(2,1))
plot(summQFD[[1]][1:N+2*N+21, 1], type = "l", main = "QFD")
points(simu_signal, col = "blue")
plot(summ_celerite[[1]][1:N + (N+23),1], type = "l")
points(simu_signal,col = "blue")

par(mfrow = c(3,1))
plot(summQFD[[1]][1:(N-1) + (N + 22),1])
plot(res$state)
plot(simu_signal[-1])
(summQFD[[1]][1:22 + N,1])
