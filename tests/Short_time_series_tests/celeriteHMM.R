library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[16400:17400,c("TIME","PDCSAP_FLUX")]
rawdata <- na.omit(rawdata)
rawdata[,2] <- rawdata[,2] - mean(rawdata[,2])
N <- nrow(rawdata)
plot(rawdata)
tt = 50 * (rawdata[,1] - min(rawdata[,1]))/(range(rawdata[,1])[2]-range(rawdata[,1])[1]) # sort of normalize the time to avoid super short period

star_data <- list(N=N, t = tt, y = rawdata[,2],
                     sigma_prior = c(-2,2),
                     period_prior = c(-1,3),
                     Q0_prior = c(-2,2),
                     dQ_prior = c(-2,2),
                     f_prior = c(1e-6,1-1e-6),
                     alpha = c(1,1), # transition
                     mu0 = 0, # usual mean
                     lambda = 0.001, # prior for usual
                     noise_prior = c(0.01,0.01),
                     rate = 0.0001, 
                     diag = rep(0,N),
                     eps_neg = 1e-2)

modellaplace <- stan_model(file = './Stan/Prototypes/CeleriteHMM/celerite_HMM_laplace.stan', 
            model_name = "celeritHMMlaplace", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fitlaplace <- sampling(modellaplace, data = star_data,
    control = list(adapt_delta = 0.99, max_treedepth=15), 
    iter = 4000)

save.image("./Res/131799991_16400-17400/res.RData")

summ_fitlaplace <- summary(fitlaplace)
states_est2 <- summ_fitlaplace[[1]][1:N + 1016,1]

pdf("Res/131799991_16400-17400/det.pdf", width = 10, height = 4)
par(mfrow = c(1,2))

firing_est <- which( states_est2 >=1.5) 
plot(rawdata, main = "celeriteHMM")
points(rawdata[firing_est,], col = "red")
detrended <- summ_fitlaplace[[1]][1:N + 3016,1]

firing_3sigma <- which(abs(detrended-mean(detrended))>= 3 * sd(detrended))
plot(rawdata, main = "3-sigma")
points(rawdata[firing_3sigma,], col = "red")
dev.off()
