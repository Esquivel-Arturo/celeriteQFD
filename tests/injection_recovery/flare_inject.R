library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("./R/misc.R") # some helper
source("./R/simuFlares.R")

rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[10000:11500,c("TIME","PDCSAP_FLUX")]
rawdata <- read.csv("./Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[1:1500,c("TIME","PDCSAP_FLUX")]

rawdata[,2] <- rawdata[,2]-mean(rawdata[,2], na.rm = T)

set.seed(12345)
keplerflare_sim <- kepler_flare(rawdata[,1], .00005, 5,rPareto,xm = 50, alpha = 1, offset = 0, upper = 300) 
plot(keplerflare_sim$flare) 
#rawdata <- na.omit(rawdata)

injected <- rawdata
injected[,2] <- injected[,2] + keplerflare_sim$flare
observed <- (!is.na(injected[,2])) * 1
injected[is.na(injected[,2]),2] <- 0

N <- nrow(injected)

QFD_data <- list(N=N, t = injected[,1],
                y = injected[,2],
                observed = observed,
                sigma_prior = c(-8,8),
                Q0_prior = c(-8,8),
                dQ_prior = c(-8,8),
                period_prior = c(-8,8),
                f_prior = c(1e-6,1-1e-6),
                alpha_quiet = c(1,.1), 
                alpha_firing = c(1,1),
                alpha_decay = c(1,.1,1),
                mu0_quiet = 0,
                lambda_quiet = .01,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3,
                diag = rep(1e-6,N)
                )

modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN-missing-handling.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.9, max_treedepth=10), iter = 2000,init_r = 15, chains = 2)
summQFD <- summary(fitQFD)

tt <- rawdata[,1]


plot(tt, injected[,2], col = "blue")
lines(tt, summQFD[[1]][1:N+2*N+21, 1], type = "l")

QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]

Viterbi_max <- apply(Viterbi_raw,2,majority)

modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celerite-missing-handling.stan', 
            model_name = "celerit2", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

celeritedata <- QFD_data
celeritedata$err_prior <- c(0.01,0.01)

fitcelerite <- sampling(modelcelerite, data = celeritedata,control = list(adapt_delta = 0.9, max_treedepth=10), iter = 2000,init_r = 2, chains = 2)
summcelerite <- summary(fitcelerite)
celerite_trend <- summcelerite[[1]][1:N + (N+23),1]
residual <- injected[,2] - celerite_trend

flares3sigma <- residual >= (mean(residual) + 3 * sd(residual))

pdf("QFD_example_031381302.pdf", width = 9, height = 10)
par(mfrow = c(4,1))

plot(injected, main = "ground truth")
points(injected[which(keplerflare_sim$states==2),], col = "red",lwd=3.0)
points(injected[which(keplerflare_sim$states==3),], col = "blue",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
                lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
                cex = 1.2)

plot(injected, main = "proposed HMM")
points(injected[which(Viterbi_max==2)+1,], col = "red",lwd=3.0)
points(injected[which(Viterbi_max==3)+1,], col = "blue",lwd=3.0)
lines(injected[,1], summQFD[[1]][1:N+2*N+21, 1], col = "#d400ff",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
                lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
                cex = 1.2)

plot(injected, main = "sigma-clipping")
points(injected[flares3sigma,], col = "red",lwd=3.0)
lines(injected[,1], summcelerite[[1]][1:N + (N+23), 1], col = "#d400ff",lwd=3.0)
legend("topleft", legend = c("Potential flare","Trend"), 
                lty = c(NA,1), pch = c(1,NA), col = c("red","#d400ff"),
                cex = 1.2)



plot(injected[,1],keplerflare_sim$flare, main = "flare channel", xlab = "TIME", ylab = "FLUX")

dev.off()

save.image("../simple_flare_injection_031381302.RData")

