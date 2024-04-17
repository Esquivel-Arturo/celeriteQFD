library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[16000:18000,c("TIME","PDCSAP_FLUX")]
rawdata <- na.omit(rawdata)
rawdata[,2] <- rawdata[,2] - mean(rawdata[,2])

celerite_data <- list(N=nrow(rawdata), t = rawdata[,1], y = rawdata[,2],
                     sigma_prior = c(-2,2),
                     period_prior = c(-1,3),
                     Q0_prior = c(-2,2),
                     dQ_prior = c(-2,2),
                     f_prior = c(1e-6,1-1e-6),
                     err_prior = c(0.01,0.01),
                     diag = 0*rawdata$TIME)


mod <- stan_model(file = "./Stan/celerite.stan", model_name = "celerit", 
            allow_undefined = TRUE,
           includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))


fit <- sampling(mod, data = celerite_data )
#fit2 <- optimizing(mod, data = celerite_data, verbose = T, iter = 10000) # should not be used for hidden variable parameterization 
plot(fit2$par[1:1605], type = "l")
plot(rawdata[,2]-fit2$par[1:1605])
summ_fit <- summary(fit)
