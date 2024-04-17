library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("./CeleriteQFD/R/misc.R") # some helper
args <- commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)

# 10 slices in total to be run on CHTC
all_slices <- list.files("./Data/TimeSeriesCSV_sliced", pattern = "csv",full.names = TRUE)
the_slice <- all_slices[args + 1]

name_slice <- paste0("./res-",args,"/", sub("./Data/TimeSeriesCSV_sliced/","" ,sub(".csv","", the_slice)))

rawdata <- read.csv(the_slice)[,c("TIME","PDCSAP_FLUX")]

rawdata[,2] <- rawdata[,2] - mean(rawdata[,2], na.rm = TRUE) # centering it
observed <- (!is.na(rawdata[,2])) * 1 # observed


rawdata[is.na(rawdata[,2]),2] <- 0
N <- nrow(rawdata)

QFD_data <- list(N=N, t = rawdata[,1],
                y = rawdata[,2],
                observed = observed,
                sigma_prior = c(-8,8),
                Q0_prior = c(-8,8),
                dQ_prior = c(-8,8),
                period_prior = c(-8,8), # kind of standard
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
                diag = rep(1e-6,N),
                err_prior <- c(0.01,0.01)
                )

modelQFD <- stan_model(file = './CeleriteQFD/Stan/Morphology/QFD/CeleriteQFDexN-missing-handling.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'CeleriteQFD/celerite2/celerite2.hpp'), '"\n'))

fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000,init_r = 15, chains = 2, thin = 2)

## get trend
QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]
QFD_trend_raw <- QFD_samples[,1:N+2*N+21]

#write.csv(QFD_samples,paste0(name_slice,"-QFD-samples.csv"))
write.csv(QFD_trend_raw,paste0(name_slice,"-QFD-trend.csv"))
write.csv(Viterbi_raw,paste0(name_slice,"-QFD-Viterbi.csv"))