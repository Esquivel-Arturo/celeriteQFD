library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rawdata <- read.csv(list.files("./Data", full.names = TRUE)[1])

#timeseries <- rawdata$PDCSAP_FLUX[-1] # first entry is NA...
timeseries <- rawdata$PDCSAP_FLUX[6000:8000]
observed <- 1 * (!is.na(timeseries))

timeseries[is.na(timeseries)] <- timeseries[1]
plot(timeseries)


star_stan <- list(T=length(timeseries),
                    K=2, 
                    y = timeseries,
                    observed = observed, 
                    alpha = c(1,1),
                    mu0=0,
                    lambda = 1e-10,
                    shape = 0.001,
                    rate = 0.001)

fit <- stan(file = './Stan/Vanilla.stan',model_name = "star_vanilla",data = star_stan, 
                iter = 3000, control = list(adapt_delta = .8))
plot(fit)

# for Win using cmdstan
library(cmdstanr)

file <- './Stan/Vanilla.stan'
mod <- cmdstan_model(file)

fit <- mod$sample(
  data = star_stan,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  max_treedepth = 10
)

fit$save_object(file = "./Stan/vanilla_star.RDS")

