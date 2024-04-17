source("./R/sampling_model.R")
N <- 1000
set.seed(12345)
res <- simuQFD(N = N, incrm_firing = 6, 
                theta_quiet = c(0.99,0.01),
                theta_firing = c(0.6, 0.4),
                theta_decay = c(0.1,0,0.9),rate_decay = 0.85, 
                sigma_decay = 2, 
                sigma_firing = 2, sigma_quiet = 2)
plot(res$timeseries)
points(which(res$state==2),res$timeseries[res$state==2], col = "red")
points(which(res$state==3),res$timeseries[res$state==3], col = "blue")

tt <- seq(0,50,length.out = N)
f <- 5 * sin(4*tt/(2*pi))

simu_signal <- f+res$timeseries

plot(simu_signal)
points(which(res$state==2),simu_signal[res$state==2], col = "red")
points(which(res$state==3),simu_signal[res$state==3], col = "blue")

res2 <- simuQFDeN(rate_firing = .3, sigma_noise = 1)
plot(res2$timeseries)
points(which(res2$state==2),res2$timeseries[res2$state==2], col = "red")
points(which(res2$state==3),res2$timeseries[res2$state==3], col = "blue")
