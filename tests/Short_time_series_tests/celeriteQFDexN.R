library(rstan)
library(dplyr)
library(ggpubr)
library(ggplot2)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)
source("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/R/misc.R") # some helper

# run QFD
#rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[16000:17400,c("TIME","PDCSAP_FLUX")]
rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[6000:8000,c("TIME","PDCSAP_FLUX")]

# handling missing data is currently not implemented, but only little are missing so I will just omit it for now

rawdata[,2] <- rawdata[,2] - mean(rawdata[,2], na.rm = TRUE)
observed <- (!is.na(rawdata[,2])) * 1
rawdata[is.na(rawdata[,2]),2] <- 0
N <- nrow(rawdata)
plot(rawdata)
tt <- rawdata[,1]

QFD_data <- list(N=N, t = rawdata[,1],
                y = rawdata[,2],
                observed = observed,
                sigma_prior = c(-8,8),
                Q0_prior = c(0,8),
                #Q0_prior = c(2,4), # Higher Q -> higher wave amplitude
                #Q0_prior = c(-8,8),# this is key, we need to set quality to be not too small
                dQ_prior = c(-8,8),
                #period_prior = c(-3,3),
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

modelQFD <- stan_model(file = '/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Stan/Morphology/QFD/CeleriteQFDexN-missing-handling.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path('/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/', 
                             'celerite2/celerite2.hpp'), '"\n'))

fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000,init_r = 15, chains = 2)
summQFD <- summary(fitQFD)

plot(tt, rawdata[,2], col = "blue")
lines(tt, summQFD[[1]][1:N+2*N+21, 1], type = "l")

par(mfrow = c(2,1))
plot(tt[-1], summQFD[[1]][1:(N-1) + (N + 22),1])
plot(tt[-1], rawdata[-1,2])
lines(tt, summQFD[[1]][1:N+2*N+21, 1], type = "l", col = "red")

plot(tt, rawdata[,2]-summQFD[[1]][1:N+2*N+21, 1])

residual <- rawdata[,2]-summQFD[[1]][1:N+2*N+21, 1]
plot(tt, residual)
flare3s <- which(residual >= (mean(residual) + 3 * sd(residual)))
points(tt[flare3s], residual[flare3s], col = "red")

#save.image("../res164-174-cQFDexN.RData")
load("../res164-174-cQFDexN.RData")

## visualize
QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]

Viterbi_max <- apply(Viterbi_raw,2,majority)

pdf("./Res/CeleriteQFD/131799991_16400-17400/det.pdf", width = 10, height = 6)
plot(rawdata)
lines(rawdata[,1], summQFD[[1]][1:N+2*N+21, 1], col = "#d400ff",lwd=3.0)
points(rawdata[which(Viterbi_max==2)+1,], col = "red",lwd=3.0)
points(rawdata[which(Viterbi_max==3)+1,], col = "blue",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
                lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
                cex = .75)
dev.off()

# celerite along

modelcelerite <- stan_model(file = '/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Stan/Prototypes/Celerite/celerite-missing-handling.stan', 
            model_name = "celerit2", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path('/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/', 
                             'celerite2/celerite2.hpp'), '"\n'))

celeritedata <- QFD_data
celeritedata$err_prior <- c(0.01,0.01)

fitcelerite <- sampling(modelcelerite, data = celeritedata,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000,init_r = 2, chains = 2)
summcelerite <- summary(fitcelerite)
celerite_trend <- summcelerite[[1]][1:N + (N+23),1]
residual <- rawdata[,2] - celerite_trend

flares3sigma <- residual >= (mean(residual) + 3 * sd(residual))

##### Fig.6 in the paper #####
save.image("../res031381302_6000-8000-cQFDexN.RData")
load("../res031381302_6000-8000-cQFDexN.RData")


#png("./Res/CeleriteQFD/031381302_6000-8000/det_compare2.png", width = 7, height = 5.5, units ="in" , res = 500)
pdf("./Res/CeleriteQFD/031381302_6000-8000/det_compare2.pdf", width = 7, height = 5.5)
par(mfrow = c(2,2),mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
par(fig = c(0,7,5,10)/10)
plot(rawdata, main = "Proposed", ylab = "Centered Flux", xlab = "Time", pch = 3)
points(rawdata[which(Viterbi_max==2)+1,], col = "red",lwd=2.0, pch = 3)
points(rawdata[which(Viterbi_max==3)+1,], col = "blue",lwd=2.0, pch = 3)
lines(rawdata[,1], summQFD[[1]][1:N+2*N+21, 1], col = "#d400ff",lwd=3.0)

rect(xleft = 1334.9, ybottom = -30, xright = 1335.25, ytop = 60, border = "orange", lwd = 2)
#rect(xleft = 1336.05, ybottom = -30, xright = 1336.3, ytop = 60, border = "orange", lwd = 2)
highlight_range <- rawdata[,1] >= 1334.9 & rawdata[,1] <= 1335.25
legend("topleft", legend = c("Firing","Decay","Trend"), 
                lty = c(NA,NA,1), pch = c(3,3,NA), col = c("red","blue","#d400ff"),
                cex = .75)

par(fig = c(7,10,5,10)/10, mar = c(3,0,2,1), new = T)
plot(rawdata[highlight_range, ], ylab = "", xlab = "Time", pch = 3)
points(rawdata[which(Viterbi_max==2)+1,], col = "red",lwd=2.0, pch = 3)
points(rawdata[which(Viterbi_max==3)+1,], col = "blue",lwd=2.0, pch = 3)
lines(rawdata[highlight_range,1], summQFD[[1]][1:N+2*N+21, 1][highlight_range], col = "#d400ff",lwd=3.0)


par(fig = c(0,7,0,5)/10, mar = c(3,3,2,2), new = T)
plot(rawdata, main = "sigma-clipping", ylab = "Centered Flux", xlab = "Time", pch = 3)
points(rawdata[flares3sigma,], col = "red",lwd=2.0, pch = 3)
lines(rawdata[,1], summcelerite[[1]][1:N + (N+23), 1], col = "#d400ff",lwd=3.0)
rect(xleft = 1334.9, ybottom = -30, xright = 1335.25, ytop = 60, border = "orange", lwd = 2)
#rect(xleft = 1336.05, ybottom = -30, xright = 1336.3, ytop = 60, border = "orange", lwd = 2)

legend("topleft", legend = c("Potential Flares","Trend"), 
                lty = c(NA,1), pch = c(3,NA), col = c("red","#d400ff"),
                cex = .75)

par(fig = c(7,10,0,5)/10, mar = c(3,0,2,1), new = T)
plot(rawdata[highlight_range, ], ylab = "", xlab = "Time", pch = 3)
points(rawdata[flares3sigma,], col = "red",lwd=2.0, pch = 3)
lines(rawdata[highlight_range,1], summcelerite[[1]][1:N + (N+23), 1][highlight_range], 
      col = "#d400ff",lwd=3.0)
dev.off()

detrended_data <- data.frame(TIME = rawdata$TIME, detrended = rawdata[,2]-summQFD[[1]][1:N+2*N+21, 1])

pdf("./Res/CeleriteQFD/131799991_16400-17400/det_flares.pdf", width = 10, height = 10)
par(mfrow = c(2,1))
plot(detrended_data[90:140,],type = "l")
points(detrended_data[90:140,])
points(detrended_data[110,],col = "red")
points(detrended_data[111:120,], col = "blue")
plot(detrended_data[500:550,],type = "l") 
points(detrended_data[500:550,])
points(detrended_data[c(520,527),],col = "red")
points(detrended_data[c(521:526,528:533),], col = "blue")
dev.off()


# Producir sensitivity y specificity plots
## Phil's input: duration, max energy, integrated energy 
# Queremos detrended curve (y_d)
# Trend 
# Viterbi sequences 
# Variables relevantes para los plots 
# Ajustar toda la estrella y guardar el output 


highlight_range <- rawdata[,1] >= 1335.05 & rawdata[,1] <= 1335.15
flare_det <- detrended_data[highlight_range,]
flare_vit <- Viterbi_raw[, highlight_range[2:2001]]
flare_vit_max <- apply(flare_vit,2,majority)

pdf("/Users/arturoesquivel/Downloads/flares.pdf", width = 10, height = 6)
par(mfrow = c(2,1))
plot(flare_det, ylab = "Brightness", xlab = "Time")
points(flare_det[which(flare_vit_max==2),], col = "red",lwd=2.0)
points(flare_det[which(flare_vit_max==3),], col = "blue",lwd=2.0)
 
   
flare_vit_probs <- matrix(0, 3, ncol(flare_vit))
for(i in 1:ncol(flare_vit)){
  flare_vit_probs[1,i] <- sum(ifelse(flare_vit[,i] == 1, 1, 0))/nrow(flare_vit)
  flare_vit_probs[2,i] <- sum(ifelse(flare_vit[,i] == 2, 1, 0))/nrow(flare_vit)
  flare_vit_probs[3,i] <- sum(ifelse(flare_vit[,i] == 3, 1, 0))/nrow(flare_vit)
}

barplot(flare_vit_probs, 
        col=c("black","red","blue") , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="State Probabilities")


highlight_range2 <- rawdata[,1] >= 1335.4 & rawdata[,1] <= 1335.47
flare_det2 <- detrended_data[highlight_range2,]
flare_vit2 <- Viterbi_raw[, highlight_range2[2:2001]]
flare_vit_max2 <- apply(flare_vit2,2,majority)
par(mfrow = c(2,1))
plot(flare_det2, ylab = "Brightness", xlab = "Time")
points(flare_det2[which(flare_vit_max2==2),], col = "red",lwd=2.0)
points(flare_det2[which(flare_vit_max2==3),], col = "blue",lwd=2.0)

flare_vit_probs2 <- matrix(0, 3, ncol(flare_vit2))
for(i in 1:ncol(flare_vit2)){
  flare_vit_probs2[1,i] <- sum(ifelse(flare_vit2[,i] == 1, 1, 0))/nrow(flare_vit2)
  flare_vit_probs2[2,i] <- sum(ifelse(flare_vit2[,i] == 2, 1, 0))/nrow(flare_vit2)
  flare_vit_probs2[3,i] <- sum(ifelse(flare_vit2[,i] == 3, 1, 0))/nrow(flare_vit2)
}

barplot(flare_vit_probs2, 
        col=c("black","red","blue") ,
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="State Probabilities")


highlight_range3 <- rawdata[,1] >= 1336.1 & rawdata[,1] <= 1336.17
flare_det3 <- detrended_data[highlight_range3,]
flare_vit3 <- Viterbi_raw[, highlight_range3[2:2001]]
flare_vit_max3 <- apply(flare_vit3,2,majority)
par(mfrow = c(2,1))
plot(flare_det3, ylab = "Brightness", xlab = "Time")
points(flare_det3[which(flare_vit_max3==2),], col = "red",lwd=2.0)
points(flare_det3[which(flare_vit_max3==3),], col = "blue",lwd=2.0)
 
flare_vit_probs3 <- matrix(0, 3, ncol(flare_vit3))
for(i in 1:ncol(flare_vit3)){
  flare_vit_probs3[1,i] <- sum(ifelse(flare_vit3[,i] == 1, 1, 0))/nrow(flare_vit3)
  flare_vit_probs3[2,i] <- sum(ifelse(flare_vit3[,i] == 2, 1, 0))/nrow(flare_vit3)
  flare_vit_probs3[3,i] <- sum(ifelse(flare_vit3[,i] == 3, 1, 0))/nrow(flare_vit3)
  }

barplot(flare_vit_probs3, 
        col=c("black","red","blue") , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="State Probabilities")
dev.off()

## Check time series for the divergent chunk 
### Compare parameters between chunks
## Try to compute relative energy 
## Distribution of relative energies and durations
## Integrated energy? Max energy? 

load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_12551-15060.rda")
summQFD <- summary(fitQFD)
# summQFD[[1]][(3*N + 17):(3*N + 36),]

tpm <- matrix(0,3,3)

tpm[1,1:2] <- summQFD[[1]][(N+6):(N + 7),1]
tpm[2,2:3] <- summQFD[[1]][(N+10):(N + 11),1]
tpm[3,1:3] <- summQFD[[1]][(N+13):(N + 15),1]

colnames(tpm) <- c("Quiet", "Firing", "Decay")
rownames(tpm) <- c("Quiet", "Firing", "Decay")
tpm_heat <- melt(tpm)
tpm_heat$interval <- rep(n,9)
all_tpm <- rbind(all_tpm, tpm_heat)

colnames(all_tpm) <- c("from_state", "to_state", "probability", "interval")
ggplot(data = data.frame(all_tpm), aes(x = to_state, y = from_state, fill= probability)) + 
  geom_tile() +
  facet_grid(interval~., scales="free_y")


rstan::traceplot(fitQFD, pars = c("theta_quiet", "theta_firing", "theta_decay"))
rstan::traceplot(fitQFD, pars = c("rho1", "w2", "tau2"))
rstan::traceplot(fitQFD, pars = c("trend[1]","trend[50]", "trend[600]", "trend[1400]", "trend[1990]"))
rstan::traceplot(fitQFD, pars = c("lsigma", "lperiod", "mu_quiet", "sigma", "rate_decay"))
rstan::traceplot(fitQFD, pars = c("eta[11]","eta[150]", "eta[400]", "eta[1403]", "eta[1992]"))


flares <- NULL
durations <- NULL
bright <- NULL

rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[12551:15060,c("TIME","PDCSAP_FLUX")]

QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]
Viterbi_max <- apply(Viterbi_raw,2,majority)

count <- 1
for(i in Viterbi_max[2:length(Viterbi_max)]){
  count <- count + 1
  if(i == 2 & Viterbi_max[count-1] == 1){
    flare_start <- rawdata$TIME[count]
    flares <- append(flares, flare_start)
    bright <- append(bright, rawdata$PDCSAP_FLUX[count])
  }
  if(i == 3 & Viterbi_max[count+1] == 1){
    flare_end <- rawdata$TIME[count]
    durations <- append(durations, flare_end - flare_start)
  }
}


firing_rates <- NULL
firing_rates_lb <- NULL
firing_rates_ub <- NULL

load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/sum_2-2510.rda")
summQFD1 <- summQFD
N <- 2509

firing_rates <- append(firing_rates, summQFD1[[1]][N + 21, 1])
firing_rates_lb <- append(firing_rates_lb, summQFD1[[1]][N + 21, 4])
firing_rates_ub <- append(firing_rates_ub, summQFD1[[1]][N + 21, 8])


load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_2511-5020.rda")
rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[2511:5020,c("TIME","PDCSAP_FLUX")]
N <- 2510
summQFD2 <- summary(fitQFD)

firing_rates <- append(firing_rates, summQFD2[[1]][N + 21, 1])
firing_rates_lb <- append(firing_rates_lb, summQFD2[[1]][N + 21, 4])
firing_rates_ub <- append(firing_rates_ub, summQFD2[[1]][N + 21, 8])


load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/sum_5021-7530.rda")
rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[5021:7530,c("TIME","PDCSAP_FLUX")]
N <- 2510
summQFD3 <- summQFD

firing_rates <- append(firing_rates, summQFD3[[1]][N + 21, 1])
firing_rates_lb <- append(firing_rates_lb, summQFD3[[1]][N + 21, 4])
firing_rates_ub <- append(firing_rates_ub, summQFD3[[1]][N + 21, 8])


load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_7531-9524.rda")
summQFD4 <- summary(fitQFD)
N=1994
rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[7531:9524,c("TIME","PDCSAP_FLUX")]

firing_rates <- append(firing_rates, summQFD4[[1]][N + 21, 1])
firing_rates_lb <- append(firing_rates_lb, summQFD4[[1]][N + 21, 4])
firing_rates_ub <- append(firing_rates_ub, summQFD4[[1]][N + 21, 8])


load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_10341-12550.rda")
N=2210
summQFD5 <- summary(fitQFD)
rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[10341:12550,c("TIME","PDCSAP_FLUX")]

firing_rates <- append(firing_rates, summQFD5[[1]][N + 21, 1])
firing_rates_lb <- append(firing_rates_lb, summQFD5[[1]][N + 21, 4])
firing_rates_ub <- append(firing_rates_ub, summQFD5[[1]][N + 21, 8])


load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_12551-15060.rda")
N=2510
summQFD6 <- summary(fitQFD)
rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[12551:15060,c("TIME","PDCSAP_FLUX")]

firing_rates <- append(firing_rates, summQFD6[[1]][N + 21, 1])
firing_rates_lb <- append(firing_rates_lb, summQFD6[[1]][N + 21, 4])
firing_rates_ub <- append(firing_rates_ub, summQFD6[[1]][N + 21, 8])


load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_17571-20076.rda")
summQFD7 <- summary(fitQFD)
N=2506

firing_rates <- append(firing_rates, summQFD7[[1]][N + 21, 6])
firing_rates_lb <- append(firing_rates_lb, summQFD7[[1]][N + 21, 4])
firing_rates_ub <- append(firing_rates_ub, summQFD7[[1]][N + 21, 8])


firing_rates_df <- data.frame("est" = firing_rates,
                              "lb" = firing_rates_lb,
                              "ub" = firing_rates_ub,
                              "interval" = as.factor(seq(1,7,1)))
ggplot(firing_rates_df, aes(x=interval, y=est, group=interval, color=interval)) + 
  geom_point()+
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2,
                position=position_dodge(0.05)) +
  labs(y = "Firing Rate Flux Increase Estimate")



#rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[16000:17400,c("TIME","PDCSAP_FLUX")]
rawdata <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[ ,c("TIME","PDCSAP_FLUX")]

# handling missing data is currently not implemented, but only little are missing so I will just omit it for now

rawdata[,2] <- rawdata[,2] - mean(rawdata[,2], na.rm = TRUE)
observed <- (!is.na(rawdata[,2])) * 1
rawdata[is.na(rawdata[,2]),2] <- 0
N <- nrow(rawdata)
plot(rawdata)

pdf("/Users/arturoesquivel/Downloads/interval_comparison.pdf", width = 10, height = 6)
ggplot(data = data.frame(all_tpm), aes(x = to_state, y = from_state, fill= probability)) + 
  geom_tile() +
  facet_grid(interval~., scales="free_y")

ggplot(firing_rates_df, aes(x=interval, y=est, group=interval, color=interval)) + 
  geom_point()+
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2,
                position=position_dodge(0.05)) +
  labs(y = "Firing Rate Flux Increase Estimate")
plot(rawdata)
dev.off()

rawdata$Flare <- as.factor(rawdata$firing)
ggplot(rawdata, aes(x=Flare, y=PDCSAP_FLUX, fill=Flare )) +
  geom_boxplot()

for(i in 1:nrow(rawdata)){
  for(j in 1:length(flares)){
    if(!is.na(rawdata$TIME[i])){
    if(rawdata$TIME[i] >= flares[j] & rawdata$TIME[i] <= flares[j] + durations[j]){
      rawdata$firing[i] <- 1
    }
  }
  }
}


## trapezoidal sum to integrate the energy of the flare, calculate it in relative energy (divide by the median)
## re-fit the model 


load(file="/Users/arturoesquivel/Downloads/model_1.rda")

data <- read.csv("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[ ,c("TIME","PDCSAP_FLUX")]

observed <- (!is.na(data[,2])) * 1
med <- median(data[,2], na.rm = TRUE)
data[is.na(data[,2]),2] <- med
data[,2] <- data[,2]/med

rawdata <- data[2:2510, ]
N <- nrow(rawdata)

rstan::traceplot(fitQFD, pars = c("theta_quiet", "theta_firing", "theta_decay"))
rstan::traceplot(fitQFD, pars = c("rho1", "w2", "tau2"))
rstan::traceplot(fitQFD, pars = c("trend[1]","trend[50]", "trend[600]", "trend[1400]", "trend[1990]"))
rstan::traceplot(fitQFD, pars = c("lsigma", "lperiod", "mu_quiet", "sigma", "rate_decay"))
rstan::traceplot(fitQFD, pars = c("eta[11]","eta[150]", "eta[400]", "eta[1403]", "eta[1992]"))
rstan::traceplot(fitQFD, pars = c("viterbi[1784]","viterbi[1790]", "viterbi[1773]"))


summQFD <- summary(fitQFD)
## visualize
QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]

Viterbi_max <- apply(Viterbi_raw,2,majority)

plot(rawdata)
lines(rawdata[,1], summQFD[[1]][1:N+2*N+21, 1]/med+1, col = "#d400ff",lwd=3.0)
points(rawdata[which(Viterbi_max==2)+1,], col = "red",lwd=3.0)
points(rawdata[which(Viterbi_max==3)+1,], col = "blue",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
       lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
       cex = .75)

# plug the one into the trend to see if we can change it from being centered around zero
# within celerite 


## trapezoidal integration for the duration of the star in relative energy
## for each of the posterior samples, distribution of the energy of each flare across posterior samples
## and distribution of energy across the flares (also posterior x flare)

load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_2511-5020.rda")
QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]

# remember to make sure we work with values centered around 1 (the median is used to standardize)
# First compute the total energy for one flare in one viterbi sequence, 
# then make it an algorithm so that I do it for all flares in the sequences
# then further extend the algorithm so that it's and stored in a list for all sequences

data[is.na(data[,2]),2] <- median(data[,2], na.rm = TRUE)
data[,2] <- data[,2]/median(data[,2], na.rm = TRUE)

test <- data[2511:5020,]
aux <- Viterbi_raw[1,]>1



flare_int <- function(dat, vit){
  aux <- vit > 1
  dur <- 0
  flares_e <- NULL
  flare <- 0
  
  for(i in 1:length(aux)){
    
    if(aux[i]){
      flare = flare + (dat[i+1,1] - dat[i,1])/2*(dat[i+1,2] + dat[i,2])
      dur <- 1
    } else{
      if(dur){
        flares_e <- append(flares_e, flare)
        flare <- 0
        dur <- 0
      }
    }
    
  }
  
  return(flares_e)
}



flares_e1 <- NULL
for(i in 1:nrow(Viterbi_raw)){
  flares_e1 <- append(flares_e1, flare_int(data[2611:2711,], Viterbi_raw[i,100:200]))
}

e_dist1 <- data.frame("relative energy" = flares_e1)
ggplot(e_dist1, aes(x=relative.energy)) + 
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 35)


flares_e2 <- NULL
for(i in 1:nrow(Viterbi_raw)){
  flares_e2 <- append(flares_e2, flare_int(data[2811:2911,], Viterbi_raw[i,300:400]))
}

e_dist2 <- data.frame("relative energy" = flares_e2)
ggplot(e_dist2, aes(x=relative.energy)) + 
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 35)


flares_e3 <- NULL
for(i in 1:nrow(Viterbi_raw)){
  flares_e3 <- append(flares_e3, flare_int(data[4211:4411,], Viterbi_raw[i,1700:1900]))
}

e_dist3 <- data.frame("relative energy" = flares_e3)
ggplot(e_dist3, aes(x=relative.energy)) + 
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 35)


flares_e <- NULL
ind <- list(2:2510, 2511:5020, 5021:7530, 7531:9524, 12551:15060, 17571:20076)
ind_char <- list("2-2510", "2511-5020", "5021-7530", "7531-9524", "12551-15060", "17571-20076")
for(j in 1:length(ind)){
  if(j==1){
    load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/samples_2-2510.rda")
  } else{ if(j==3){
    load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/samples_5021-7530.rda")
  } else{
    load(file=paste0("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_",ind_char[[j]],".rda"))
    QFD_samples <- as.data.frame(fitQFD)
  }
  }
  N <- range(ind[[j]])[2]-range(ind[[j]])[1]+1
  Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]
  for(i in 1:nrow(Viterbi_raw)){
    flares_e <- append(flares_e, flare_int(data[ind[[j]],], Viterbi_raw[i,]))
  }
}


flares_e_maj <- NULL
ind <- list(2:2510, 2511:5020, 5021:7530, 7531:9524, 12551:15060, 17571:20076)
ind_char <- list("2-2510", "2511-5020", "5021-7530", "7531-9524", "12551-15060", "17571-20076")
for(j in 1:length(ind)){
  if(j==1){
    load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/samples_2-2510.rda")
  } else{ if(j==3){
    load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/samples_5021-7530.rda")
  } else{
    load(file=paste0("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_",ind_char[[j]],".rda"))
    QFD_samples <- as.data.frame(fitQFD)
  }
  }
  N <- range(ind[[j]])[2]-range(ind[[j]])[1]+1
  Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]
  flare_vit_max <- apply(Viterbi_raw, 2, majority)
  flares_e_maj <- append(flares_e_maj, flare_int(data[ind[[j]],], flare_vit_max))
  
}

e_dist <- data.frame("relative energy" = c(flares_e, flares_e_maj),
                     "Distribution" = c(rep("Viterbi Posterior Samples", length(flares_e)),
                                        rep("Viterbi Majority", length(flares_e_maj))))

pdf("/Users/arturoesquivel/Downloads/energy_distributions.pdf", width = 10, height = 6)
par(mfrow = c(2,1))

ggplot(e_dist1, aes(x=relative.energy)) + 
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 35) + 
  ggtitle("Flare 1 Energy Distribution")

ggplot(e_dist2, aes(x=relative.energy)) + 
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 35) + 
  ggtitle("Flare 2 Energy Distribution")

ggplot(e_dist3, aes(x=relative.energy)) + 
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 35) + 
  ggtitle("Flare 3 Energy Distribution")

ggplot(e_dist, aes(x=relative.energy, fill=Distribution)) + 
  geom_histogram(aes(y = after_stat(density)), alpha = 0.35, position = "identity", bins = 35) + 
  ggtitle("Energy Distribution (All Flares)")

# ggplot(e_dist, aes(x=relative.energy, fill=Distribution)) + 
#   geom_density(alpha = 0.35)

dev.off()









detrended_data <- data.frame(TIME = rawdata$TIME, Brightness = rawdata$PDCSAP_FLUX-summQFD[[1]][1:N+2*N+21, 1])
highlight_range <- data[2511:5020,1] >= 1328.918 & data[2511:5020,1] <= 1329
flare_det <- data[2511:5020, ][highlight_range,]
flare_vit <- Viterbi_raw[, highlight_range]
flare_vit_max <- apply(flare_vit,2,majority)

pdf("/Users/arturoesquivel/Downloads/flares2.pdf", width = 10, height = 6)
par(mfrow = c(2,1))
plot(flare_det, ylab = "Brightness", xlab = "TIME", main = "Flare 1 Decoding Across Posterior Samples")
points(flare_det[which(flare_vit_max==2),], col = "red",lwd=2.0)
points(flare_det[which(flare_vit_max==3),], col = "blue",lwd=2.0)


flare_vit_probs <- matrix(0, 3, ncol(flare_vit))
for(i in 1:ncol(flare_vit)){
  flare_vit_probs[1,i] <- sum(ifelse(flare_vit[,i] == 1, 1, 0))/nrow(flare_vit)
  flare_vit_probs[2,i] <- sum(ifelse(flare_vit[,i] == 2, 1, 0))/nrow(flare_vit)
  flare_vit_probs[3,i] <- sum(ifelse(flare_vit[,i] == 3, 1, 0))/nrow(flare_vit)
}

barplot(flare_vit_probs, 
        col=c("black","red","blue") , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="State Probabilities")


highlight_range2 <- data[2511:5020,1] >= 1329.23 & data[2511:5020,1] <= 1329.29
flare_det2 <- data[2511:5020, ][highlight_range2,]
flare_vit2 <- Viterbi_raw[, highlight_range2]
flare_vit_max2 <- apply(flare_vit2,2,majority)
par(mfrow = c(2,1))
plot(flare_det2, ylab = "Brightness", xlab = "Time", main = "Flare 2 Decoding Across Posterior Samples")
points(flare_det2[which(flare_vit_max2==2),], col = "red",lwd=2.0)
points(flare_det2[which(flare_vit_max2==3),], col = "blue",lwd=2.0)

flare_vit_probs2 <- matrix(0, 3, ncol(flare_vit2))
for(i in 1:ncol(flare_vit2)){
  flare_vit_probs2[1,i] <- sum(ifelse(flare_vit2[,i] == 1, 1, 0))/nrow(flare_vit2)
  flare_vit_probs2[2,i] <- sum(ifelse(flare_vit2[,i] == 2, 1, 0))/nrow(flare_vit2)
  flare_vit_probs2[3,i] <- sum(ifelse(flare_vit2[,i] == 3, 1, 0))/nrow(flare_vit2)
}

barplot(flare_vit_probs2, 
        col=c("black","red","blue") ,
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="State Probabilities")


highlight_range3 <- data[2511:5020,1] >= 1331.24 & data[2511:5020,1] <= 1331.37
flare_det3 <- data[2511:5020,][highlight_range3,]
flare_vit3 <- Viterbi_raw[, highlight_range3]
flare_vit_max3 <- apply(flare_vit3,2,majority)
par(mfrow = c(2,1))
plot(flare_det3, ylab = "Brightness", xlab = "Time", main = "Flare 3 Decoding Across Posterior Samples")
points(flare_det3[which(flare_vit_max3==2),], col = "red",lwd=2.0)
points(flare_det3[which(flare_vit_max3==3),], col = "blue",lwd=2.0)

flare_vit_probs3 <- matrix(0, 3, ncol(flare_vit3))
for(i in 1:ncol(flare_vit3)){
  flare_vit_probs3[1,i] <- sum(ifelse(flare_vit3[,i] == 1, 1, 0))/nrow(flare_vit3)
  flare_vit_probs3[2,i] <- sum(ifelse(flare_vit3[,i] == 2, 1, 0))/nrow(flare_vit3)
  flare_vit_probs3[3,i] <- sum(ifelse(flare_vit3[,i] == 3, 1, 0))/nrow(flare_vit3)
}

barplot(flare_vit_probs3, 
        col=c("black","red","blue") , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="State Probabilities")
dev.off()


# STart running, compute canada 
# tidy up plots 


load(file="/Users/arturoesquivel/Downloads/tess2018206045859-s0001-0000000089257479-0120-s_lc/model3.rda")
data <- read.csv("/Users/arturoesquivel/Downloads/tess2018206045859-s0001-0000000089257479-0120-s_lc/3.csv")[ ,c("TIME","PDCSAP_FLUX")]

# Tune hyperparameters so that th oscilation can be acounted for 
# Re-run second star with new hyperparameters 

data[,2] <- data[,2] - mean(data[,2], na.rm = TRUE)
observed <- (!is.na(data[,2])) * 1
data[is.na(data[,2]),2] <- 0
N <- nrow(data)

summQFD <- summary(fitQFD)

## visualize
QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]

Viterbi_max <- apply(Viterbi_raw,2,majority)

plot(data)
lines(data[,1], summQFD[[1]][1:N+2*N+21, 1], col = "#d400ff",lwd=3.0)
points(data[which(Viterbi_max==2)+1,], col = "red",lwd=3.0)
points(data[which(Viterbi_max==3)+1,], col = "blue",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
       lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
       cex = .75)





flares_e <- NULL
ind_char <- list("1", "2", "3", "5", "6", "8")
for(j in 1:length(ind_char)){

  load(file=paste0("/Users/arturoesquivel/Downloads/tess2018206045859-s0001-0000000089257479-0120-s_lc/model",ind_char[[j]],".rda"))
  data <- read.csv(paste0("/Users/arturoesquivel/Downloads/tess2018206045859-s0001-0000000089257479-0120-s_lc/",ind_char[[j]],".csv"))[ ,c("TIME","PDCSAP_FLUX")]
  data[is.na(data[,2]),2] <- median(data[,2], na.rm = TRUE)
  data[,2] <- data[,2]/median(data[,2], na.rm = TRUE)
  
  QFD_samples <- as.data.frame(fitQFD)
  N <- nrow(data) #N <- range(ind[[j]])[2]-range(ind[[j]])[1]+1
  Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]
  for(i in 1:nrow(Viterbi_raw)){
    flares_e <- append(flares_e, flare_int(data, Viterbi_raw[i,]))
  }
}


flares_e_maj <- NULL
states_time <- NULL
states_from <- NULL
states_to <- NULL
ind_char <- list("1", "2", "3", "5", "6", "8")
for(j in 1:length(ind_char)){
  load(file=paste0("/Users/arturoesquivel/Downloads/tess2018206045859-s0001-0000000089257479-0120-s_lc/model",ind_char[[j]],".rda"))
  data <- read.csv(paste0("/Users/arturoesquivel/Downloads/tess2018206045859-s0001-0000000089257479-0120-s_lc/",ind_char[[j]],".csv"))[ ,c("TIME","PDCSAP_FLUX")]
  data[is.na(data[,2]),2] <- median(data[,2], na.rm = TRUE)
  data[,2] <- data[,2]/median(data[,2], na.rm = TRUE)
  
  QFD_samples <- as.data.frame(fitQFD)
  N <- nrow(data) #N <- range(ind[[j]])[2]-range(ind[[j]])[1]+1
  Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]
  flare_vit_max <- apply(Viterbi_raw, 2, majority)
  flares_e_maj <- append(flares_e_maj, flare_int(data, flare_vit_max))
  
  for(i in 2:length(flare_vit_max)){
    if(flare_vit_max[i] != flare_vit_max[i-1]){
     states_time <- append(states_time, data[i,1])
     states_from <- append(states_from, flare_vit_max[i-1])
     states_to <- append(states_to, flare_vit_max[i])
    }
  }
  
}

e_dist <- data.frame("relative energy" = c(flares_e, flares_e_maj),
                     "Distribution" = c(rep("Viterbi Posterior Samples", length(flares_e)),
                                        rep("Viterbi Majority", length(flares_e_maj))))

pdf("/Users/arturoesquivel/Downloads/energy_distributions2.pdf", width = 10, height = 6)
par(mfrow = c(2,1))

e_dist2 <- e_dist[e_dist[,1] != max(e_dist[,1]), ]
ggplot(e_dist3, aes(x=relative.energy, fill=Distribution)) + 
  geom_histogram(aes(y = after_stat(density)), alpha = 0.35, position = "identity", bins = 35) + 
  ggtitle("Energy Distribution (All Flares)")

# ggplot(e_dist, aes(x=relative.energy, fill=Distribution)) + 
#   geom_density(alpha = 0.35)

dev.off()

# Look at the rhats for the problematic parameters ## all rhats less than 1.02, n_effs above 200
# Look into the one drive for the energy recovery plots material 
# automatize, get the rest of the quantities of interest
# Maybe cover relative energies to seconds

pl <- data.frame("Time" = data[2:N, 1],
                 "PDCSAP_FLUX" = data[2:N, 2],
                 "State" = as.factor(flare_vit_max),
                 "Trend" = summQFD[[1]][1:N+2*N+21, 1][2:N]
                 )

plot(data)
lines(data[,1], summQFD[[1]][1:N+2*N+21, 1], col = "#d400ff",lwd=3.0)
points(data[which(Viterbi_max==2)+1,], col = "red",lwd=3.0)
points(data[which(Viterbi_max==3)+1,], col = "blue",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
       lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
       cex = .75)

ggplot(pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = State, size = State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  
# annotate(geom = "rect",
#          xmin = 0, xmax = 1250, ymin = -Inf, ymax = Inf,
#          alpha = 1/6, linewidth = 0) +
#   geom_line(aes(color = chain),
#             linewidth = .15) + 
  theme_bw(base_size = 16) +
  labs(color = NULL, size = NULL) +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.08, .85)) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1))))
  #scale_linetype_manual(name = "Distribution", values=c("dashed", "solid")) +
    # scale_color_manual(values = c("No Cav" = "#003f5c", "Mild" = "#bc5090", "Severe" = "#ffa600", "Mixture" = "black"), 
    #                    name = 'State',
    #                    labels = c("No Cav" = "No Cav", "Mixture", "Severe", "Mixture")
    # ) +
    # labs(x = "F-V", y = "Density", title = "B: P-splines")  +
    # theme_bw(base_size = 13) +
    # theme(plot.title = element_text(face = "bold"))


# Light curves have no periodic pattern
# tess2018206045859-s0001-0000000176935123-0120-s_lc, running on screen 3
# tess2018206045859-s0001-0000000165096120-0120-s_lc 
# tess2018206045859-s0001-0000000139200489-0120-s_lc # only model8 missing for this one, running on screen 1

## tess2018206045859-s0001-0000000197736350-0120-s_lc is running on screen 2


N <- 1490
sim_pl <- data.frame("Time" = injected$TIME[2:N],
                     "PDCSAP_FLUX" = injected$PDCSAP_FLUX[2:N],
                     "GT" = as.factor(t(gtstate_df))[2:N],
                     "QFD_Trend" = summary(fitQFD)[[1]][1:N+2*N+21, 1][2:N],
                     "Sigma_Trend" = summcelerite[[1]][1:N + (N+23),1][2:N],
                     "QFD_State" = as.factor(t(QFD_df)),
                     "Sigma_State" = as.factor(t(sigma_df))[2:N],
                     "Flux" = t(flare_df)[2:N]
)



pl_1 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = GT, size = GT )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay")) +
  scale_size_manual(values=c(1.5,2,2),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay")) + 
  scale_shape_manual(values=c(19,19,19),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay")) + 
  theme_bw(base_size = 22) +
  labs(color = NULL, size = NULL, title = "Ground Truth", y ="PDCSAP") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.2, .855),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0))))


pl_2 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  theme_bw(base_size = 22) +
  labs(color = NULL, size = NULL, title = "Proposed HMM", y ="PDCSAP") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.26, .855),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1))))

pl_3 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  theme_bw(base_size = 22) +
  labs(color = NULL, size = NULL, title = "Sigma-clipping", y ="PDCSAP") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.20, .855),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0, 1))))


pl_4 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = Flux)) + 
  theme_bw(base_size = 22) +
  labs(color = NULL, size = NULL, title = "Flare Channel") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.08, .85))



pdf("/Users/arturoesquivel/Downloads/QFD_example_031381302.pdf", width = 10, height = 11)
ggarrange(pl_1, pl_2, pl_3, pl_4, ncol=1, nrow=4, common.legend = FALSE)#, legend="right")
dev.off()

## Periodic Stars:
### tess2018206045859-s0001-0000000234526939-0120-s_lc
### tess2018206045859-s0001-0000000238176832-0120-s_lc
### tess2018206045859-s0001-0000000273369281-0120-s_lc


#### Do we want to show them centered around 0 or 1? 

pl2_1 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  geom_rect(aes(xmin = 1326.15, xmax = 1326.45, ymin = -25, ymax = 50),
            fill = "transparent", color = "#56B4E9", size = 1.5) + 
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "Proposed HMM", y ="PDCSAP") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.38, .85),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1))))


pl2_2 <- sim_pl |> filter(Time >= 1326.15 & Time <= 1326.45) |>
ggplot(aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "") +
  theme(plot.title = element_text(face = "bold"),
        panel.background=element_rect(colour = "#56B4E9", size = 3),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        axis.title.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1))))


pl2_3 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  geom_rect(aes(xmin = 1326.15, xmax = 1326.45, ymin = -25, ymax = 50),
            fill = "transparent", color = "#56B4E9", size = 1.5) + 
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "Sigma-clipping", y ="PDCSAP") +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.27, .85),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0, 1))))


pl2_4 <- sim_pl |> filter(Time >= 1326.15 & Time <= 1326.45) |>
  ggplot(aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "") +
  theme(plot.title = element_text(face = "bold"),
        panel.background=element_rect(colour = "#56B4E9", size = 3),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        axis.title.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0, 1))))


pdf("/Users/arturoesquivel/Downloads/sim_det_compare.pdf", width = 8, height = 6)
ggarrange(pl2_1, pl2_2, pl2_3, pl2_4, ncol=2, nrow=2, common.legend = FALSE, widths = c(1, 0.5))
dev.off()


### reproduce box plots, data already there
### probability plots 
### window simulated and real data plot 

N <- 2001
zoom_pl <- data.frame("Time" = rawdata$TIME[2:N],
                     "PDCSAP_FLUX" = rawdata$PDCSAP_FLUX[2:N],
                     "QFD_Trend" = summary(fitQFD)[[1]][1:N+2*N+21, 1][2:N],
                     "Sigma_Trend" = summcelerite[[1]][1:N + (N+23),1][2:N],
                     "QFD_State" = as.factor(Viterbi_max),
                     "Sigma_State" = as.factor(flares3sigma + 1)[2:N]
)

pl3_1 <- ggplot(zoom_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  geom_rect(aes(xmin = 1334.9, xmax = 1335.25, ymin = -30, ymax = 60),
            fill = "transparent", color = "#56B4E9", size = 1.5) + 
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "Proposed HMM", y ="PDCSAP") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.38, .85),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1))))


pl3_2 <- zoom_pl |> filter(Time >= 1334.9 & Time <= 1335.25) |>
  ggplot(aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "") +
  theme(plot.title = element_text(face = "bold"),
        panel.background=element_rect(colour = "#56B4E9", size = 3),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        axis.title.y = element_blank()) +
  ylim(-30, 60) + 
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1))))


pl3_3 <- ggplot(zoom_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  geom_rect(aes(xmin = 1334.9, xmax = 1335.25, ymin = -30, ymax = 60),
            fill = "transparent", color = "#56B4E9", size = 1.5) + 
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "Sigma-clipping", y ="PDCSAP") +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.27, .85),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0, 1))))


pl3_4 <- zoom_pl |> filter(Time >= 1334.9 & Time <= 1335.25) |>
  ggplot(aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "") +
  theme(plot.title = element_text(face = "bold"),
        panel.background=element_rect(colour = "#56B4E9", size = 3),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        axis.title.y = element_blank()) +
  ylim(-30, 60) + 
  guides(color = guide_legend(override.aes = list(linetype = c(0,0, 1))))


pdf("/Users/arturoesquivel/Downloads/det_compare2.pdf", width = 8, height = 6)
ggarrange(pl3_1, pl3_2, pl3_3, pl3_4, ncol=2, nrow=2, common.legend = FALSE, widths = c(1, 0.5))
dev.off()

# for 5021-7530 chuk
flare_vit <- Viterbi_raw[, 2250:2300]
flare_vit_probs <- matrix(0, 3, ncol(flare_vit))
for(i in 1:ncol(flare_vit)){
     flare_vit_probs[1,i] <- sum(ifelse(flare_vit[,i] == 1, 1, 0))/nrow(flare_vit)
     flare_vit_probs[2,i] <- sum(ifelse(flare_vit[,i] == 2, 1, 0))/nrow(flare_vit)
     flare_vit_probs[3,i] <- sum(ifelse(flare_vit[,i] == 3, 1, 0))/nrow(flare_vit)
   }
flare_vit_probs <- as.data.frame(flare_vit_probs)
flare_vit_probs$State <- c(1,2,3)
probs_data <- flare_vit_probs |> pivot_longer(cols = V1:V121, names_to = "Time", values_to = "p")
probs_data$State <- as.factor(probs_data$State)
probs_data$Time <- rep(seq(1,length(probs_data$Time)/3),3)

(probs_pl <- ggplot(data = probs_data, aes(x = Time, y = p, fill = State)) +
  geom_bar(stat="identity", width = .85) + 
  scale_fill_manual(values = c("#000000", "#D55E00", "#E69F00"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay")) +
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, fill = NULL, title = "Viterbi Decoding Posterior Proportions") +
  theme(plot.title = element_text(face = "bold"),
        panel.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.direction = "horizontal",
        axis.title.x=element_blank(),
        axis.title.y = element_blank()))

full_data$state <- as.factor(full_data$state)

pdf("/Users/arturoesquivel/Downloads/031381302-det.pdf", width = 9, height = 5)
ggplot(full_data[2:9524, ], aes(x = TIME)) +
  geom_point(aes( y = PDCSAP_FLUX, col = state, size = state )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  theme_bw(base_size = 22) +
  labs(color = NULL, size = NULL, title = "Proposed HMM Fit", x = "Time", y ="PDCSAP") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.30, .85),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1))))
dev.off()


(flux_pl <- ggplot(full_data[2621:2670, ], aes(x = TIME)) +
  geom_point(aes( y = PDCSAP_FLUX, col = state, size = state )) + #, shape = State)) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
  theme_bw(base_size = 15) +
  labs(color = NULL, size = NULL, title = "Light Curve") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.29, .85),
        legend.direction = "horizontal",
        axis.title.y=element_blank(),
        legend.background = element_rect(fill='transparent')) +
  guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1)))))


pdf("/Users/arturoesquivel/Downloads/full_star_probs.pdf", width = 9, height = 6)
ggarrange(probs_pl, flux_pl, ncol=1, nrow=2, common.legend = TRUE)
dev.off()


(probs_pl2 <- ggplot(data = probs_data, aes(x = Time, y = p, fill = State)) +
    geom_bar(stat="identity", width = .85) + 
    scale_fill_manual(values = c("#000000", "#D55E00", "#E69F00"),
                      labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay")) +
    theme_bw(base_size = 15) +
    labs(color = NULL, size = NULL, fill = NULL, title = "Viterbi Decoding Posterior Proportions") +
    theme(plot.title = element_text(face = "bold"),
          panel.background=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.title.y = element_blank()))

(flux_pl2 <- ggplot(full_data[7271:7321, ], aes(x = TIME)) +
    geom_point(aes( y = PDCSAP_FLUX, col = state, size = state )) + #, shape = State)) +
    scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                       labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
    scale_size_manual(values=c(1.5,2,2,1.5),
                      labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
    scale_shape_manual(values=c(19,19,19,NA),
                       labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
    geom_line(aes(y = Trend, col = "Trend", size = "Trend")) +#, shape = "Trend")) +
    theme_bw(base_size = 15) +
    labs(color = NULL, size = NULL, title = "Light Curve") +
    theme(plot.title = element_text(face = "bold"),
          panel.background = element_blank(),
          panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
          panel.grid.minor = element_blank(),
          legend.position = c(0.29, .85),
          legend.direction = "horizontal",
          axis.title.y=element_blank(),
          legend.background = element_rect(fill='transparent')) +
    guides(color = guide_legend(override.aes = list(linetype = c(0,0,0, 1)))))

pdf("/Users/arturoesquivel/Downloads/full_star_probs2.pdf", width = 9, height = 6)
ggarrange(probs_pl2, flux_pl2, ncol=1, nrow=2, common.legend = TRUE)
dev.off()

fit_states <- as.numeric(full_data[2:9524, ]$state) > 1
location_flares <- which(fit_states)

flare_count <- (length(location_flares)!=0)+sum((location_flares[-1]-location_flares[-length(location_flares)])!=1)

starting <- location_flares[-1]
ending <- location_flares[-length(location_flares)]
tt <- which((location_flares[-1]-location_flares[-length(location_flares)])!=1)

durations <- c(2, 5, 10, 3, 25, 10, 4, 24, 8, 12)
mean(durations)


for(j in 1:length(ind)){
  if(j==1){
    load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/samples_2-2510.rda")
  } else{ if(j==3){
    load(file="/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/samples_5021-7530.rda")
  } else{
    load(file=paste0("/Users/arturoesquivel/Library/CloudStorage/OneDrive-Personal/Documentos/GitHub/flares/Results/tess2018206045859-s0001-0000000031381302-0120-s_lc/model_",ind_char[[j]],".rda"))
    QFD_samples <- as.data.frame(fitQFD)
  }
  }
  firing_rate <- append(firing_rate, QFD_samples["theta_quiet[2]"]$`theta_quiet[2]`)
  
}


