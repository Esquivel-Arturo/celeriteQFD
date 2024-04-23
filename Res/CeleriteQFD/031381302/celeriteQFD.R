library(rstan)
library(dplyr)
library(ggpubr)
library(ggplot2)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)
source("./celeriteQFD/R/misc.R") # Helper functions

# run QFD
rawdata <- read.csv("./celeriteQFD/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[6000:8000,c("TIME","PDCSAP_FLUX")]

# Center flux and mark missing observations
rawdata[,2] <- rawdata[,2] - mean(rawdata[,2], na.rm = TRUE)
observed <- (!is.na(rawdata[,2])) * 1
rawdata[is.na(rawdata[,2]),2] <- 0
N <- nrow(rawdata)
plot(rawdata)
tt <- rawdata[,1]

# Fit celeriteQFD
QFD_data <- list(N=N, t = rawdata[,1],
                y = rawdata[,2],
                observed = observed,
                sigma_prior = c(-8,8),
                Q0_prior = c(0,8), # This values seem to work for most stars
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

modelQFD <- stan_model(file = './celeriteQFD/Stan/Morphology/QFD/CeleriteQFDexN-missing-handling.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path('./celeriteQFD', 
                             'celerite2/celerite2.hpp'), '"\n'))

fitQFD <- sampling(modelQFD, 
                   data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), 
                   iter = 2000, 
                   init_r = 15, 
                   chains = 2)
# Model summary
summQFD <- summary(fitQFD)

# Inspect trace plots
rstan::traceplot(fitQFD, pars = c("theta_quiet", "theta_firing", "theta_decay"))
rstan::traceplot(fitQFD, pars = c("rho1", "w2", "tau2"))
rstan::traceplot(fitQFD, pars = c("trend[1]","trend[50]", "trend[600]", "trend[1400]", "trend[1990]"))
rstan::traceplot(fitQFD, pars = c("lsigma", "lperiod", "mu_quiet", "sigma", "rate_decay"))
rstan::traceplot(fitQFD, pars = c("eta[11]","eta[150]", "eta[400]", "eta[1403]", "eta[1992]"))

## visualize model fit 
QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]

# Select the most frequent state
Viterbi_max <- apply(Viterbi_raw,2,majority)

# Plot data 
data <- data.frame("TIME" = rawdata$TIME[2:N],
                   "PDCSAP_FLUX" = rawdata$PDCSAP_FLUX[2:N],
                   "Trend" = summQFD[[1]][2:N+2*N+21, 1],
                   "state" = as.factor(Viterbi_max))

ggplot(data, aes(x = TIME)) +
  geom_point(aes( y = PDCSAP_FLUX, col = state, size = state )) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = Trend, col = "Trend", size = "Trend")) +
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


# Fit celerite alone

modelcelerite <- stan_model(file = './celeriteQFD/Stan/Prototypes/Celerite/celerite-missing-handling.stan', 
            model_name = "celerit2", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path('./celeriteQFD/', 
                             'celerite2/celerite2.hpp'), '"\n'))

celeritedata <- QFD_data
celeritedata$err_prior <- c(0.01, 0.01)

fitcelerite <- sampling(modelcelerite, 
                        data = celeritedata,control = list(adapt_delta = 0.99, max_treedepth=15), 
                        iter = 2000,
                        init_r = 2, 
                        chains = 2)

summcelerite <- summary(fitcelerite)
celerite_trend <- summcelerite[[1]][1:N + (N+23),1]
# Detrendng
residual <- rawdata[,2] - celerite_trend
# 1-3sigma
flares3sigma <- residual >= (mean(residual) + 3 * sd(residual))


##### Fig. 2 in the paper #####
# After using the code above to model the entries 2511:5020 of the data (full_data)
flare_vit <- Viterbi_raw[, 1756:1886]
flare_vit_probs <- matrix(0, 3, ncol(flare_vit))
# State proportions for each observation
for(i in 1:ncol(flare_vit)){
  flare_vit_probs[1,i] <- sum(ifelse(flare_vit[,i] == 1, 1, 0))/nrow(flare_vit)
  flare_vit_probs[2,i] <- sum(ifelse(flare_vit[,i] == 2, 1, 0))/nrow(flare_vit)
  flare_vit_probs[3,i] <- sum(ifelse(flare_vit[,i] == 3, 1, 0))/nrow(flare_vit)
}
flare_vit_probs <- as.data.frame(flare_vit_probs)
flare_vit_probs$State <- c(1,2,3)
probs_data <- flare_vit_probs |> pivot_longer(cols = V1:V131, names_to = "Time", values_to = "p")
probs_data$State <- as.factor(probs_data$State)
probs_data$Time <- rep(seq(1,length(probs_data$Time)/3),3)

# Bar plot
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

# Light curve and color coded observations
(flux_pl <- ggplot(full_data[4267:4397, ], aes(x = TIME)) +
    geom_point(aes( y = PDCSAP_FLUX, col = state, size = state )) + 
    scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                       labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
    scale_size_manual(values=c(1.5,2,2,1.5),
                      labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
    scale_shape_manual(values=c(19,19,19,NA),
                       labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
    geom_line(aes(y = Trend, col = "Trend", size = "Trend")) +
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

# Output plots
pdf("./celeriteQFD/Res/CeleriteQFD/031381302/full_star_probs.pdf", width = 9, height = 6)
ggarrange(probs_pl, flux_pl, ncol=1, nrow=2, common.legend = TRUE)
dev.off()


##### Fig. 3 in the paper #####
# After using the code above to model the entries 5021:7530 of the data (full_data)
flare_vit <- Viterbi_raw[, 2250:2300]
flare_vit_probs <- matrix(0, 3, ncol(flare_vit))
# State proportions for each observation
for(i in 1:ncol(flare_vit)){
  flare_vit_probs[1,i] <- sum(ifelse(flare_vit[,i] == 1, 1, 0))/nrow(flare_vit)
  flare_vit_probs[2,i] <- sum(ifelse(flare_vit[,i] == 2, 1, 0))/nrow(flare_vit)
  flare_vit_probs[3,i] <- sum(ifelse(flare_vit[,i] == 3, 1, 0))/nrow(flare_vit)
}
flare_vit_probs <- as.data.frame(flare_vit_probs)
flare_vit_probs$State <- c(1,2,3)
probs_data <- flare_vit_probs |> pivot_longer(cols = V1:V51, names_to = "Time", values_to = "p")
probs_data$State <- as.factor(probs_data$State)
probs_data$Time <- rep(seq(1,length(probs_data$Time)/3),3)
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
    geom_point(aes( y = PDCSAP_FLUX, col = state, size = state )) + 
    scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                       labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
    scale_size_manual(values=c(1.5,2,2,1.5),
                      labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
    scale_shape_manual(values=c(19,19,19,NA),
                       labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
    geom_line(aes(y = Trend, col = "Trend", size = "Trend")) +
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

pdf("./celeriteQFD/Res/CeleriteQFD/031381302/full_star_probs2.pdf", width = 9, height = 6)
ggarrange(probs_pl2, flux_pl2, ncol=1, nrow=2, common.legend = TRUE)
dev.off()



### Fig. 4 in the paper ### 

# After running 1 time the code in tests/injection_recovery/flare_injection_CHTC to simulate and inject flares
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

# Ground Truth
pl_1 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = GT, size = GT )) + 
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

# celeriteQFD fit
pl_2 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +
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

# celerite/sigma-clipping fit
pl_3 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) +
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +
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

# Injected channel
pl_4 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = Flux)) + 
  theme_bw(base_size = 22) +
  labs(color = NULL, size = NULL, title = "Flare Channel") +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "#999999",size = 0.5, linetype = 2),
        panel.grid.minor = element_blank(),
        legend.position = c(0.08, .85))



pdf("./celeriteQFD/Res/CeleriteQFD/031381302/QFD_example_031381302.pdf", width = 10, height = 11)
ggarrange(pl_1, pl_2, pl_3, pl_4, ncol=1, nrow=4, common.legend = FALSE)
dev.off()



## Fig 6 ##

pl2_1 <- ggplot(sim_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) + 
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +
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
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +
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
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + 
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +
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
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + 
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +
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


pdf("./celeriteQFD/Res/CeleriteQFD/031381302/sim_det_compare.pdf", width = 8, height = 6)
ggarrange(pl2_1, pl2_2, pl2_3, pl2_4, ncol=2, nrow=2, common.legend = FALSE, widths = c(1, 0.5))
dev.off()


pl <- data.frame("Time" = data[2:N, 1],
                 "PDCSAP_FLUX" = data[2:N, 2],
                 "State" = as.factor(flare_vit_max),
                 "Trend" = summQFD[[1]][1:N+2*N+21, 1][2:N]
)


## Fig 8 ##
N <- 2001
zoom_pl <- data.frame("Time" = rawdata$TIME[2:N],
                     "PDCSAP_FLUX" = rawdata$PDCSAP_FLUX[2:N],
                     "QFD_Trend" = summary(fitQFD)[[1]][1:N+2*N+21, 1][2:N],
                     "Sigma_Trend" = summcelerite[[1]][1:N + (N+23),1][2:N],
                     "QFD_State" = as.factor(Viterbi_max),
                     "Sigma_State" = as.factor(flares3sigma + 1)[2:N]
)

pl3_1 <- ggplot(zoom_pl, aes(x = Time)) +
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +
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
  geom_point(aes( y = PDCSAP_FLUX, col = QFD_State, size = QFD_State )) +
  scale_color_manual(values = c("#000000", "#D55E00", "#E69F00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) +
  scale_size_manual(values=c(1.5,2,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  scale_shape_manual(values=c(19,19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Firing", "3" = "Decay", "Trend")) + 
  geom_line(aes(y = QFD_Trend, col = "Trend", size = "Trend")) +
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
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + 
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) +
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
  geom_point(aes( y = PDCSAP_FLUX, col = Sigma_State, size = Sigma_State )) + 
  scale_color_manual(values = c("#000000", "#D55E00", "#CC79A7"),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) +
  scale_size_manual(values=c(1.5,2,1.5),
                    labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  scale_shape_manual(values=c(19,19,NA),
                     labels = c("1" = "Quiet", "2" = "Flare", "Trend")) + 
  geom_line(aes(y = Sigma_Trend, col = "Trend", size = "Trend")) ++
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


pdf("./celeriteQFD/Res/CeleriteQFD/031381302/det_compare2.pdf", width = 8, height = 6)
ggarrange(pl3_1, pl3_2, pl3_3, pl3_4, ncol=2, nrow=2, common.legend = FALSE, widths = c(1, 0.5))
dev.off()

