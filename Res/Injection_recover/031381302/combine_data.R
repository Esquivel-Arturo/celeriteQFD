library(reshape2)
library(ggplot2)
error_decoding <- function(decoding, states){
  res <- data.frame(matrix(NA,1,10))
  colnames(res) <- c("n_detected","n_injected","TP","FP","FN","SEN","PPV", "fSEN", "fSPC", "fPPV")
  
  states <- states > 1
  decoding <- decoding > 1
  
  location_inj <- which(states)
  # Count different flare occurrences, without considering specific durations and associated timesteps 
  res$n_injected <- (length(location_inj)!=0)+sum((location_inj[-1]-location_inj[-length(location_inj)])!=1)
  # Count number of different flares detected
  location_decoding <- which(decoding)
  res$n_detected <- (length(location_decoding)!=0)+sum((location_decoding[-1]-location_decoding[-length(location_decoding)])!=1)
  
  location_TP <- which( states & decoding)
  # Number of true flares detected, for any time step detected within the duration of a flare
  res$TP <- count_success(location_TP, location_inj)
  location_FP <- which(!states & decoding)
  # Number of false flares detected, for any time step detected outside the duration of a flare
  res$FP <- count_success(location_FP, location_decoding)
  
  res$FN <- res$n_injected-res$TP
  res$SEN <- res$TP/res$n_injected
  res$PPV <- res$TP/(res$TP+res$FP) 
  
  # Sensitivity and PPV in terms of all time steps
  res$fSEN <- length(location_TP)/length(location_inj)
  
  location_TN <- which(!states & !decoding)
  location_ninj <- which(!states)
  res$fSPC <- length(location_TN)/length(location_ninj)
  
  res$fPPV <- length(location_TP)/(length(location_TP) + length(location_FP))
  
  return(res)
  
}


count_success <- function(location_both, location_base){
  n_true <- (length(location_base)!=0)+sum((location_base[-1]-location_base[-length(location_base)])!=1)
  if(n_true==1 & length(location_both)>0) return(1)
  starting <- location_base[-1]
  ending <- location_base[-length(location_base)]
  tt <- which((location_base[-1]-location_base[-length(location_base)])!=1)
  
  starting <- starting[tt]
  starting <- c(location_base[1], starting)
  ending <- ending[tt]
  ending <- c(ending, location_base[length(location_base)])
  
  success <- 0
  for(i in 1:n_true){
    temp <- (location_both >= starting[i]) & (location_both<=ending[i])
    if(sum(temp)>0){
      success <- success + 1
    }
    
  }
  
  return(success)
  
}

clean_CHTC_data <- function(thedir){
  loss <- read.csv(paste(thedir,"inj_rec_loss.csv",sep = "/"),row.names = 1)
  QFD <- read.csv(paste(thedir,"inj_rec_QFD.csv",sep = "/"),row.names = 1)
  QFD <- as.numeric(QFD)
  gt <- read.csv(paste(thedir,"inj_rec_gtstate.csv",sep = "/"),row.names = 1)
  gt <- as.numeric(gt)
  sigmarule <- read.csv(paste(thedir,"inj_rec_3sigma.csv",sep = "/"),row.names = 1)
  sigmarule <- as.numeric(sigmarule)
  
  loss$fSEN <- c(0,0)
  loss$fSPC <- c(0,0)
  loss$fPPV <- c(0,0)
  
  loss[1,9:18] <- error_decoding(QFD,gt[-1])
  loss[2,9:18] <- error_decoding(sigmarule, gt)
  loss$res_dir <- thedir
  
  return(loss)
  
}

all_res_dir <- list.dirs()[-c(1,2)]

all_res_df <- lapply(all_res_dir, clean_CHTC_data)

all_res_df1 <- Reduce(rbind, all_res_df)
# Store results
write.csv(all_res_df1,"test_flares_full.csv", row.names = F)


### Sensitivity/PPV plots (Fig. 5 on the paper)###

colnames(all_res_df1)[14:15] <- c("Sensitivity","PPV")
# per observation results (space used at the end to distinguish)
colnames(all_res_df1)[16:18] <- c("Sensitivity ", "Specificity", "PPV ")
plot_data <- melt(all_res_df1[,], measure.vars = colnames(all_res_df1)[9:18])


(bp1 <- ggplot(plot_data[plot_data$variable %in% c("Sensitivity","PPV"),], 
               aes(x = variable, y = value)) +
    geom_boxplot(position = position_dodge(1), size = 1, outlier.colour = NULL, outlier.alpha = 0.5, 
                 aes(colour = method)) + 
    ylab("") + 
    xlab("") + 
    theme_bw(base_size = 25) +
    labs(title = "Small Flares") +
    theme(plot.title = element_text(face = "bold"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),        
          legend.position = "bottom",
          legend.direction = "horizontal",
          panel.grid.major = element_blank()) +
    scale_color_manual(values=c("#56B4E9", "#0072B2"),
                       labels = c("1-3-sigma" = expression(paste("1-3", sigma)),
                                  "QFD" = "HMM"),
                       name = ""))
# Large flares plot
setwd("./celeriteQFD/Res/Injection_recover/031381302/large_flares/")
all_res_dir <- list.dirs()[-c(1,2)]

all_res_df <- lapply(all_res_dir, clean_CHTC_data)

all_res_df2 <- Reduce(rbind, all_res_df)
colnames(all_res_df2)[14:15] <- c("Sensitivity","PPV")
colnames(all_res_df2)[16:18] <- c("Sensitivity ", "Specificity", "PPV ")
plot_data2 <- melt(all_res_df2[,], measure.vars = colnames(all_res_df2)[9:18])

(bp2 <- ggplot(plot_data2[plot_data2$variable %in% c("Sensitivity","PPV"),], 
               aes(x = variable, y = value)) +
    geom_boxplot(position = position_dodge(1), size = 1, outlier.colour = NULL, outlier.alpha = 0.5, 
                 aes(colour = method)) + 
    ylab("") + 
    xlab("") + 
    theme_bw(base_size = 25) +
    labs(title = "Large Flares") +
    theme(plot.title = element_text(face = "bold"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),        
          legend.position = "bottom",
          legend.direction = "horizontal",
          panel.grid.major = element_blank()) +
    scale_color_manual(values=c("#56B4E9", "#0072B2"),
                       labels = c("1-3-sigma" = expression(paste("1-3", sigma)),
                                  "QFD" = "HMM"),
                       name = ""))

pdf("./box_plots.pdf", width = 8, height = 7)
ggarrange(bp1, bp2, ncol=1, nrow=2, common.legend = TRUE)
dev.off()

# For per observation results
(bpf1 <- ggplot(plot_data[plot_data$variable %in% c("Sensitivity ","PPV "),], 
                aes(x = variable, y = value)) +
    geom_boxplot(position = position_dodge(1), size = 1, outlier.colour = NULL, outlier.alpha = 0.5, 
                 aes(colour = method)) + 
    ylab("") + 
    xlab("") + 
    theme_bw(base_size = 25) +
    labs(title = "Small Flares") +
    theme(plot.title = element_text(face = "bold"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),        
          legend.position = "bottom",
          legend.direction = "horizontal",
          panel.grid.major = element_blank()) +
    scale_color_manual(values=c("#56B4E9", "#0072B2"),
                       labels = c("1-3-sigma" = expression(paste("1-3", sigma)),
                                  "QFD" = "HMM"),
                       name = ""))


(bpf2 <- ggplot(plot_data2[plot_data2$variable %in% c("Sensitivity ","PPV "),], 
                aes(x = variable, y = value)) +
    geom_boxplot(position = position_dodge(1), size = 1, outlier.colour = NULL, outlier.alpha = 0.5, 
                 aes(colour = method)) + 
    ylab("") + 
    xlab("") + 
    theme_bw(base_size = 25) +
    labs(title = "Large Flares") +
    theme(plot.title = element_text(face = "bold"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),        
          legend.position = "bottom",
          legend.direction = "horizontal",
          panel.grid.major = element_blank()) +
    scale_color_manual(values=c("#56B4E9", "#0072B2"),
                       labels = c("1-3-sigma" = expression(paste("1-3", sigma)),
                                  "QFD" = "HMM"),
                       name = ""))



