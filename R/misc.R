majority <- function(x){
    as.numeric(names(which.max(table(x))))
}

error_decoding <- function(decoding, states){
    res <- data.frame(matrix(NA,1,7))
    colnames(res) <- c("n_detected","n_injected","TP","FP","FN","SEN","SPC")
    
    states <- states > 1
    decoding <- decoding > 1
    
    location_inj <- which(states)
    res$n_injected <- (length(location_inj)!=0)+sum((location_inj[-1]-location_inj[-length(location_inj)])!=1)
    location_decoding <- which(decoding)
    res$n_detected <- (length(location_decoding)!=0)+sum((location_decoding[-1]-location_decoding[-length(location_decoding)])!=1)
    
    location_TP <- which( states & decoding)
    res$TP <- count_success(location_TP, location_inj)
    location_FP <- which(!states & decoding)
    
    res$FP <- count_success(location_FP, location_decoding)
    res$FN <- res$n_injected-res$TP
    res$SEN <- res$TP/res$n_injected
    res$SPC <- res$TP/res$n_detected
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
    
    loss[1,9:15] <- error_decoding(QFD,gt[-1])
    loss[2,9:15] <- error_decoding(sigmarule, gt)
    loss$res_dir <- thedir
    
    return(loss)
    
}