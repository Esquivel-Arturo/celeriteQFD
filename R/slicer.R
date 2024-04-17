# function to slice timeseries data smart

remove_start_na <- function(df){
    if(is.na(df$PDCSAP_FLUX[1])){
        return(remove_start_na(df[-1,]))
    }

    else{
        return(df)
    }
}

remove_ending_na <- function(df){
    if(is.na(df$PDCSAP_FLUX[length(df$PDCSAP_FLUX)])){
        return(remove_ending_na(df[-length(df$PDCSAP_FLUX),]))
    }

    else{
        return(df)
    }
}

slicer <- function(file_name, path_to_store, target_length = 2500,length_bounds = c(1500,3000), max_NA = 500, verbose = T){
    if(verbose){
        cat(file_name, "\n")
    }
    data <- tryCatch( read.csv(file_name) , error = function(e){
        if(verbose){
            cat("Error: ", e$message, "\n")
        }
        return(NULL)
    })
    if(is.null(data)){
        return(file_name)
    }
    n <- nrow(data) # number of rows
    
    # initial guess
    remaining <- n %% target_length
    if(remaining<length_bounds[1]){
        n_slices <- floor(n/target_length)
        # if add all remainings to each slice not causing problem and will move less than moving previous to the last
        if(ceiling(remaining/n_slices) <= (length_bounds[2]-target_length) ){
            # how many steps to add
            steps_add_to_slice_in_full <- ceiling(remaining/n_slices)
            slices_to_add_in_full <- floor(remaining/steps_add_to_slice_in_full)

            # after adding there are still some remainings
            still_remain <- remaining - (slices_to_add_in_full*steps_add_to_slice_in_full)
            
            length_slices <- rep(target_length,n_slices)
            length_slices[1:slices_to_add_in_full] <- target_length + steps_add_to_slice_in_full
            length_slices[slices_to_add_in_full+1] <- target_length + still_remain
            
        }
        ## if allocate remaining to everyone causing problem
        else{
            ## see if we can take something out from the previous slices
            if(ceiling((length_bounds[1]-remaining) / n_slices) < length_bounds[2]-target_length){
                length_slices <- rep(target_length,n_slices)
                length_slices <- length_slices - ceiling((length_bounds[1]-remaining) / n_slices)
                length_slices[n_slices + 1] <- remaining + ceiling((length_bounds[1]-remaining) / n_slices) * n_slices
            }
            else{
                length_slices <- rep(target_length,n_slices)
                length_slices[n_slices + 1] <- remaining
            }
                
        }
        idx_high <- cumsum(length_slices)
        idx_low <- idx_high - length_slices + 1
    }
    else{
        n_slices <- ceiling(n/target_length)
        idx_low <- (1:n_slices-1)*target_length+1
        idx_high <- (1:n_slices)*target_length
        idx_high[n_slices] <- n
    }

    logs <- data.frame(ori_file = file_name, idx_low = idx_low, idx_high = idx_high, output = NA, note = NA)
    for(i in 1:length(idx_low)){
        temp <- data[idx_low[i]:idx_high[i],]
        if(sum(is.na(temp$PDCSAP_FLUX))>max_NA){
            logs[i,]$note <- "too many NA"
        }
        else{
            temp <- remove_start_na(temp)
            temp <- remove_ending_na(temp)
            write.csv(temp, paste0(path_to_store, sub( ".csv","", basename(file_name)), "_", i, ".csv"), row.names = FALSE)
            logs[i,]$output <- paste0(path_to_store, sub( ".csv","", basename(file_name)), "_", i, ".csv")
        }
        write.csv(logs, paste0(path_to_store, sub( ".csv","", basename(file_name)), "_log.csv"))
        
    } 
    return(TRUE)

}
