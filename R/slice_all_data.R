all_csvs <- list.files("../../Data/TimeSeriesCSV",pattern = "*.csv", full.names = TRUE)
path_to_store <- "../../Data/TimeSeriesCSV_sliced/"
source("R/slicer.R")

the_log <- lapply(all_csvs, slicer, path_to_store = path_to_store, target_length = 2500, length_bounds = c(1500,3000), max_NA = 500)

the_log <- Reduce(f = c, x = the_log)
write.csv(data.frame(the_log), "../../Data/TimeSeriesCSV_sliced/slicer_log_overall.csv")

all_logs <- list.files("../../Data/TimeSeriesCSV_sliced",pattern = "*log.csv", full.names = TRUE) |>
    lapply(read.csv) |>
    Reduce(f = rbind)

write.csv(all_logs, "../../Data/TimeSeriesCSV_sliced/slicer_log_by_star.csv")

all_logs <- list.files("../../Data/TimeSeriesCSV_sliced",pattern = "*log.csv", full.names = TRUE) |>
    lapply(file.remove)
