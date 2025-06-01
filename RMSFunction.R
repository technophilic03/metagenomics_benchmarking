args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(tidyverse)})

real_data_file <- args[1]	# rtmp($REALDATAFILE)
sim_data_file <- args[2]	# stmp($SIMDATAFILE)
original_name <- args[3] # $ORIGINALSNAME

real_data <- read_table(real_data_file, col_names = FALSE, skip = 1)
sim_data <- read_table("stmp", col_names = FALSE, skip = 1)
sim_data_header <- read_table("stmp", col_names = FALSE, n_max = 1)

real <- real_data |>
  mutate(Species = paste(X1, X2), .before = 1) |>
  select(-X1, -X2) |>
  rename(real_abundance = X3)

header <- sim_data_header[, -(1:2)]

col_mapping <- setNames(
  paste0("X", 3 + seq_along(header)),
  header
)

sum_cols <- setNames(
  lapply(col_mapping, function(col) expr(sum(!!sym(col)))),
  paste0("sum_", names(col_mapping))
)


sim <- sim_data |>
  mutate(Species = paste(X1, X2), .before = 1) |>
  group_by(Species) |>
  summarise(!!!sum_cols)

group_data <- merge(real, sim, by = "Species")

df_rms <- data.frame(matrix(ncol = 0, nrow = 2))
rownames(df_rms) <- c("RRMSE", "AVGRE")

# rrmse = rmse / mean of observation
# observed: SIMDATA/true microbial abd
# predicted: REALDATA/ground truth/real_abundance

# excellent <10%< good 10<20%, fair 20<30, poor >30
get_rrmse <- function(observed, predicted) {
  squared_diff <- (observed - predicted) ^ 2
  mse <- mean(squared_diff)
  rmse <- sqrt(mse)
  mean_observed <- mean(observed)
  rrmse <- (rmse / mean_observed) * 100
  return(rrmse)
}

get_avgre <- function(observed, predicted) {
  if (0 %in% observed) {
    valid_indices <- observed != 0
    observed <- observed[valid_indices]
    predicted <- predicted[valid_indices]
  }
  avgre <- mean(abs(observed - predicted) / observed)
  return(avgre)
  
}

if (length(group_data) > 1) {
  for (i in 3:6) {
    rrmse <- get_rrmse(group_data[, i], group_data$real_abundance)
    avgre <- get_avgre(group_data[, i], group_data$real_abundance)
    df_rms[[colnames(group_data)[i]]] <- c(rrmse, avgre)
  }
} else {
  r2 <- c(NA)
  return(NA)
}


write.csv(df_rms, paste("RMS_", original_name, ".csv"), row.names = TRUE)
