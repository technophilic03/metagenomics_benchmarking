args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(tidyverse)})
real_data_file <- args[1]	# rtmp($REALDATAFILE)
sim_data_file <- args[2]	# stmp($SIMDATAFILE)
total_genomes <- args[3] # TOTALGENOMES
original_name <- args[4] # $ORIGINALSNAME


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

missing_in_real <- anti_join(sim, real, by = "Species")
missing_in_sim <- anti_join(real, sim, by = "Species")

df_roc <- data.frame(matrix(ncol = 0, nrow = 2))
rownames(df_roc) <- c("TPR", "FPR") # True positive rate and false positive rate

TP <- 0
FP <- 0
FN <- 0
FP <- FP + nrow(missing_in_real) # FP when a sequence found in the simulation isn't in the real data.
FN <- FN + nrow(missing_in_sim) # FN when sequences in real data but missing in simulation
TN <- as.numeric(total_genomes)  - nrow(sim) # TN = total genomes - total simulated entries

if (length(group_data) > 1) {
  for (i in 3:6) {
    for (j in 1:nrow(group_data)) {
      if (group_data[j, i] >= group_data$real_abundance[j] * 0.5) {
        TP <- TP + 1 # TP + 1 when sim_read >= real_reads/2
      } else if (group_data[j, i] == 0 &&
                 group_data$real_abundance[j] > 0) {
        FP <- FP + 1 # FP +1 when simulated reads = 0 but real reads >0
      }
    }
    
    get_TPR <- TP / (TP + FN)
    get_FPR <- FP / (FP + TN)
    df_roc[[colnames(group_data)[i]]] <- c(get_TPR, get_FPR)
    t_df_roc <- as.data.frame(t(df_roc)) # transposed
    
  }
  
} else{
  return(NA)
}



write.csv(t_df_roc, paste("ROC_", original_name, ".csv"), row.names = TRUE)
