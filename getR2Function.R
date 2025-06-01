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

df_r2 <- data.frame(matrix(ncol = 0, nrow = 1))
rownames(df_r2) <- c("R2")

if (length(group_data) > 1) {
  for (i in 3:6) {
    model <- lm(group_data[, i] ~ real_abundance, data = group_data)
    r2 <- summary(model)$r.squared
    df_r2[[colnames(group_data)[i]]] <- r2
  }
} else {
  r2 <- c(NA)
  return(NA)
}

write.csv(df_r2, paste("R2_", original_name), row.names = TRUE)
