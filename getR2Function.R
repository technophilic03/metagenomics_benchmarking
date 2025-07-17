args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(tidyverse)})
real_data_file <- args[1]	# rtmp($REALDATAFILE)
sim_data_file <- args[2]	# stmp($SIMDATAFILE)
original_name <- args[3] # $ORIGINALSNAME

real_data <- read.csv(real_data_file)
sim_data <- read.csv(sim_data_file)

real <- real_data %>%
  rename(Species = species, real_abundance = abundance)

sim <- sim_data %>%
  rename(Species = species) %>%
  group_by(Species) %>%
  summarise(across(everything(), sum))

group_data <- merge(real, sim, by = "Species")

df_r2 <- data.frame(matrix(ncol = 0, nrow = 1))
rownames(df_r2) <- c("R2")


df_r2 <- data.frame(matrix(ncol = 0, nrow = 1))
rownames(df_r2) <- c("R2")

if (ncol(group_data) > 2) {
  for (i in 3:ncol(group_data)) {
    model <- lm(group_data[, i] ~ real_abundance, data = group_data)
    r2 <- summary(model)$r.squared
    df_r2[[colnames(group_data)[i]]] <- r2
  }
} else {
  r2 <- c(NA)
  return(NA)
}

write.csv(df_r2, paste("R2_", original_name, ".csv", sep=""), row.names = TRUE)
