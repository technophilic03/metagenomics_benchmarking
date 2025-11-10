# Generate Analysis Metrics Files
# Sean Lu


suppressPackageStartupMessages(library(tidyverse))

# Define profilers to process
profilers <- c("Bracken", "Centrifuge", "Kraken2", "MetaScope", 
               "MetaScope_blast", "MetaScope_priors", "mOTUs", "PathoScope2")

# Define total genomes for each profiler (for ROC calculations)
total_genomes_map <- list(
  "Centrifuge" = 8425,
  "MetaScope" = 22951,
  "MetaScope_blast"= 22951,
  "MetaScope_priors" = 22951,
  "PathoScope2" = 21418,
  "Bracken" = 19987,
  "Kraken2" = 19987,
  "mOTUs" = 34341
)

# Ground Truth
processed_real_data <- read.csv("MeSS_code/ground_truth/ground_truth.csv") |>
  select(species, tax_abundance, sample) |>  
  mutate(real_abundance = tax_abundance / 100) |>
  select(species, real_abundance, sample) |>
  filter(!grepl("_500_", sample)) |> 
  unique()


# Get Sample Name Helper Function
getSampleName <- function(sample_path) {
  patterns <- c("100k_err", "100k_noerr", "1mil_err", "1mil_noerr", "10mil_err", "10mil_noerr")
  idx <- which(str_detect(sample_path, fixed(patterns)))
  if (length(idx)) patterns[idx[1]] else NA_character_
}
  
# Get R2 
getR2Function <- function(sim_data_file, analysis_res_dir, profiler) {

  sim_data_raw <- read.csv(sim_data_file)
  
  sim_data <- data.frame(
    species = sim_data_raw$Name,
    sim_data_raw[, 10:ncol(sim_data_raw)]
  )
  sim_data[is.na(sim_data)] <- 0
  
  processed_sim_data <- sim_data |>
    group_by(species) |>
    summarise(across(everything(), sum)) |>
    mutate(across(-species, ~ .x / sum(.x))) |>
    pivot_longer(cols = -c(species), values_to = "sim_read_counts", names_to = "sample") |>
    mutate(sample = tolower(sample)) |>
    mutate(sample = gsub("_noerr", "", sample))
  
  
  grouped_data <- full_join(processed_real_data, processed_sim_data, by = c("species", "sample"))
  grouped_data$real_abundance[is.na(grouped_data$real_abundance)] <- 0
  grouped_data$sim_read_counts[is.na(grouped_data$sim_read_counts)] <- 0
  
  
  output_r2 <- grouped_data |> group_by(sample) |>
    do({
      model <- lm(sim_read_counts ~ real_abundance, data = .)
      data.frame(r2 = summary(model)$r.squared)
    })
  
  get_rmse <- function(observed, predicted) {
    squared_diff <- (observed - predicted)^2
    rmse <- sqrt(mean(squared_diff))
    return(rmse)
  }
  get_rrmse <- function(observed, predicted) {
    squared_diff <- (observed - predicted)^2
    rmse <- sqrt(mean(squared_diff))
    mean_observed <- mean(observed)
    rrmse <- (rmse / mean_observed)
    return(rrmse)
  }
  
  get_nrmse <- function(observed, predicted) {
    squared_diff <- (observed - predicted)^2
    rmse <- sqrt(mean(squared_diff))
    nrmse <- rmse / sd(observed)
    return(nrmse)
  }
  
  get_avgre <- function(observed, predicted) {
    valid_indices <- !(predicted == 0)
    observed <- observed[valid_indices]
    predicted <- predicted[valid_indices]
    
    avgre <- mean(abs(observed - predicted) / predicted)
    return(avgre)
  }
  
  output_rms <- grouped_data |> 
    group_by(sample) |>
    summarize(
      rmse = get_rmse(sim_read_counts, real_abundance),
      rrmse = get_rrmse(sim_read_counts, real_abundance),
      nrmse = get_nrmse(sim_read_counts, real_abundance),
      avgre = get_avgre(sim_read_counts, real_abundance))
  
  output_roc <- grouped_data |> 
    mutate(TP = (real_abundance > 0 & sim_read_counts > 0), 
           FP = (real_abundance == 0 & sim_read_counts > 0),
           FN = (real_abundance > 0 & sim_read_counts == 0)) |>
    group_by(sample) |> 
    summarize(TP = sum(TP),
              FP = sum(FP), 
              FN = sum(FN), 
              TN = total_genomes - FN, 
              TPR = TP / (TP + FN), 
              FPR = FP / (FP + TN), 
              MCC = (TP * TN - FP * FN) / (sqrt((TP + TN)*(TP + FN)*(TN + FP)*(TN+FN))))
  
  analysis_out <- left_join(left_join(output_r2, output_rms, by = "sample"), output_roc, by = "sample") 
  print(analysis_out)
  
  out_name <- paste0(profiler, "_", getSampleName(sim_data_file), ".csv")
  out_path <- file.path(analysis_res_dir, out_name)
  write.csv(analysis_out, file = out_path, row.names = FALSE)
  message("R2 result saved to ", out_path, "\n")
}


# Process each profiler
for (profiler in profilers) {
  message("Processing profiler: ", profiler)

  # Define profiler-specific paths
  sim_files <- list.files(
    path = paste0("profilers/", profiler, "/combined_files/updated_species_names/updated"),
    pattern = "\\.csv$",
    full.names = TRUE
  )
  
  analysis_res_dir <- paste0(
    "results/analysis_results"
  )
  dir.create(analysis_res_dir, recursive = TRUE, showWarnings = FALSE)
  
  # # Get total genomes for this profiler
  total_genomes <- total_genomes_map[[profiler]]

  if (is.null(total_genomes)) {
    warning("Total genomes not defined for ", profiler)
    next
  }

  
  # Process all files for this profiler
  for (sim_file in sim_files) {
    getR2Function(sim_file, analysis_res_dir, profiler)
  }
  message("\nCompleted processing for ", profiler, "\n")
}
