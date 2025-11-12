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

# Ground truth settings
real_files <- list.files(
  path = "MeSS_code/ground_truth/updated/",
  pattern = "\\.csv$",
  full.names = TRUE
)

process_real_temp_dir <- "MeSS_code/ground_truth/updated/temp/"

# Create ground truth temp folder
dir.create(process_real_temp_dir, recursive = TRUE, showWarnings = FALSE)

process_real_data <- function(real_data_file) {
  message("Processing ground truth (real data)")
  real_data_raw <- read.csv(real_data_file)
  
  real_data <- data.frame(
    species = real_data_raw$species,
    abundance = real_data_raw$tax_abundance, 
    sample = real_data_raw$sample
  )
  
  grouped_data <- real_data |> 
    mutate(abundance = abundance / 100)
  
  write.csv(
    grouped_data,
    paste0(process_real_temp_dir, "processed_real_data.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  message("Real data saved to ", process_real_temp_dir, "\n")
}

process_sim_data <- function(sim_data_file, process_sim_temp_dir) {
  sim_data_raw <- read.csv(sim_data_file)
  
  sim_data <- data.frame(
    species = sim_data_raw$Name,
    sim_data_raw[, 10:ncol(sim_data_raw)]
  )
  
  out_name <- paste0("processed_", basename(sim_data_file))
  out_path <- file.path(process_sim_temp_dir, out_name)
  
  write.csv(sim_data, out_path, row.names = FALSE, quote = FALSE)
}

getR2Function <- function(sim_data_file, analysis_res_dir) {
  real_data <- read.csv("MeSS_code/ground_truth/updated/temp/processed_real_data.csv")
  sim_data <- read.csv(sim_data_file)
  
  processed_real_data <- real_data |>
    mutate(real_abundance = abundance / 100) |>
    select(species, real_abundance, sample) |>
    filter(!grepl("_500_", sample))
    
  
  processed_sim_data <- sim_data |>
    group_by(Name) |>
    summarise(across(everything(), sum)) |>
    mutate(across(-Name, ~ .x / sum(.x))) |>
    pivot_longer(cols = -c(species), values_to = "sim_read_counts", names_to = "sample")
  
  grouped_data_r2 <- full_join(processed_real_data, processed_sim_data, by = c("species", "sample"))
  grouped_data_r2$real_abundance[is.na(grouped_data_r2$real_abundance)] <- 0
  grouped_data_r2$sim_read_counts[is.na(grouped_data_r2$sim_read_counts)] <- 0
  
  
  output_r2 <- grouped_data_r2 |> group_by(sample) |>
    do({
      model <- lm(sim_read_counts ~ real_abundance, data = .)
      data.frame(r2 = summary(model)$r.squared)
    })
  
  out_name <- paste0("R2_", basename(sim_data_file))
  out_path <- file.path(analysis_res_dir, out_name)
  write.csv(output_r2, file = out_path, row.names = FALSE)
  message("R2 result saved to ", analysis_res_dir, "\n")
}

RMSFunction <- function(real_data_file, sim_data_file, analysis_res_dir) {
  real_data <- read.csv(real_data_file)
  sim_data <- read.csv(sim_data_file)
  
  processed_real_data <- real_data |>
    rename(real_abundance = abundance)
  
  processed_sim_data <- sim_data |>
    group_by(species) |>
    summarise(across(everything(), sum))
  
  grouped_data_rms <- full_join(processed_real_data, processed_sim_data, by = "species")
  grouped_data_rms$real_abundance[is.na(grouped_data_rms$real_abundance)] <- 0
  grouped_data_rms[, 3:ncol(grouped_data_rms)][is.na(grouped_data_rms[, 3:ncol(grouped_data_rms)])] <- 0
  
  output_rms <- data.frame(matrix(ncol = 0, nrow = 3))
  rownames(output_rms) <- c("RRMSE", "AVGRE", "NRMSE")
  
  get_rrmse <- function(observed, predicted) {
    valid_indices <- !(observed == 0 & predicted == 0)
    observed <- observed[valid_indices]
    predicted <- predicted[valid_indices]
    
    squared_diff <- (observed - predicted)^2
    rmse <- sqrt(mean(squared_diff))
    mean_observed <- mean(observed)
    rrmse <- (rmse / mean_observed) * 100
    return(rrmse)
  }
  
  get_nrmse <- function(observed, predicted) {
    valid_indices <- !(observed == 0 & predicted == 0)
    observed <- observed[valid_indices]
    predicted <- predicted[valid_indices]
    
    squared_diff <- (observed - predicted)^2
    normalizer <- ((observed + predicted) / 2)
    
    nrmse <- sqrt(mean(squared_diff / normalizer)) * 100
    return(nrmse)
  }
  
  get_avgre <- function(observed, predicted) {
    valid_indices <- observed != 0 & !(observed == 0 & predicted == 0)
    observed <- observed[valid_indices]
    predicted <- predicted[valid_indices]
    
    avgre <- mean(abs(observed - predicted) / observed) * 100
    return(avgre)
  }
  
  for (i in 3:ncol(grouped_data_rms)) {
    rrmse <- get_rrmse(grouped_data_rms[, i], grouped_data_rms$real_abundance)
    nrmse <- get_nrmse(grouped_data_rms[, i], grouped_data_rms$real_abundance)
    avgre <- get_avgre(grouped_data_rms[, i], grouped_data_rms$real_abundance)
    output_rms[[colnames(grouped_data_rms)[i]]] <- c(rrmse, avgre, nrmse)
  }
  
  out_name <- paste0("RRMSE_", basename(sim_data_file))
  out_path <- file.path(analysis_res_dir, out_name)
  write.csv(output_rms, file = out_path, row.names = TRUE)
  message("RRMSE result saved to ", analysis_res_dir, "\n")
}

ROCFunction <- function(real_data_file, sim_data_file, total_genomes, analysis_res_dir) {
  real_data <- read.csv(real_data_file)
  sim_data <- read.csv(sim_data_file)

  processed_real_data <- real_data %>%
    group_by(species) %>%
    summarise(real_abundance = sum(abundance), .groups = "drop")

  processed_sim_data <- sim_data |>
    group_by(species) |>
    summarise(across(everything(), sum))

  grouped_data_roc <- full_join(processed_real_data, processed_sim_data, by = "species")
  grouped_data_roc$real_abundance[is.na(grouped_data_roc$real_abundance)] <- 0
  grouped_data_roc[, 3:ncol(grouped_data_roc)][is.na(grouped_data_roc[, 3:ncol(grouped_data_roc)])] <- 0

  output_roc <- data.frame(matrix(ncol = 0, nrow = 6))
  rownames(output_roc) <- c("TP", "TN", "FP", "FN", "TPR", "FPR")

  for (i in 3:ncol(grouped_data_roc)) {
    TP <- 0
    FP <- 0
    FN <- 0
    TN <- as.numeric(total_genomes) - nrow(grouped_data_roc)

    for (j in 1:nrow(grouped_data_roc)) {
      predicted <- grouped_data_roc$real_abundance[j]
      observed <- grouped_data_roc[j, i]

      if (predicted > 0 & observed > 0) {
        TP <- TP + 1
      } else {
        if (predicted == 0 & observed > 0) {
          FP <- FP + 1
        }
      }

      if (predicted > 0 & observed == 0) {
        FN <- FN + 1
      }
    }

    get_TPR <- TP / (TP + FN)
    get_FPR <- FP / (FP + TN)
    output_roc[[colnames(grouped_data_roc)[i]]] <- c(TP, TN, FP, FN, get_TPR, get_FPR)
  }

  t_output_roc <- as.data.frame(t(output_roc))

  out_name <- paste0("ROC_", basename(sim_data_file))
  out_path <- file.path(analysis_res_dir, out_name)
  write.csv(t_output_roc, file = out_path, row.names = TRUE)
  message("ROC result saved to ", analysis_res_dir, "\n")
}


# MCC
MCCFunction <- function(real_data_file, sim_data_file, total_genomes, analysis_res_dir) {
  real_data <- read.csv(real_data_file)
  sim_data <- read.csv(sim_data_file)

  processed_real_data <- real_data |>
    rename(real_abundance = abundance)

  processed_sim_data <- sim_data |>
    group_by(species) |>
    summarise(across(everything(), sum))

  grouped_data_roc <- full_join(processed_real_data, processed_sim_data, by = "species")
  grouped_data_roc$real_abundance[is.na(grouped_data_roc$real_abundance)] <- 0
  grouped_data_roc[, 3:ncol(grouped_data_roc)][is.na(grouped_data_roc[, 3:ncol(grouped_data_roc)])] <- 0

  output_MCC <- data.frame(matrix(ncol = 0, nrow = 1))
  rownames(output_MCC) <- c("MCC")

  for (i in 3:ncol(grouped_data_roc)) {
    TP <- 0
    FP <- 0
    FN <- 0
    TN <- as.numeric(total_genomes) - nrow(grouped_data_roc)

    for (j in 1:nrow(grouped_data_roc)) {
      predicted <- grouped_data_roc$real_abundance[j]
      observed <- grouped_data_roc[j, i]

      if (predicted > 0 & observed > 0) {
        TP <- TP + 1
      } else {
        if (predicted == 0 & observed > 0) {
          FP <- FP + 1
        }
      }

      if (predicted > 0 & observed == 0) {
        FN <- FN + 1
      }
    }

    get_MCC <- (TP * TN - FP * FN) / (sqrt((TP + TN)*(TP + FN)*(TN + FP)*(TN+FN)))
    output_MCC[[colnames(grouped_data_roc)[i]]] <- c(get_MCC)
  }

  t_output_MCC <- as.data.frame(t(output_MCC))

  out_name <- paste0("MCC_", basename(sim_data_file))
  out_path <- file.path(analysis_res_dir, out_name)
  write.csv(t_output_MCC, file = out_path, row.names = TRUE)
  message("MCC result saved to ", analysis_res_dir, "\n")
}

process_file <- function(sim_file, process_sim_temp_dir, analysis_res_dir, 
                         real_data_dir, total_genomes) {
  message("  Processing: ", basename(sim_file))
  
  process_sim_data(sim_file, process_sim_temp_dir)
  
  sim_processed_path <- file.path(
    process_sim_temp_dir,
    paste0("processed_", basename(sim_file))
  )
  
  getR2Function(sim_processed_path, analysis_res_dir)
  #RMSFunction(real_data_dir, sim_processed_path, analysis_res_dir)
  #ROCFunction(real_data_dir, sim_processed_path, total_genomes, analysis_res_dir)
  #MCCFunction(real_data_dir, sim_processed_path, total_genomes, analysis_res_dir)
}




# Process ground truth data (only need to run once)
process_real_data(real_files[1])
real_data_dir <- paste0(process_real_temp_dir, "processed_real_data.csv")

profilers <- c("Bracken", "Centrifuge", "Kraken2", "MetaScope", 
               "MetaScope_blast", "MetaScope_priors", "PathoScope2", "mOTUs")



#profilers <- c("mOTUs")
# Process each profiler
for (profiler in profilers) {
  message("Processing profiler: ", profiler)

  # Define profiler-specific paths
  sim_files <- list.files(
    path = paste0("profilers/", profiler, "/combined_files/updated_species_names/updated/"),
    pattern = "\\.csv$",
    full.names = TRUE
  )
  
  process_sim_temp_dir <- paste0(
    "profilers/", profiler, "/combined_files/updated_species_names/analysis_results/temp/"
  )
  
  analysis_res_dir <- paste0(
    "profilers/", profiler, "/combined_files/updated_species_names/analysis_results/"
  )
  
  # Create directories
  dir.create(process_sim_temp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(analysis_res_dir, recursive = TRUE, showWarnings = FALSE)
  
  # # Get total genomes for this profiler
  total_genomes <- total_genomes_map[[profiler]]

  if (is.null(total_genomes)) {
    warning("Total genomes not defined for ", profiler)
    next
  }

  
  # Process all files for this profiler
  for (sim_file in sim_files) {
    process_file(
      sim_file = sim_file,
      process_sim_temp_dir = process_sim_temp_dir,
      analysis_res_dir = analysis_res_dir,
      real_data_dir = real_data_dir, 
      total_genomes = total_genomes
    )
  }
  message("\nCompleted processing for ", profiler, "\n")
}


# cleanup temp files
# unlink(process_real_temp_dir, recursive = TRUE)
# for (profiler in profilers) {
#   temp_dir <- paste0("profilers/", profiler, "/combined_files/updated_species_names/analysis_results/temp/")
#   unlink(temp_dir, recursive = TRUE)
# }