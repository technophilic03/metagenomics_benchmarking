suppressPackageStartupMessages(library(tidyverse))

real_files <- list.files(path      = "MeSS_code/ground_truth/updated/",
                         pattern   = "\\.csv$",
                         full.names = TRUE)

sim_files <- list.files(path      = "profilers/Kraken2/combined_files/updated_species_names/updated/",
                        pattern   = "\\.csv$",
                        full.names = TRUE)

process_real_temp_dir <- "MeSS_code/ground_truth/updated/temp/"
process_sim_temp_dir <- "profilers/Kraken2/combined_files/updated_species_names/analysis_results/temp/"
analysis_res_dir <- "profilers/Kraken2/combined_files/updated_species_names/analysis_results/"

# create folder if not exist
dirs_to_create <- c(process_real_temp_dir, process_sim_temp_dir, analysis_res_dir)
sapply(dirs_to_create, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

process_real_data <- function(real_data_file) {
  message("processing ground truth (real data)")
  real_data_raw <- read.csv(real_data_file)
  
  real_data <- data.frame(species = real_data_raw$species,
                          abundance = real_data_raw$tax_abundance)
  
  grouped_data <- aggregate(abundance ~ species, real_data, sum)
  
  write.csv(
    grouped_data,
    paste0(process_real_temp_dir, "processed_real_data.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  message(
    "real_data saved to",
    process_real_temp_dir,
    "\n"
  )
}

process_sim_data <- function(sim_data_file) {
  sim_data_raw <- read.csv(sim_data_file)
  
  sim_data <- data.frame(species = sim_data_raw$Name, abundance_cols = sim_data_raw[, 10:ncol(sim_data_raw)])
  
  out_dir <- process_sim_temp_dir
  out_name <- paste0("processed_", basename(sim_data_file))
  out_path <- file.path(out_dir, out_name)
  
  write.csv(sim_data, out_path, row.names = FALSE, quote = FALSE)
}

getR2Function <- function(real_data_file, sim_data_file) {
  # real_data_file = processed_real_data.csv
  # sim_data_file = processed_sim_data.csv
  real_data <- read.csv(real_data_file)
  sim_data <- read.csv(sim_data_file)
  
  processed_real_data <- real_data |>
    rename(real_abundance = abundance)
  
  processed_sim_data <- sim_data |>
    group_by(species) |> # group_by(Name); after process_sim_data function, species == Name
    summarise(across(everything(), sum))
  
  grouped_data_r2 <- full_join(processed_real_data, processed_sim_data, by = "species")
  grouped_data_r2$real_abundance[is.na(grouped_data_r2$real_abundance)] <- 0
  grouped_data_r2[, 3:ncol(grouped_data_r2)][is.na(grouped_data_r2[, 3:ncol(grouped_data_r2)])] <- 0
  
  output_r2 <- data.frame(matrix(ncol = 0, nrow = 1))
  rownames(output_r2) <- c("R2")
  
  for (i in 3:ncol(grouped_data_r2)) {
    model <- lm(real_abundance ~ grouped_data_r2[, i], data = grouped_data_r2)
    r2 <- summary(model)$r.squared
    output_r2[[colnames(grouped_data_r2)[i]]] <- r2
  }
  
  out_dir <- analysis_res_dir
  out_name <- paste0("R2_", basename(sim_data_file))
  out_path <- file.path(out_dir, out_name)
  write.csv(output_r2, file = out_path, row.names = FALSE)
  message("R2 result saved to ", out_dir, "\n")
}

RMSFunction <- function(real_data_file, sim_data_file) {
  # real_data_file = processed_real_data.csv
  # sim_data_file = processed_sim_data.csv
  real_data <- read.csv(real_data_file)
  sim_data <- read.csv(sim_data_file)
  
  processed_real_data <- real_data |>
    rename(real_abundance = abundance)
  
  processed_sim_data <- sim_data |>
    group_by(species) |> # group_by(Name); after process_sim_data function, species == Name
    summarise(across(everything(), sum))
  
  grouped_data_rms <- full_join(processed_real_data, processed_sim_data, by = "species")
  grouped_data_rms$real_abundance[is.na(grouped_data_rms$real_abundance)] <- 0
  grouped_data_rms[, 3:ncol(grouped_data_rms)][is.na(grouped_data_rms[, 3:ncol(grouped_data_rms)])] <- 0
  
  output_rms <- data.frame(matrix(ncol = 0, nrow = 3))
  rownames(output_rms) <- c("RRMSE", "AVGRE", "NRMSE")
  
  # observed -> sim
  # predicted -> real
  get_rrmse <- function(observed, predicted) {
    ## exclude both 0
    valid_indices <- !(observed == 0 & predicted == 0)
    observed <- observed[valid_indices]
    predicted <- predicted[valid_indices]
    
    squared_diff <- (observed - predicted) ^ 2
    rmse <- sqrt(mean(squared_diff))
    mean_observed <- mean(observed)
    rrmse <- (rmse / mean_observed) * 100 # mean_observed becomes small if many 0 in observed;
    return(rrmse)
  }
  
  get_nrmse <- function (observed, predicted) {
    ## exclude both 0
    valid_indices <- !(observed == 0 & predicted == 0)
    observed <- observed[valid_indices]
    predicted <- predicted[valid_indices]
    
    squared_diff <- (observed - predicted) ^ 2
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
  
  
  out_dir <- analysis_res_dir
  out_name <- paste0("RRMSE_", basename(sim_data_file))
  out_path <- file.path(out_dir, out_name)
  write.csv(output_rms, file = out_path, row.names = TRUE)
  message("RRMSE result saved to ", out_dir, "\n")
}


ROCFunction <- function(real_data_file,
                        sim_data_file,
                        total_genomes) {
  
  # total_genomes = 21418 for Metascope/Pathoscope2
  # total_genomes = 19987 for Bracken/Kraken
  
  
  # real_data_file = processed_real_data.csv
  # sim_data_file = processed_sim_data.csv
  real_data <- read.csv(real_data_file)
  sim_data <- read.csv(sim_data_file)
  
  processed_real_data <- real_data |>
    rename(real_abundance = abundance)
  
  processed_sim_data <- sim_data |>
    group_by(species) |> # group_by(Name); after process_sim_data function, species == Name
    summarise(across(everything(), sum))
  
  grouped_data_roc <- full_join(processed_real_data, processed_sim_data, by = "species")
  grouped_data_roc$real_abundance[is.na(grouped_data_roc$real_abundance)] <- 0
  grouped_data_roc[, 3:ncol(grouped_data_roc)][is.na(grouped_data_roc[, 3:ncol(grouped_data_roc)])] <- 0
  
  output_roc <- data.frame(matrix(ncol = 0, nrow = 2))
  rownames(output_roc) <- c("TPR", "FPR")
  
  for (i in 3:ncol(grouped_data_roc)) {
    TP <- 0
    FP <- 0
    FN <- 0
    TN <- as.numeric(total_genomes) - nrow(grouped_data_roc)
    
    for (j in 1:nrow(grouped_data_roc)) {
      predicted <- grouped_data_roc$real_abundance[j] # ground_truth
      observed <- grouped_data_roc[j, i] # sim
      
      if (predicted >= observed * 0.5) {
        TP <- TP + 1 # TP when sim >= real * 0.5
      } else {
        if (predicted == 0 & observed > 0){
          FP <- FP + 1
        } 
      }
      
      if (predicted > 0 & observed == 0){
        FN <- FN + 1
      }
    } 
    

    
    get_TPR <- TP / (TP + FN)
    get_FPR <- FP / (FP + TN)
    output_roc[[colnames(grouped_data_roc)[i]]] <- c(get_TPR, get_FPR)
  }
  
  
  t_output_roc <- as.data.frame(t(output_roc))
  
  out_dir <- analysis_res_dir
  out_name <- paste0("ROC_", basename(sim_data_file))
  out_path <- file.path(out_dir, out_name)
  write.csv(t_output_roc, file = out_path, row.names = TRUE)
  message("ROC result saved to ", out_dir, "\n")
}

#
### Run Function

process_real_data(real_files[1]) # only need to run once

real_data_dir <- "MeSS_code/ground_truth/updated/temp/processed_real_data.csv"

process_file <- function(sim_file) {
   message("Processing: ", basename(sim_file))

   process_sim_data(sim_file)

   sim_processed_path <- file.path(
     process_sim_temp_dir,
     paste0("processed_", basename(sim_file))
   )

   getR2Function(
     real_data_dir,
     sim_processed_path
   )

   RMSFunction(
     real_data_dir,
     sim_processed_path
   )

   ROCFunction(
     real_data_dir,
     sim_processed_path,
     total_genomes = 19987
   )
 }

 walk(sim_files, process_file)
 
 # cleanup temp files at the end
 
 
 
 
