suppressPackageStartupMessages(
    library(tidyverse)
)

real_files <- list.files(
    path      = "MeSS_code/ground_truth/updated/",
    pattern   = "\\.csv$",
    full.names= TRUE
)

sim_files <- list.files(
    path      = "profilers/Bracken/combined_files/name_update/updated/",
    pattern   = "\\.csv$",
    full.names= TRUE
)


process_real_data <- function(real_data_file) {
    message("processing groud truth (real data)")
    real_data_raw <- read.csv(real_data_file)

    real_data <- data.frame(
        species = real_data_raw$species,
        abundance = real_data_raw$tax_abundance
    )

    grouped_data <- aggregate(abundance ~ species, real_data, sum)

    write.csv(
        grouped_data,
        "MeSS_code/ground_truth/updated/temp/processed_real_data.csv",
        row.names = FALSE,
        quote = FALSE
        )
    message("real_data saved to MeSS_code/ground_truth/updated/temp/processed_real_data.csv\n")
}

process_sim_data <- function(sim_data_file) {
    sim_data_raw <- read.csv(sim_data_file)

    sim_data <- data.frame(
        species = sim_data_raw$Name,
        abundence_cols = sim_data_raw[, 10:ncol(sim_data_raw)]
    )
    
    out_dir<- "profilers/Bracken/combined_files/name_update/updated/analysis_results/temp/"
    out_name <- paste0("processed_", basename(sim_data_file))
    out_path <- file.path(out_dir, out_name)
    
    write.csv(
        sim_data,
        out_path,
        row.names = FALSE,
        quote = FALSE
    )
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
        model <- lm(real_abundance ~ grouped_data_r2[, i], 
                    data = grouped_data_r2)
        r2 <- summary(model)$r.squared
        output_r2[[colnames(grouped_data_r2)[i]]] <- r2
    }

    out_dir<- "profilers/Bracken/combined_files/name_update/updated/analysis_results/"
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
    
    output_rms <- data.frame(matrix(ncol = 0, nrow = 2))
    rownames(output_rms) <- c("RRMSE", "AVGRE")

    # observed -> sim
    # predicted -> real
    get_rrmse <- function(observed, predicted) {
        squared_diff <- (observed - predicted) ^ 2
        mse <- mean(squared_diff)
        rmse <- sqrt(mse)
        mean_observed <- mean(observed)
        rrmse <- (rmse / mean_observed) * 100 # mean_observed becomes small if many 0 in observed; 
        # TODO: should we exclude 0?
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

    if (ncol(grouped_data_rms) > 2) {
        for (i in 3:ncol(grouped_data_rms)) {
            rrmse <- get_rrmse(grouped_data_rms[, i], grouped_data_rms$real_abundance)
            avgre <- get_avgre(grouped_data_rms[, i], grouped_data_rms$real_abundance)
            output_rms[[colnames(grouped_data_rms)[i]]] <- c(rrmse, avgre)
        }
    } else {
        return(NA)
    }
    
    out_dir<- "profilers/Bracken/combined_files/name_update/updated/analysis_results/"
    out_name <- paste0("RRMSE_", basename(sim_data_file))
    out_path <- file.path(out_dir, out_name)
    write.csv(output_rms, file = out_path, row.names = TRUE)
    message("RRMSE result saved to ", out_dir, "\n")
}


ROCFunction <- function(real_data_file, sim_data_file, total_genomes) {
    
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
    

    TP <- 0
    
    FP <- sum(grouped_data_roc$real_abundance == 0 & rowSums(grouped_data_roc[, 3:ncol(grouped_data_roc)]) > 0)
    FN <- sum(grouped_data_roc$real_abundance > 0 & rowSums(grouped_data_roc[, 3:ncol(grouped_data_roc)]) == 0)
    
    TN <- as.numeric(total_genomes) - nrow(grouped_data_roc)
    
    output_roc <- data.frame(matrix(ncol = 0, nrow = 2))
    rownames(output_roc) <- c("TPR", "FPR")
    
    if (ncol(grouped_data_roc) > 2) {
        for (i in 3:ncol(grouped_data_roc)) {
            TP_col <- TP
            FP_col <- FP
            
            for (j in 1:nrow(grouped_data_roc)) {
                if (grouped_data_roc[j, i] >= grouped_data_roc$real_abundance[j] * 0.5) {
                    TP_col <- TP_col + 1 # TP when sim >= real * 0.5
                    # TODO: do we want to apply a threshold?
                } 
            }
            
            get_TPR <- TP_col / (TP_col + FN)
            get_FPR <- FP_col / (FP_col + TN)
            df_roc[[colnames(grouped_data_roc)[i]]] <- c(get_TPR, get_FPR)
        }
    } else {
        return(NA)
    }

    t_output_roc <- as.data.frame(t(output_roc))
    
    out_dir<- "profilers/Bracken/combined_files/name_update/updated/analysis_results/"
    out_name <- paste0("ROC_", basename(sim_data_file))
    out_path <- file.path(out_dir, out_name)
    write.csv(output_ROC, file = out_path, row.names = TRUE)
    message("ROC result saved to ", out_dir, "\n")
}

#
### Run Function

# process_real_data(real_files[1]) # only need to run once


process_file <- function(sim_file) {
    message("Processing: ", basename(sim_file))
    
    process_sim_data(sim_file)
    
    sim_processed_path <- file.path(
        "profilers/Bracken/combined_files/name_update/updated/analysis_results/temp/", 
        paste0("processed_", basename(sim_file))
    )
    
    getR2Function(
        "MeSS_code/ground_truth/updated/temp/processed_real_data.csv",
        sim_processed_path
    )
    
    RMSFunction(
        "MeSS_code/ground_truth/updated/temp/processed_real_data.csv",
        sim_processed_path
    )
    
    ROCFunction(
        "MeSS_code/ground_truth/updated/temp/processed_real_data.csv",
        sim_processed_path
    )
}

walk(sim_files, process_file)

