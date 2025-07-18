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

    real_data <- read.csv(real_data_file)
    sim_data <- read.csv(sim_data_file)

    
    processed_real_data <- real_data |>
        rename(real_abundance = abundance)

    processed_sim_data <- sim_data |>
        group_by(species) |> # group_by(Name); after process_sim_data function, species == Name
        summarise(across(everything(), sum))

    grouped_data_r2 <- processed_real_data |>
        left_join(processed_sim_data, by = "species")
    
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

RMSFunction <- function(real_data, sim_data, original_name) {
    real_data <- read.csv(real_data)
    sim_data <- read.csv(sim_data)

    real <- real_data |>
        rename(Species = species, real_abundance = abundance)

    sim <- sim_data |>
        rename(Species = species) |>
        group_by(Species) |>
        summarise(across(everything(), sum))

    group_data <- merge(real, sim, by = "Species")

    df_rms <- data.frame(matrix(ncol = 0, nrow = 2))
    rownames(df_rms) <- c("RRMSE", "AVGRE")

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

    if (ncol(group_data) > 2) {
        for (i in 3:ncol(group_data)) {
            rrmse <- get_rrmse(group_data[, i], group_data$real_abundance)
            avgre <- get_avgre(group_data[, i], group_data$real_abundance)
            df_rms[[colnames(group_data)[i]]] <- c(rrmse, avgre)
        }
    } else {
        return(NA)
    }

    result_file <- paste0("RMS_", original_name, ".csv")
    write.csv(df_rms, result_file, row.names = TRUE)
    return(df_rms)

}

ROCFunction <- function(real_data, sim_data, total_genomes, original_name) {
    real_data <- read.csv(real_data)
    sim_data <- read.csv(sim_data)

    real <- real_data |>
        rename(Species = species, real_abundance = abundance)

    sim <- sim_data |>
        rename(Species = species) |>
        group_by(Species) |>
        summarise(across(everything(), sum))

    group_data <- merge(real, sim, by = "Species")

    missing_in_real <- anti_join(sim, real, by = "Species")
    missing_in_sim <- anti_join(real, sim, by = "Species")

    df_roc <- data.frame(matrix(ncol = 0, nrow = 2))
    rownames(df_roc) <- c("TPR", "FPR")

    TP <- 0
    FP <- 0
    FN <- 0
    FP <- FP + nrow(missing_in_real) # FP when a sequence found in the simulation isn't in the real data
    FN <- FN + nrow(missing_in_sim) # FN when sequences in real data but missing in simulation
    TN <- as.numeric(total_genomes) - nrow(sim) # TN = total genomes - total simulated entries

    if (ncol(group_data) > 2) {
        for (i in 3:ncol(group_data)) {
            TP_col <- TP
            FP_col <- FP

            for (j in 1:nrow(group_data)) {
                if (group_data[j, i] >= group_data$real_abundance[j] * 0.5) {
                    TP_col <- TP_col + 1 # TP + 1 when sim_read >= real_reads/2
                } else if (group_data[j, i] == 0 && group_data$real_abundance[j] > 0) {
                    FP_col <- FP_col + 1 # FP +1 when simulated reads = 0 but real reads >0
                }
            }

            get_TPR <- TP_col / (TP_col + FN)
            get_FPR <- FP_col / (FP_col + TN)
            df_roc[[colnames(group_data)[i]]] <- c(get_TPR, get_FPR)
        }
    } else {
        return(NA)
    }

    t_df_roc <- as.data.frame(t(df_roc))

    result_file <- paste0("ROC_", original_name, ".csv")
    write.csv(t_df_roc, result_file, row.names = TRUE)
    return(t_df_roc)
}


### Run Function
process_real_data(real_files[1])

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
}

walk(sim_files, process_file)

