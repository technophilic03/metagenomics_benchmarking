library(tidyverse)
run_analysis <- function(real_data_file, sim_data_file, analysis_methods = c("R2", "RRMSE", "ROC"), total_genomes = NULL) {

    cat("Processing data files...\n")
    process_real_data(real_data_file)
    process_sim_data(sim_data_file)

    original_name <- tools::file_path_sans_ext(basename(sim_data_file))

    results <- list()

    for (method in analysis_methods) {
        cat("Running", method, "analysis...\n")

        if (method == "R2") {
            results$R2 <- getR2Function("processed_real_data.csv",
                                       "processed_sim_data.csv",
                                       original_name)
        } else if (method == "RRMSE") {
            results$RRMSE <- RMSFunction("processed_real_data.csv",
                                             "processed_sim_data.csv",
                                             original_name)
        } else if (method == "ROC") {
            results$SIMOBS <- ROCFunction("processed_real_data.csv",
                                              "processed_sim_data.csv",
                                              total_genomes <- 1000,
                                              original_name)
        }
    }

    cat("Analysis complete!\n")
    return(results)
}

process_real_data <- function(real_data_file) {
    real_data_raw <- read.csv(real_data_file)

    real_data <- data.frame(
        species = real_data_raw$species,
        abundance = as.numeric(real_data_raw$tax_abundance)
    )

    grouped_data <- aggregate(abundance ~ species, real_data, sum)

    write.csv(
        grouped_data,
        "processed_real_data.csv",
        row.names = FALSE,
        quote = FALSE
        )
}

process_sim_data <- function(sim_data_file) {
    sim_data_raw <- read.csv(sim_data_file)

    sim_data <- data.frame(
        species = sim_data_raw$Species,
        abundence_cols <- sim_data_raw[, 10:ncol(sim_data_raw)]
    )

    write.csv(
        sim_data,
        "processed_sim_data.csv",
        row.names = FALSE,
        quote = FALSE
    )
}

getR2Function <- function(real_data, sim_data, original_name) {

    real_data <- read.csv(real_data)
    sim_data <- read.csv(sim_data)

    real <- real_data |>
        rename(Species = species, real_abundance = abundance)

    sim <- sim_data |>
        rename(Species = species) |>
        group_by(Species) |>
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

    result_file <- paste0("R2_", original_name, ".csv")
    write.csv(df_r2, result_file, row.names = TRUE)
    return(df_r2)
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

