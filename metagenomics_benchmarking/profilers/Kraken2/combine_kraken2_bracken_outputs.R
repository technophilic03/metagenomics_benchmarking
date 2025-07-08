# R Script to combine outputs - Kraken 2
# Metagenomics benchmarking
# Aubrey Odom
# 6/11/25

# Description ----

# This script will take all Kraken outputs,
# consolidate them, and create a formatted
# table to use in downstream analysis

## IMPORTANT
## MUST RUN combine_outputs_kraken2.sh FIRST!!!

# Setup ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(LegATo)
})

# Set file paths ----
stem <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark"
out_path <- file.path(stem, "profilers/Kraken2/combined_files")
all_batches <- list.dirs(path = file.path(stem, "profilers/Kraken2/outputs"),
                         full.names = FALSE, recursive = FALSE)

# Load functions
source(file.path(stem, "/profilers/Kraken2/combine_kraken2_bracken_functions.R"))

for (this_batch in all_batches) {
  # Read in combined raw output ----
  sample_names <- read.table(file.path(out_path,
                                       paste0(this_batch, "_samplenames.txt")),
                             sep = " ", header = FALSE) |>
    unlist() |>
    basename() |>
    stringr::str_remove_all("_taxa.txt")
  
  all_lines <- readLines(file.path(out_path,
                                   paste0(this_batch, "_combined.txt"))) |>
    stringr::str_remove_all("#") |>
    stringr::str_remove_all("\\[") |>
    stringr::str_remove_all("\\]") |>
    paste(collapse = "\n")
  
  # Split taxa strings ----
  df_all <- read.table(text = all_lines, sep = "\t", header = TRUE) |>
    rename_with(~ sample_names, starts_with("Sample."))
  
  # Kraken ----
  get_new_counts_tabs <- function(this_sample) {
    kraken_res <- df_all |>
      dplyr::select("Classification", dplyr::all_of(this_sample)) |>
      as_tibble()
    
    # Get the Kraken table for that sample
    kraken_new_calc <- kraken_res |>
      kraken_exclusive_counts() |>
      dplyr::mutate(!!this_sample := exclusive_count)
    return(kraken_new_calc)
    }
  
  all_sample_names <- colnames(df_all)[-1]
  all_kraken_counts <- lapply(all_sample_names, get_new_counts_tabs) |>
    set_names(all_sample_names)
  taxonomy_cols <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  joined_kr_counts <- all_kraken_counts |>
    lapply(function(x) select(x, -"taxon", -"counts", -"exclusive_count")) |>
    purrr::reduce(full_join, by = taxonomy_cols) |>
    filter(if_any(where(is.numeric), ~ . != 0))
  
  # Bracken ----
  # this_sample <- colnames(df_all)[3]
  get_new_br_tab <- function(this_sample) {
    br_tab <- read_delim(file.path(stem, "profilers/Bracken/outputs",
                                   this_batch,
                                   paste0(this_sample, ".bracken")),
                         show_col_types = FALSE) |>
      select(name, taxonomy_id, new_est_reads) |>
      mutate(across(name, ~ str_remove_all(., "\\[|\\]")))
    
    # Get the bracken table for that sample
    bracken_new_calc <- add_taxonomy_bracken(all_kraken_counts[[this_sample]],
                                             br_tab) |>
      dplyr::select(-"taxonomy_id", -!!this_sample) |>
      dplyr::rename(!!this_sample := "new_est_reads")
    
    return(bracken_new_calc)
  }
  
  all_bracken_counts <- lapply(all_sample_names, get_new_br_tab) |>
    set_names(all_sample_names)
  
  joined_br_counts <- all_bracken_counts |>
    purrr::reduce(full_join, by = taxonomy_cols) |>
    filter(if_any(where(is.numeric), ~ . != 0)) |>
    mutate(across(where(is.numeric), ~ replace_na(., 0)))
  
  # Replace NA taxa with the next available higher rank ----
  replace_na_with_higher_rank <- function(df_row) {
    new_row <- df_row
    ranks <- taxonomy_cols
    for (i in seq_along(ranks)) {
      if (i == 1) next # Skip Superkingdom
      higher_rank <- ranks[i - 1]
      current_rank <- ranks[i]
      j <- 0
      while (is.na(new_row[, current_rank])) {
        #print(new_row[, "taxid"])
        if (!is.na(df_row[, higher_rank])) {
          new_row[, current_rank] <- paste0(substr(higher_rank, 1, 1), "_",
                                            df_row[, higher_rank])
          
        } else {
          j <- j + 1
          higher_rank <- ranks[i - (1 + j)]
          if ((i - (1 + j)) < 1) new_row[, current_rank] <- "Unknown"
        }
      }
    }
    return(new_row)
  }
  
  br_counts_filled <- plyr::adply(joined_br_counts, 1, replace_na_with_higher_rank)
  kr_counts_filled <- plyr::adply(joined_kr_counts, 1, replace_na_with_higher_rank)
  
  # Construct final merged table ----
  
  table_br_final <- br_counts_filled |>
    rename_with(stringr::str_to_title) |>
    mutate(Name = Species,
           ti = as.numeric(factor(Name))) |>
    relocate("Name", "ti", .after = "Species")
  
  table_kr_final <- kr_counts_filled |>
    rename_with(stringr::str_to_title) |>
    mutate(Name = Species,
           ti = as.numeric(factor(Name))) |>
    relocate("Name", "ti", .after = "Species")
  
  # Write output table ----
  write.csv(table_br_final,
            file.path(stem, "profilers/Bracken/combined_files" ,
                      paste0(this_batch, "_bracken_final_summary.csv")),
            quote = FALSE, row.names = FALSE)
  
  write.csv(table_kr_final,
            file.path(stem, "profilers/Kraken2/combined_files" , 
                      paste0(this_batch, "_kraken2_final_summary.csv")),
            quote = FALSE, row.names = FALSE)
}
