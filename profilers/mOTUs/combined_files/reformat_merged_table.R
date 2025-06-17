# R Script to combine outputs - mOTUs
# Metagenomics benchmarking
# Aubrey Odom
# 5/21/25

# Description ----

# This script will take all mOTUs outputs,
# consolidate them, and create a formatted
# table to use in downstream analysis

# Setup ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(LegATo)
})

# Set file paths ----
stem <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark"
batch_name <- "Mock_lluch"
dat_path <- file.path(stem, "profilers/mOTUs/outputs")
out_path <- file.path(stem, "profilers/mOTUs/combined_files")

# Read in data ----

file_path <- file.path(out_path, "merged_mOTUs_table.txt")

all_colnames <- readLines(file_path, n = 3) |>
  magrittr::extract(3) |>
  strsplit("\t") |>
  magrittr::extract2(1)

motus_data <- read.table(
  file = file_path,
  comment.char = "#",  # Skip lines starting with #
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",          # Avoid interpreting quotes in taxonomy names
  fill = TRUE          
) |>
  magrittr::set_colnames(all_colnames)

# replace with mock data ----
## Low read counts for test data....
motus_data[, 4:ncol(motus_data)] <- floor(runif(nrow(motus_data), 0, 2000))

# Separate out taxonomy ----
new_taxa_cols <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
new_output <- motus_data |>
  tidyr::separate(col = "consensus_taxonomy",
                  into = new_taxa_cols,
                  sep =  "\\|") |>
  mutate(across(all_of(new_taxa_cols), ~ sub(".__+", "", .))) |>
  mutate(Name = Species, ti = 0) |>
  relocate(Name, ti, .after = "Species") |>
  dplyr::select(-c("NCBI_tax_id", "#mOTU")) |>
  rename_with(str_to_sentence)
  
# Write to a CSV ----
write.csv(new_output, file.path(out_path, paste0(batch_name, ".csv")),
          quote = FALSE, row.names = FALSE)
