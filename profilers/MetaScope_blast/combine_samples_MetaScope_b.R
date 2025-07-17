
# R Script to combine outputs - MetaScope
# Metagenomics benchmarking
# Aubrey Odom
# 4/27/25

# Setup ----
suppressPackageStartupMessages({
  library(MetaScope)
  library(tidyverse)
  library(LegATo)
})

# Set file paths ----
stem <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/profilers/MetaScope"
test_dir <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/test_data"
dat_path <- file.path(stem, "outputs")
out_path <- file.path(stem, "combined_outputs")
taxa_db <- "/restricted/projectnb/pathoscope/data/blastdb/2024_accession_taxa/accessionTaxa.sql"
all_batches <- list.dirs(path = dat_path, full.names = FALSE, recursive = FALSE)

for (this_batch in all_batches) {
  batch_path <- file.path(dat_path, this_batch)
  # Obtain files ----
  all_files <- list.files(batch_path, 
                          pattern = ".metascope_id.csv$",
                          full.names = TRUE)
  
  # Check annot file works ----
  these_filenames <- all_files |>
    basename() |>
    stringr::str_remove(".metascope_id.csv")
  tibble(Name = these_filenames,
         Batch = this_batch) |>
    write.csv(row.names = FALSE,
              file.path(test_dir, "temp_annot.csv"))
  
  # Use convert_animalcules ----
  combined_dat <- convert_animalcules(
    meta_counts = all_files,
    annot_path = file.path(test_dir, "temp_annot.csv"),
    which_annot_col = "Name",
    end_string = ".metascope_id.csv",
    qiime_biom_out = FALSE,
    accession_path = taxa_db
  )
  
  clean_dat <- combined_dat |>
    LegATo::clean_MAE()
  
  # saveRDS(clean_dat, file.path(out_path, "combined_metascope_samples.RDS"))
  
  # Format into single table ----
  mgx <- LegATo::parse_MAE_SE(clean_dat)
  
  rowdat_2 <- mgx$tax %>%
    select(superkingdom,
           phylum,
           class,
           order,
           family,
           genus,
           species) %>%
    mutate(name = species) %>%
    rename_with(str_to_title) %>%
    mutate(ti = as.numeric(factor(Name)))
  
  all_comb <- rowdat_2 %>%
    bind_cols(mgx$counts)
  
  # Write output ----
  write.csv(all_comb, file.path(out_path, paste0(this_batch, "_combined.csv")),
            quote = FALSE, row.names = FALSE)
}