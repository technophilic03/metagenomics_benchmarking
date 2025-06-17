
# R Script to combine outputs - PathoScope 2
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
batch_name <- "test_lluch"
stem <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark"
test_dir <- file.path(stem, "test_data")
dat_path <- file.path(stem, "profilers/PathoScope2/outputs")
out_path <- file.path(stem, "profilers/PathoScope2/combined_outputs")
taxa_db <- "/restricted/projectnb/pathoscope/data/blastdb/2024_accession_taxa/accessionTaxa.sql"

# Obtain files ----
ending <- "-sam-report.tsv"
all_files <- list.files(dat_path, 
                        pattern = paste0(ending, "$"),
                        full.names = TRUE)

# Check annot file works ----
these_filenames <- all_files |>
  basename() |>
  stringr::str_remove(ending)
read.csv(file.path(test_dir, "annot.csv")) |>
  filter(Name %in% these_filenames) |>
  write.csv(row.names = FALSE, quote = FALSE,
              file.path(test_dir, "temp_annot.csv"))

# Use convert_animalcules ----
combined_dat <- convert_animalcules_patho(
  patho_counts = dat_path,
  annot_path = file.path(test_dir, "temp_annot.csv"),
  which_annot_col = "Name",
  end_string = ending
)

clean_dat <- combined_dat |>
  LegATo::clean_MAE()

saveRDS(clean_dat, file.path(out_path, "combined_patho_samples.RDS"))

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
write.csv(all_comb, file.path(out_path, paste0(batch_name, ".csv")),
          quote = FALSE, row.names = FALSE)
