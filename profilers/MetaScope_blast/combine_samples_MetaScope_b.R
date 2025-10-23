suppressPackageStartupMessages({
  library(MetaScope)
  library(tidyverse)
  library(LegATo)
})

# Set file paths ----
stem <- "~/RProjects/metagenomics_benchmarking/profilers/MetaScope_blast"
test_dir <- file.path(stem, "test_data")
temp_dir <- file.path(stem, "outputs/temp")
dat_path <- file.path(stem, "outputs")
out_path <- file.path(stem, "combined_files")
taxa_db <- "/projects/f_wj183_1/reflib/2024_accession_taxa/accessionTaxa.sql"

# Create temp_dir if it doesn't exist
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
}

all_batches <- list.dirs(path = dat_path, full.names = FALSE, recursive = FALSE)

for (this_batch in all_batches) {
  batch_path <- file.path(dat_path, this_batch)
  
  # Create temp directory for this batch
  temp_dir_folder <- file.path(temp_dir, this_batch)
  if (!dir.exists(temp_dir_folder)) {
    dir.create(temp_dir_folder, recursive = TRUE)
  }
  
  # Obtain files ----
  all_files <- list.files(batch_path, 
                          pattern = ".metascope_blast_reassigned.csv$",
                          full.names = TRUE)
  
  for (blast_output in all_files) {
    name <- basename(blast_output)
    tmp <- read.csv(blast_output)
    
    # Select columns and save
    tmp_filtered <- tmp |> select(TaxonomyID, readsEM) |>
      rename(read_count = readsEM) |>
      dplyr::mutate(TaxonomyID = ifelse(TaxonomyID == 0, NA, TaxonomyID))
      
    write.csv(tmp_filtered, 
              file = file.path(temp_dir_folder, name),
              row.names = FALSE)
  }
  
  # Obtain temp files ----
  all_temp_files <- list.files(temp_dir_folder, 
                               pattern = ".metascope_blast_reassigned.csv$",
                               full.names = TRUE)
  
  # Check annot file works ----
  these_temp_filenames <- all_temp_files |>
    basename() |>
    stringr::str_remove(".metascope_blast_reassigned.csv")
  
  tibble(Name = these_temp_filenames,
         Batch = this_batch) |>
    write.csv(row.names = FALSE,
              file = file.path(test_dir, "temp_annot.csv"))
  
  # Use convert_animalcules ----
  combined_dat <- convert_animalcules(
    meta_counts = all_temp_files,
    annot_path = file.path(test_dir, "temp_annot.csv"),
    which_annot_col = "Name",
    end_string = ".metascope_blast_reassigned.csv",
    qiime_biom_out = FALSE,
    accession_path = taxa_db
  )
  
  clean_dat <- combined_dat |>
    LegATo::clean_MAE()
  
  # Format into single table ----
  mgx <- LegATo::parse_MAE_SE(clean_dat)
  
  rowdat_2 <- mgx$tax %>%
    select(superkingdom, phylum, class, order, family, genus, species) %>%
    mutate(name = species) %>%
    rename_with(str_to_title) %>%
    mutate(ti = as.numeric(factor(Name)))
  
  all_comb <- rowdat_2 %>%
    bind_cols(mgx$counts)
  
  # Write output ----
  write.csv(all_comb, 
            file = file.path(out_path, paste0(this_batch, "_combined.csv")),
            quote = FALSE, 
            row.names = FALSE)
}
