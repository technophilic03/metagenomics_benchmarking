# R Script to combine outputs - Centrifuge
# Metagenomics benchmarking
# Aubrey Odom
# 6/9/25

# Setup ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(taxize)
  library(plyr)
  library(rlang)
  library(purrr)
})

stem <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/profilers/Centrifuge"
Sys.setenv(ENTREZ_KEY = "01d22876be34df5c28f4aedc479a2674c809")
num_tries <- 5
all_batches <- list.dirs(path = file.path(stem, "outputs"), full.names = FALSE, recursive = FALSE)

for (this_batch in all_batches) {
  # Read input files ----
  end_string <- "_summary.tsv"
  sample_paths <- list.files(file.path(stem, "outputs", this_batch), pattern = end_string)
  sample_basenames <- stringr::str_remove(sample_paths, end_string)
  sample_summaries <- lapply(file.path(stem, "outputs", this_batch, sample_paths),
                             function(x) readr::read_tsv(x, show_col_types = FALSE) |>
                               select("name", "taxID", "numReads")) |>
    set_names(sample_basenames) |>
    # Rename numReads column to the sample name (given as list object name)
    imap(~ dplyr::rename(.x, !!sym(.y) := numReads))
  
  # Merge all summaries ----
  merged_samples <- reduce(sample_summaries, full_join,
                           by = c("name", "taxID")) |>
    # Replace NA with 0 for read counts
    mutate(across(where(is.numeric), ~ replace_na(., 0))) |>
    arrange(name)
  
  # Match Taxonomy IDs ----
  
  # Run function
  taxon_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species", "strain")
  sub_tax_ranks <- c("superkingdom",
                     "phylum",
                     "class",
                     "order",
                     "family",
                     "genus",
                     "species")
  class_taxon <- function(taxon, num_tries) {
    na_table <- data.frame(name = "Unknown", rank = taxon_ranks, id = 0)
    if (is.na(taxon)) return(na_table)
    success <- FALSE
    attempt <- 0
    e <- "NCBI request not granted. Re-attempting request."
    while (!success ) {
      try({
        attempt <- attempt + 1
        if (attempt <= num_tries) {
          tryCatch({
            classification_table <- taxize::classification(
              taxon, db = "ncbi", max_tries = num_tries)[[1]]},
            error = function(w) stop(e)
          )
        }
        if (attempt > num_tries) {
          message("UID ", taxon, " not found. Continuing search for next UID.")
          return(na_table)
        }
        success <- TRUE
      })
    }
    return(classification_table)
  }
  
  all_ncbi <- plyr::llply(merged_samples$taxID, .fun = class_taxon,
                          .progress = "text",
                          num_tries = num_tries)
  
  taxonomy_table <- plyr::llply(all_ncbi, MetaScope:::mk_table, taxon_ranks) |>
    dplyr::bind_rows() |> as.data.frame() |>
    magrittr::set_colnames(taxon_ranks) |>
    select(all_of(sub_tax_ranks))
  taxonomy_table$taxid <- merged_samples$taxID
  # Remove any brackets
  taxonomy_table$species <- gsub("\\[|\\]", "", taxonomy_table$species)
  
  # Replace NA taxa with the next available higher rank ----
  replace_na_with_higher_rank <- function(df_row) {
    new_row <- df_row
    ranks <- sub_tax_ranks
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
  
  # Replace NA values with the next available higher rank
  tax_tab_replaced <- plyr::adply(taxonomy_table, 1, replace_na_with_higher_rank)
  
  # Construct final merged table ----
  
  table_final <- tax_tab_replaced |>
    dplyr::left_join(dplyr::rename(merged_samples, "taxid" = "taxID"),
                     by = "taxid") |>
    select(-"taxid") |>
    rename_with(stringr::str_to_title) |>
    mutate(ti = as.numeric(factor(Name))) |>
    relocate("ti", .after = "Name")
  
  # Write output table ----
  write.csv(table_final,
            file.path(stem, "combined_outputs" ,
                      paste0(this_batch, "centrifuge_summary.csv")),
            quote = FALSE, row.names = FALSE)
}