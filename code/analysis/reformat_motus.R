library(tidyverse)

stem = "profilers/mOTUs/combined_files/updated_species_names/transformed"
motus_paths <- list.files(stem)
 
convert_counts <- function(path) {
  res <- read.csv(file.path(stem, path))
   if (grepl("100k",  path, fixed = TRUE)) {
     total_reads <- 1e5
   } else if (grepl("1mil",  path, fixed = TRUE)) {
     total_reads <- 1e6
   } else if (grepl("10mil", path, fixed = TRUE)) {
     total_reads <- 1e7
   } else stop("Could not infer total reads from `path` (looked for '100k', '1mil', '10mil').")
   
   num_after10 <- names(res)[sapply(res, is.numeric) & seq_along(res) > 9]
   
   output <- res |>
     dplyr::filter(dplyr::if_any(dplyr::all_of(num_after10), ~ . != 0)) |>
     dplyr::mutate(dplyr::across(dplyr::all_of(num_after10), ~ . * total_reads))
   
   
   write.csv(output, 
             file.path("profilers/mOTUs/combined_files/updated_species_names/transformed_counts", path), 
             row.names = FALSE)
   return(output)
 }
lapply(motus_paths, convert_counts)

