RUN_SCRIPT <- TRUE

if (RUN_SCRIPT) {
  suppressPackageStartupMessages(
    library(tidyverse)
  )
  
  original_files <- list.files(
    path      = "profilers/Bracken/combined_files/",
    pattern   = "\\.csv$",
    full.names= TRUE
  )
  
  ref_files <- list.files(
    path      = "profilers/Bracken/combined_files/name_update/diff/",
    pattern   = "\\.csv$",
    full.names= TRUE
  )
  
  update <- function(original_df, reference_df) {
    updated_df <- left_join(original_df, reference_df, by = join_by(Name == old_name))
    updated_df$Name <- ifelse(
      is.na(updated_df$current_name), 
      updated_df$Name, 
      updated_df$current_name
    )
    updated_df <- updated_df |> 
      select(-current_name)
    return(updated_df)
  }
  
  process_file <- function(original, reference) {
    message("Processing ", basename(original))
    
    df <- read.csv(original)
    ref_df <- read.csv(reference)
    
    result_df <- update(df, ref_df)
    
    out_dir  <- "profilers/Bracken/combined_files/name_update/updated/"
    out_name <- paste0("updated_", basename(original))
    out_path <- file.path(out_dir, out_name)
    
    write.csv(result_df, file = out_path, row.names = FALSE)
  }
  
  map2(original_files, ref_files, process_file)
}