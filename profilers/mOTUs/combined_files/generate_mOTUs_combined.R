# This script is to take mOTUs outputs and reformat them to have the same format
# with the following headers
# Superkingdom,Phylum,Class,Order,Family,Genus,Species,Name,ti,healthy_stool_10_1 ...
# For each folder

library(tidyverse)

# Helper function to read and format each motu file
read_and_format <- function(file) {
  # sample name from file name (without extension)
  ## Split and rebuild since motus have extra parameters in the file name 
  file_name <- str_extract(basename(file), ".*(?=.motus)")
  parts <- str_split(file_name, "_", simplify = TRUE)
  sample_name <- str_c(parts[1:4], collapse = "_")
  
  df <- read.table(file, sep = "\t", fill = TRUE, comment.char = "#", quote = "")
  colnames(df) <- c("consensus_taxonomy", sample_name)
  
  return(df)
}

generate_mOTUs_combined <- function(res_paths) {
  # Read tables, skip lines 3 and merge
  tables <- map(res_paths, read_and_format)
  otu_table <- reduce(tables, full_join, by = "consensus_taxonomy")
  
  # Convert consensus taxonomy (from mOTU) into updated name taxonomy
  otu_table <- otu_table |> 
    mutate(species = )
}

res_paths <- list.files("profilers/mOTUs/outputs/mess_sample_1mil_err", pattern = ".motus", 
                        full.names = TRUE)
generate_mOTUs_combined(res_paths)

