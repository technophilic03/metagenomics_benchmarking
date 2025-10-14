library(tidyverse)

paths <- list.files(path="profilers/mOTUs/outputs", pattern = ".motus", full.names = TRUE)

read_and_format <- function(file) {
  # sample name from file name (without extension)
  sample_name <- str_extract(basename(file), ".*(?=.motus)")
  
  df <- read.table(file, sep = "\t", fill = TRUE, comment.char = "#", quote = "")
  colnames(df) <- c("consensus_taxonomy", sample_name)
  
  return(df)
}

tables <- map(paths, read_and_format)
# Full join all by taxid
otu_table <- reduce(tables, full_join, by = "consensus_taxonomy")
