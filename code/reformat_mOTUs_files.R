library(tidyverse)


profiler <- "mOTUs"


input_files <- list.files(
  path = paste0("profilers/", profiler, "/combined_files/updated_species_names/updated/"),
  pattern = "\\.csv$",
  full.names = TRUE
)


output_dir <- paste0("profilers/", profiler, "/combined_files/transformed/")


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


create_dummy_taxonomy <- function(name) {
  tax_levels <- list(
    Domain = "d_Unknown",
    Phylum = "p_Unknown",
    Class = "c_Unknown",
    Order = "o_Unknown",
    Family = "f_Unknown",
    Genus = "g_Unknown",
    Species = "s_Unknown"
  )
  
  return(tax_levels)
}


process_file <- function(input_file) {


  data <- read.csv(input_file, stringsAsFactors = FALSE)
  
  # Process each row
  processed_data <- data %>%
    rowwise() %>%
    mutate(
      tax_dummy = list(create_dummy_taxonomy(consensus_taxonomy))
    ) %>%
    mutate(
      Domain = tax_dummy$Domain,
      Phylum = tax_dummy$Phylum,
      Class = tax_dummy$Class,
      Order = tax_dummy$Order,
      Family = tax_dummy$Family,
      Genus = tax_dummy$Genus,
      Species = consensus_taxonomy,
      Name = consensus_taxonomy
    ) %>%
    ungroup() %>%
    select(-tax_dummy, -consensus_taxonomy, -X)
  
  processed_data$ti <- 1000 + seq_len(nrow(processed_data))
  
  # Rename sample columns: capitalize first letter
  sample_cols <- grep("^healthy_stool|^diseased_", names(processed_data), value = TRUE, ignore.case = TRUE)
  new_sample_names <- sapply(sample_cols, function(x) {
    # Capitalize first letter
    paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  })
  
  # Rename columns
  for (i in seq_along(sample_cols)) {
    names(processed_data)[names(processed_data) == sample_cols[i]] <- new_sample_names[i]
  }
  
  # Get all sample columns (now capitalized)
  sample_cols_final <- grep("^Healthy_stool|^Diseased_", names(processed_data), value = TRUE)
  
  # Reorder columns to match output format
  final_cols <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Name", "ti",
                  sample_cols_final)
  
  processed_data <- processed_data[, final_cols]
  
  return(processed_data)
}




  for (input_file in input_files) {
    tryCatch({
      processed_data <- process_file(input_file)
      cat("Processing:", basename(input_file))
      
      base_name <- basename(input_file)
      output_file <- file.path(output_dir, base_name)
      
      write.csv(processed_data, output_file, row.names = FALSE, quote = TRUE)
      cat("Saved to:", output_file)
      
    }, error = function(e) {
    })
  }
