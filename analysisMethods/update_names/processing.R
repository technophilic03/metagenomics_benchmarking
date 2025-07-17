suppressPackageStartupMessages({
  library(rvest)
  library(tidyverse)
})

csv_files <- list.files(
  path      = "profilers/Bracken/combined_files",
  pattern   = "\\.csv$",       # only .csv
  full.names= TRUE             # give full path
)

df <- read.csv(input) |>
  select(Name)


fetch_record <- function(name_from_profilers) {
  url <- URLencode(
    paste0("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name=", name_from_profilers)
  )
  
  request <- read_html(url)
  current_name <- request |> 
    html_element("body") |> 
    html_element("form") |> 
    html_element("ul") |> 
    html_element("strong") |> 
    html_text2()
  
  if (!is.na(current_name) && current_name == "1)") {
    current_name <- name_from_profilers
  }else if (grepl("\\([A-Za-z]+ [0-9]+\\)", current_name)) {
    bracketed <- regmatches(current_name,
                            regexpr("\\[[A-Za-z]+\\]", current_name))
    if (length(bracketed) > 0 && bracketed != name_from_profilers) {
      current_name <- bracketed
    }else{
      current_name <- name_from_profilers
    }
  }
  
  
  res <- data.frame(
    old_name = name_from_profilers,
    current_name = current_name
  )
}

process_file <- function(input_path) {
  message("Processing ", basename(input_path))
  
  df <- read.csv(input_path) |>
    select(Name)
  
  result <- lapply(df$Name, fetch_record)
  result_df <- do.call(rbind, result)
  diff <- result_df |> 
    filter(old_name != current_name)
  
  out_dir  <- "profilers/Bracken/combined_files/name_update/"
  out_name <- paste0("diff_", basename(input_path))
  out_path <- file.path(out_dir, out_name)
  
  write.csv(diff, file = out_path, row.names = FALSE)
}
# purrr::walk(csv_files, process_file)