suppressPackageStartupMessages({
  library(rvest)
  library(tidyverse)
  library(future)
  library(furrr)
  library(progressr)
  library(xml2)
})

# availableCores()
plan(multisession, workers = availableCores() - 1)

all_profilers <- c("mOTUs")

for (profiler in all_profilers) {
  csv_files <- list.files(
    path      = file.path("profilers", profiler, "combined_files/"),
    pattern   = "err.*\\.csv$",
    full.names= TRUE
  )
  
  fetch_record <- function(name_from_profilers) {
    
    # Store ORIGINAL first
    original_name <- name_from_profilers
    
    # mOTUs only
    clean_name <- name_from_profilers %>%
      str_remove("\\s*\\[.*?\\]$")
    
    url <- URLencode(
      paste0("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name=", clean_name)
    )
    
    result <- tryCatch({
      Sys.sleep(0.3)
      request <- read_html(url) 
      
      current_name <- request |> 
        html_element("body") |> 
        html_element("ul") |> 
        html_element("li") |> 
        html_element("a") |> 
        html_element("strong") |> 
        html_text2()
      
      ## special cases
      if (is.na(current_name)) {
        current_name <- clean_name
      }
      
      data.frame(
        old_name = original_name,  
        current_name = current_name
      )
    }, error = function(e) {
      message("ERROR: Failed to fetch '", clean_name, "' - ", e$message)
      
      data.frame(
        old_name = original_name,   
        current_name = clean_name
      )
    })
    
    return(result)
  }
  
  
  process_file <- function(input_path) {
    message("now in ", profiler)
    message("Processing ", basename(input_path))
    
    df <- read.csv(input_path) |>
      select(consensus_taxonomy) |> 
      distinct()
    
    # df <- read.csv(input_path) |>
    #   select(Name)
    
    ## progress
    with_progress({
      p <- progressor(steps = nrow(df))
      
      result_df <- future_map_dfr(
        
        df$consensus_taxonomy, 
        # df$species,
        
        function(name) {
          result <- fetch_record(name)
          p()  
          return(result)
        }
      )
    })
    
    diff <- result_df |> 
      filter(old_name != current_name)
    
    out_dir  <- file.path("profilers", profiler, "combined_files/updated_species_names/diff")
    
    # create folder if not exist
    dirs_to_create <- c(out_dir)
    sapply(dirs_to_create, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))
    
    out_name <- paste0("diff_", basename(input_path))
    out_path <- file.path(out_dir, out_name)
    
    write.csv(diff, file = out_path, row.names = FALSE)
  }
  
  for (file in csv_files) {
    process_file(file)
  }
}

## close multisession workers
plan(sequential)

