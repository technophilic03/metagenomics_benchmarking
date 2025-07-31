RUN_SCRIPT <- TRUE

if (RUN_SCRIPT) {
  suppressPackageStartupMessages({
    library(rvest)
    library(tidyverse)
    library(future)
    library(furrr)
    library(progressr)
  })
  
  # availableCores()
  plan(multisession, workers = 4)
  
  csv_files <- list.files(
    path      = "profilers/MetaScope_priors/combined_files/",
    pattern   = "\\.csv$",
    full.names= TRUE
  )
  
  fetch_record <- function(name_from_profilers) {
    url <- URLencode(
      paste0("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name=", name_from_profilers)
    )
    
    Sys.sleep(0.1)
    
    request <- read_html(url)
    current_name <- request |> 
      html_element("body") |> 
      html_element("form") |> 
      html_element("ul") |> 
      html_element("strong") |> 
      html_text2()
    
    ## special cases
    if (!is.na(current_name) && current_name == "1)") {
      current_name <- name_from_profilers
    }
    
    data.frame(
      old_name = name_from_profilers,
      current_name = current_name
    )
  }
  
  process_file <- function(input_path) {
    message("Processing ", basename(input_path))
    
    df <- read.csv(input_path) |>
      select(Name) |> 
      distinct()
    
    # df <- read.csv(input_path) |>
    #   select(Name)
    
    ## progress
    with_progress({
      p <- progressor(steps = nrow(df))
      
      result_df <- future_map_dfr(
        
        df$Name, 
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
    
    out_dir  <- "profilers/MetaScope_priors/combined_files/updated_species_names/"
    
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
  
  ## close multisession workers
  plan(sequential)
}
