suppressPackageStartupMessages({
  library(rvest)
  library(tidyverse)
})

input <- 'profilers/Bracken/combined_files/mess_sample_100k_err_bracken_final_summary.csv'

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
  
  
  res <- data.frame(
    old_name = name_from_profilers,
    current_name = current_name
  )
  
  
}

result <- lapply(df$Name, fetch_record)
result_df <- do.call(rbind, result)
diff <- result_df %>%
  filter(old_name != current_name)

write.csv2(diff, file = "output.csv", sep = ";")
