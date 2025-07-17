library(tidyverse)
ground_truth <- read.csv("MeSS_code/ground_truth.csv")

priors_df_merged <- ground_truth |> 
  select(species, tax_abundance, sample) |>
  group_by(sample) |> 
  group_split()


mk_priors_df_from_list <- function(i) {
  sample_name <- i$sample[1]
  df <- i |> 
    select(-sample)
  write.csv(df, paste0("profilers/MetaScope_priors/priors_df/", sample_name, ".csv"),
            row.names = FALSE)
}

lapply(priors_df_merged, mk_priors_df_from_list)
