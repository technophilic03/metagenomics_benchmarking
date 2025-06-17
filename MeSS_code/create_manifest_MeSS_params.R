
# Script to create sample combinations
# Aubrey Odom
# 5/1/25

library(tidyverse)
stem <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/MeSS_code"
options(scipen = 999)

# For notation
convert_to_readable <- function(x) {
  case_when(
    x >= 1e9 ~ paste0(round(x / 1e9, 1), "bil"),
    x >= 1e6 ~ paste0(round(x / 1e6, 1), "mil"),
    x >= 1e3 ~ paste0(round(x / 1e3, 1), "k"),
    TRUE ~ as.character(x)
  )
}


# First, define params
err <- rep(c(TRUE, FALSE), times = 3)
start_base <- 3 * 10000000 
bases <- start_base * c(1, 10, 100) |>
  rep(each = 2)

# Create matrix
res_table <- tibble(errors = err, total_bases = bases) |>
  mutate(err_label = ifelse(errors, yes = "err", no = "noerr"),
         bases_readable = convert_to_readable(total_bases/(150*2)),
         name = paste("mess_sample",
                      bases_readable, err_label,
                      sep = "_")) |>
  select(-err_label, -bases_readable) |>
  relocate(name)

# output matrix
write.table(res_table, file.path(stem, "manifest_MeSS_params.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = FALSE)
