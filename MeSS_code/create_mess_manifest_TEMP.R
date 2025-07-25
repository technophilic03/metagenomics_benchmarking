
# Temporary script
# Identify which samples did not run
# Create manifest

# Aubrey Odom

suppressPackageStartupMessages({
  library(tidyverse)
})

stem <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/profilers/MetaScope/outputs"
all_files <- list.files(stem, recursive = TRUE, pattern = "id.csv$")

# Look at processed samples

all_samps_p <- stringr::str_split(all_files, pattern = "/") |>
  unlist() |>
  matrix(ncol = 2, byrow = TRUE,
         dimnames = list(c(), c("folder", "sample"))) |>
  as.data.frame() |>
  dplyr::mutate(sample = stringr::str_remove(sample, ".metascope_id.csv"),
                "completed" = 1)


# look at Mess files

mess_stem <- "~/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/MeSS_code"
mess_manifest <- read_table(file.path(mess_stem, "mess_manifest.txt"),
                            col_names = FALSE) |>
  mutate(sample = stringr::str_split_i(X1, "out_samples/", i = 2),
         sample = stringr::str_remove(sample, "_R1.fq.gz")) |>
  separate(sample, into = c("folder", NA, "sample"), sep = "/")

# Join to id missing samples
missing_samples <- mess_manifest |>
  left_join(all_samps_p, by = c("folder", "sample")) |>
  filter(is.na(completed)) |>
  select(X1)

# Export
write.table(missing_samples, file.path(mess_stem, "mess_manifest_TEMP.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)



