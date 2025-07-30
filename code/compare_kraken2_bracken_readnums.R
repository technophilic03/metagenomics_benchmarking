
# Script to compare read numbers from profiler outputs
# Aubrey Odom
# 5/2/25

# Setup -----

suppressPackageStartupMessages({
  library(tidyverse)
  library(R.utils)
})
stem <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/profilers"
data_dir <- "/restricted/projectnb/pathoscope/data/Mock_Lluch_2015"

# Read in data ----

## Kraken species
k_species <- read.csv(file.path(stem, "Kraken2/combined_files", 
                                "test_lluch_Species.csv"))

## Kraken genera
k_genus <- read.csv(file.path(stem, "Kraken2/combined_files", 
                                "test_lluch_Genus.csv"))

## Bracken
b_species <- read.table(file.path(stem, "Bracken/combined_files",
                                  "test_lluch.txt"), sep = "\t",
                        header = TRUE)

## MetaScope
m_species <- read.csv(file.path(stem, "MetaScope/combined_outputs", 
                              "test_lluch.csv"))

## PathoScope 2
p_species <- read.csv(file.path(stem, "PathoScope2/combined_outputs", 
                                "test_lluch.csv"))

# Get raw sample counts ----

# Function to count reads in a fastq.gz file
count_reads <- function(file) {
  count <- system(paste("zcat", file, "| grep -c '^@'"), intern = TRUE)
  return(as.numeric(count))
}

# Get list of all fastq.gz files in the directory
files <- list.files(path = data_dir, pattern = "*.fastq.gz", full.names = TRUE)

# Create a data frame with file names and read counts
results <- data.frame(
  file = base::basename(files),
  read_count = sapply(files, count_reads)
)


# Perform comparisons ----

# Compare total read counts of Kraken 2 genera vs. species
# Create a df with rows - samples and columns - which profiler alg

k_format1 <- k_species |>
  select(starts_with("Run")) |>
  colSums()
k_format2 <- tibble(sample = names(k_format1),
                    k_read_count = k_format1) |>
  mutate(sample = tolower(sample))

k_format1_g <- k_genus |>
  select(starts_with("Run")) |>
  colSums()

k_comp <- tibble(sample = names(k_format1_g) |> tolower(),
                    k_read_count = k_format1_g) |>
  left_join(k_format2, by = "sample", suffix = c(".genus", ".species"))

k_comp

# These are not the same!
# Goal: compare total read counts of Kraken 2 species vs. Bracken output

b_format1 <- b_species |>
  select(ends_with("bracken_num")) |>
  rename_with(function(x) stringr::str_remove_all(x, ".bracken_num")) |>
  colSums() |>
  as.data.frame() |>
  rownames_to_column("sample") |>
  magrittr::set_colnames(c("sample", "b_read_count.species")) |>
  mutate(sample = tolower(sample))

# MetaScope
m_format1 <- m_species |>
  select(starts_with("Run")) |>
  colSums()
m_format2 <- tibble(sample = names(m_format1),
                    m_read_count.species = m_format1) |>
  mutate(sample = tolower(sample))

# PathoScope
p_format1 <- p_species |>
  select(starts_with("Run")) |>
  colSums()
p_format2 <- tibble(sample = names(p_format1),
                    p_read_count.species = p_format1) |>
  mutate(sample = tolower(sample))

# Look at all of the information in the outputs

# Function to count reads in a fastq.gz file
count_reads <- function(file) {
  count <- system(paste("zcat", file, "| grep -c '^@'"), intern = TRUE)
  return(as.numeric(count))
}

# Get list of all fastq.gz files in the directory
files <- list.files(path = data_dir, pattern = "*.fastq.gz", full.names = TRUE)

# Create a data frame with file names and read counts
all_read_counts <- data.frame(
  file = basename(files),
  read_count = plyr::aaply(files, 1, count_reads, .progress = "text")
)

all_read_counts_fix <- all_read_counts |>
  mutate(file = stringr::str_remove(file, ".fastq.gz"),
         sample = sub("(_R[12])$", "", file),
         read_type = sub(".*_R([12])$", "R\\1", file)) |>
  select(-file) |>
  pivot_wider(names_from = read_type, values_from = read_count) |>
  mutate(sample = tolower(sample) |> stringr::str_replace_all("-", "\\."))

# Combine all
all_combined <- k_comp |>
  left_join(b_format1, by = "sample") |>
  left_join(m_format2, by = "sample") |>
  left_join(p_format2, by = "sample") |>
  left_join(all_read_counts_fix, by = "sample")

all_combined |> View()
