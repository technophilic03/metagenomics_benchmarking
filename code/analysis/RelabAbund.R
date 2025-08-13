library(tidyverse)

# Read in Ground Truth
ground_truth <- read.csv("MeSS_code/ground_truth/updated/updated_ground_truth.csv") |>
  distinct(sample, species, .keep_all = TRUE)

# Read in Profiler Results
read_in_prof_res <- function(path_to_files, pipeline) {
  file_paths <- list.files(path_to_files, full.name = TRUE)
  file_names <- basename(file_paths)
  res_list <-  map2(file_paths, file_names, ~{
    
    extracted_tag <- str_split(.y, "_", simplify = TRUE)
    tag <- str_c(extracted_tag[4], extracted_tag[5], sep = "_")
    
    read.csv(.x) |>
      pivot_longer(
        cols=contains("Healthy_stool"),
        names_to="sample_name",
        values_to="read_count") |>
    mutate(pipeline = pipeline, 
           file_name = tag,
           sample_name = str_to_lower(sample_name),
           full_sample_name = paste0(sample_name, "_", tag))
  })
  return(res_list)
}

brk_res <- read_in_prof_res("profilers/Bracken/combined_files/updated_species_names/updated", "Bracken")
ctf_res <- read_in_prof_res("profilers/Centrifuge/combined_files/updated_species_names/updated", "Centrifuge")
#clrk_res <- read_in_prof_res("profilers/Clark/combined_files/updated_species_names/updated", "Clark")
krk2_res <- read_in_prof_res("profilers/Kraken2/combined_files/updated_species_names/updated", "Kraken2")
ms_res <- read_in_prof_res("profilers/MetaScope/combined_files/updated_species_names/updated", "MetaScope")
msp_res <- read_in_prof_res("profilers/MetaScope_priors/combined_files/updated_species_names/updated", "MetaScope Priors")
#mOTUs_res <- read_in_prof_res("profilers/mOTUs/combined_files/updated_species_names/updated", "mOTUs")
ps2_res <- read_in_prof_res("profilers/PathoScope2/combined_files/updated_species_names/updated", "PathoScope2")

generate_ground_truth_categories <- function(df) {
  # Step 1: Calculate relative read counts grouped by full_sample_name
  df_relab <- df |>
    group_by(full_sample_name) |>
    mutate(rel_read_count = read_count / sum(read_count, na.rm = TRUE)) |>
    ungroup()
  
  # Step 2: Join with ground truth using sample_name and species
  res <- df_relab |>
    left_join(
      ground_truth |> mutate(is_ground_truth = species),
      by = join_by(sample_name == sample, Name == species)
    ) |>
    mutate(is_ground_truth = coalesce(is_ground_truth, "Incorrect Call"))
  
  return(res)
}

  
summary_df <- lapply(c(brk_res,ctf_res,krk2_res, ms_res, msp_res,ps2_res), generate_ground_truth_categories) |> 
  bind_rows()

replicated_ground_truth <- map_dfr(c("100k_err", "100k_noerr", "10mil_err",
                                     "10mil_noerr", "1mil_err", "1mill_noerr"),
                                   ~ground_truth |> mutate(rep = .x)) |>
  dplyr::mutate(rel_read_count = tax_abundance / 100, 
                full_sample_name = paste0(sample, "_", rep), 
                is_ground_truth = species,
                pipeline = "Ground Truth")

summary_df_ground_truth <- bind_rows(summary_df, replicated_ground_truth)

species_levels <- setdiff(unique(summary_df_ground_truth$is_ground_truth), "Incorrect Call")
summary_df_ground_truth$is_ground_truth <- factor(summary_df_ground_truth$is_ground_truth, 
                             levels = c("Incorrect Call", species_levels))
summary_df_ground_truth$pipeline <- factor(summary_df_ground_truth$pipeline,
                                           levels = c("Ground Truth",
                                                      "Centrifuge",
                                                      "Kraken2",
                                                      "Bracken",
                                                      "PathoScope",
                                                      "MetaScope",
                                                      "MetaScope Priors"))


p1 <- ggplot(data = summary_df_ground_truth, 
             aes(fill = is_ground_truth, y = rel_read_count, x = full_sample_name)) +
  geom_bar(position ="stack", stat = "identity")+
  #scale_fill_manual(values = wheel_colors, name = "Species") + 
  ylab("Relative Abundance") + 
  xlab("") +
  #scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
  facet_grid(~pipeline) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

p1

         