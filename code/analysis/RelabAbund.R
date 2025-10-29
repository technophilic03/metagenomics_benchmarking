library(tidyverse)
library(paletteer)

# Read in Ground Truth
ground_truth <- read.csv("MeSS_code/ground_truth/updated/updated_ground_truth.csv") |>
  distinct(sample, species, .keep_all = TRUE) |>
  filter(!grepl("healthy_stool_500", sample))

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
           sample_name = str_to_lower(str_replace(sample_name, "_noerr", "")),
           full_sample_name = paste0(sample_name, "_", tag)) |>
    rename(any_of(c(Name = "consensus_taxonomy")))
  })
  return(res_list)
}

brk_res <- read_in_prof_res("profilers/Bracken/combined_files/updated_species_names/updated", "Bracken")
ctf_res <- read_in_prof_res("profilers/Centrifuge/combined_files/updated_species_names/updated", "Centrifuge")
krk2_res <- read_in_prof_res("profilers/Kraken2/combined_files/updated_species_names/updated", "Kraken2")
ms_res <- read_in_prof_res("profilers/MetaScope/combined_files/updated_species_names/updated", "MetaScope")
msp_res <- read_in_prof_res("profilers/MetaScope_priors/combined_files/updated_species_names/updated", "MetaScope Priors")
msb_res <- read_in_prof_res("profilers/MetaScope_blast/combined_files/updated_species_names/updated", "MetaScope Blast")
mOTUs_res <- read_in_prof_res("profilers/mOTUs/combined_files/updated_species_names/updated", "mOTUs")
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

profiler_res <- c(ms_res, msp_res,msb_res, ps2_res, ctf_res, krk2_res, brk_res, mOTUs_res)
summary_df <- lapply(profiler_res, generate_ground_truth_categories) |> 
  bind_rows() |>
  filter(rel_read_count > 0)

replicated_ground_truth <- map_dfr(c("100k_err", "100k_noerr", "10mil_err",
                                     "10mil_noerr", "1mil_err", "1mil_noerr"),
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
                                                      "mOTUs",
                                                      "Kraken2",
                                                      "Bracken",
                                                      "PathoScope2",
                                                      "MetaScope",
                                                      "MetaScope Priors",
                                                      "MetaScope Blast"
                                                      ))

top_species <- ground_truth |>
  dplyr::filter(sample == "healthy_stool_100_1") |> 
  dplyr::arrange(desc(tax_abundance)) |> 
  dplyr::pull(species)
top_10_species <- c("Incorrect Call", top_species[1:10])
top_20_species <- c("Incorrect Call", top_species[1:20])

wheel_colors = paletteer::paletteer_d("khroma::soil", 21)
wheel_colors <- append(wheel_colors, c("grey85"), after = 0)


summary_df_ground_truth_10 <- summary_df_ground_truth |>
  dplyr::mutate(
    is_ground_truth = as.character(is_ground_truth),
    is_ground_truth = ifelse(is_ground_truth %in% top_10_species, is_ground_truth, "other")
  )
summary_df_ground_truth_10$is_ground_truth <- factor(
  summary_df_ground_truth_10$is_ground_truth, levels = c(top_10_species, "other")
)

summary_df_ground_truth_20 <- summary_df_ground_truth |>
  dplyr::mutate(
    is_ground_truth = as.character(is_ground_truth),
    is_ground_truth = ifelse(is_ground_truth %in% top_20_species, is_ground_truth, "other")
  )
summary_df_ground_truth_20$is_ground_truth <- factor(
  summary_df_ground_truth_20$is_ground_truth, levels = c(top_20_species, "other")
)

p1 <- ggplot(data = summary_df_ground_truth_20, 
             aes(fill = is_ground_truth, y = rel_read_count, x = full_sample_name)) +
  geom_bar(position ="stack", stat = "identity")+
  scale_fill_manual(values = wheel_colors, name = "Species") + 
  ylab("Relative Abundance") + 
  xlab("") +
  #scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
  facet_grid(cols = vars(pipeline)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom", 
        text = element_text(size = 8, family = "Arial"))

p1

ggsave("results/figures/relab_plot.png", p1, dpi = 600, units = "mm", width = 240, height = 120)



write.csv(summary_df_ground_truth, "code/analysis/summary_df.csv", row.names = FALSE)
         