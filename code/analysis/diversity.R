library(vegan)
library(tidyverse)


summary_df <- read.csv("code/analysis/summary_df.csv")

summary_df_filtered <- summary_df |>
  select(Name, read_count, pipeline, full_sample_name, rel_read_count, species) |>
  mutate(read_count = case_when(
    pipeline %in% c("mOTUs", "Ground Truth") & str_detect(full_sample_name, "100k") ~ round(rel_read_count * 100000),
    pipeline %in% c("mOTUs", "Ground Truth") & str_detect(full_sample_name, "1mil") ~ round(rel_read_count * 1000000),
    pipeline  %in% c("mOTUs", "Ground Truth") & str_detect(full_sample_name, "10mil") ~ round(rel_read_count * 10000000),
    TRUE ~ read_count  # keep original value if no condition matches)
  )) |>
  mutate(Name = case_when(
    pipeline == "Ground Truth" ~ species,
    TRUE ~ Name
  )) |>
  select(-rel_read_count, -species)|>
  filter(Name != "") |>
  filter(Name != "Unassigned") |>
  filter(read_count > 0) |>
  group_by(Name, pipeline, full_sample_name) |> 
  summarise(read_count = sum(read_count), .groups = "drop") |>
  filter(pipeline != "MetaScope Priors") |>
  filter(pipeline != "MetaScope Blast")




calculate_diversity <- function(df, index) {
  df <- df |>
    filter(Name != "") |>
    filter(Name != "Unassigned") |>
    select(-pipeline) |> 
    pivot_wider(names_from = Name, values_from = read_count, values_fn = sum, values_fill = 0) |>
    column_to_rownames("full_sample_name")
  vegan::diversity(df, index = index)
}



# Calculate Diversity and Plot Diversity + RMSE
alpha_diversity_plots <- function(index, df) {
  
  # Use df passed into the function
  split_list <- df |>
    group_by(pipeline) |>
    group_split()
  
  group_names <- df |>
    group_by(pipeline) |>
    group_keys() |>
    pull(pipeline)
  
  names(split_list) <- group_names
  

  shannon_list <- lapply(split_list, function(x) calculate_diversity(x, index = index))
  
  shannon_diff <- lapply(shannon_list, function(v) v - shannon_list[[3]][names(v)])
  
  shannon_diff_df <- enframe(shannon_diff, name = "pipeline", value = "differences") |>
    unnest_longer(differences, values_to = "difference", indices_to = "difference_id") |>
    mutate(total_species = str_split_i(difference_id, "_", 3), 
           species_distribution = str_split_i(difference_id, "_", 4), 
           total_reads = str_split_i(difference_id, "_", 5), 
           err_model = str_split_i(difference_id, "_", 6)) |>
    filter(pipeline != "Ground Truth")
  shannon_diff_df$pipeline <- factor(shannon_diff_df$pipeline, 
                                     levels = c("Centrifuge",
                                                "Kraken2",
                                                "Bracken",
                                                "mOTUs",
                                                "PathoScope2",
                                                "MetaScope"))
  
  p1 <- ggplot(shannon_diff_df, aes(x = pipeline, y = difference, fill = pipeline)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter(aes(shape = total_species, color = pipeline), alpha = 0.5) +
    theme_minimal() +
    labs(y = "Difference from ground truth", x = "Pipeline",
         title = paste(index, "— Distribution of deviations from ground truth")) + 
    coord_flip() +
    geom_hline(yintercept = 0)
  
  p2 <- ggplot(shannon_diff_df, aes(x = difference_id, y = difference, color = pipeline, group = pipeline)) +
    geom_line(alpha = 0.6) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)) +
    labs(y = "Difference", x = "Sample", title = paste(index, "— Per-sample deviations by pipeline"))
  
  shannon_rmse <- shannon_diff_df |>
    group_by(pipeline) |>
    summarise(RMSE = sqrt(mean(difference^2)))
  
  p3 <- ggplot(shannon_rmse, aes(x = pipeline, y = RMSE, fill = pipeline)) +
    geom_col() +
    theme_minimal() +
    labs(title = paste(index, "— RMSE of pipelines"), y = "RMSE", x = "Pipeline")
  
  list(p1, p2, p3)
}



## By definition MetaScope blast decreases alpha diversity
alpha_plots <- lapply(c("shannon", "simpson", "inv"), alpha_diversity_plots, df = summary_df_filtered)


alpha_plots[[1]][1]
alpha_plots[[1]][2]
alpha_plots[[1]][3]

alpha_plots[[2]][1]
alpha_plots[[2]][2]
alpha_plots[[2]][3]

alpha_plots[[3]][1]
alpha_plots[[3]][2]
alpha_plots[[3]][3]


################################################################################
shannon_list <- lapply(split_list, function(x) calculate_diversity(x, index = "shannon"))

shannon_diff <- lapply(shannon_list, function(v) v - shannon_list[[3]][names(v)])

shannon_diff_df <- enframe(shannon_diff, name = "pipeline", value = "differences") |>
  unnest_longer(differences, values_to = "difference", indices_to = "difference_id") |>
  mutate(total_species = str_split_i(difference_id, "_", 3), 
         species_distribution = str_split_i(difference_id, "_", 4), 
         total_reads = str_split_i(difference_id, "_", 5), 
         err_model = str_split_i(difference_id, "_", 6)) |>
  filter(pipeline != "Ground Truth") 
shannon_diff_df$pipeline <- factor(shannon_diff_df$pipeline, 
                                   levels = c("Centrifuge",
                                              "Kraken2",
                                              "Bracken",
                                              "mOTUs",
                                              "PathoScope2",
                                              "MetaScope"))

simpson_list <- lapply(split_list, function(x) calculate_diversity(x, index = "simpson"))

simpson_diff <- lapply(simpson_list, function(v) v - simpson_list[[3]][names(v)])

simpson_diff_df <- enframe(simpson_diff, name = "pipeline", value = "differences") |>
  unnest_longer(differences, values_to = "difference", indices_to = "difference_id") |>
  mutate(total_species = str_split_i(difference_id, "_", 3), 
         species_distribution = str_split_i(difference_id, "_", 4), 
         total_reads = str_split_i(difference_id, "_", 5), 
         err_model = str_split_i(difference_id, "_", 6)) |>
  filter(pipeline != "Ground Truth")
simpson_diff_df$pipeline <- factor(simpson_diff_df$pipeline, 
                                   levels = c("Centrifuge",
                                              "Kraken2",
                                              "Bracken",
                                              "mOTUs",
                                              "PathoScope2",
                                              "MetaScope"))

inv_list <- lapply(split_list, function(x) calculate_diversity(x, index = "inv"))

inv_diff <- lapply(inv_list, function(v) v - inv_list[[3]][names(v)])

inv_diff_df <- enframe(inv_diff, name = "pipeline", value = "differences") |>
  unnest_longer(differences, values_to = "difference", indices_to = "difference_id") |>
  mutate(total_species = str_split_i(difference_id, "_", 3), 
         species_distribution = str_split_i(difference_id, "_", 4), 
         total_reads = str_split_i(difference_id, "_", 5), 
         err_model = str_split_i(difference_id, "_", 6)) |>
  filter(pipeline != "Ground Truth")
inv_diff_df$pipeline <- factor(inv_diff_df$pipeline, 
                                   levels = c("Centrifuge",
                                              "Kraken2",
                                              "Bracken",
                                              "mOTUs",
                                              "PathoScope2",
                                              "MetaScope"))

shannon_diff_df <- shannon_diff_df |> dplyr::mutate(alpha_metric = "Shannon")
simpson_diff_df <- simpson_diff_df |> dplyr::mutate(alpha_metric = "Simpson")
inv_diff_df <- inv_diff_df |> dplyr::mutate(alpha_metric = "Inverse simpson")

full_alpha_df <- rbind(shannon_diff_df, simpson_diff_df, inv_diff_df)
full_alpha_df$alpha_metric <- factor(full_alpha_df$alpha_metric, 
                                     levels = c("Shannon", "Simpson", "Inverse simpson"))

alpha_plot_full <- ggplot(full_alpha_df, aes(x = pipeline, y = difference, fill = pipeline)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(aes(shape = total_species, color = pipeline), alpha = 0.5) +
  facet_wrap(vars(alpha_metric), scales = "free_x") +
  labs(y = "Difference from ground truth", x = "Pipeline",
       title = paste(index, "— Distribution of deviations from ground truth")) + 
  coord_flip() +
  geom_hline(yintercept = 0)

ggsave("results/figures/alpha_div_plot.png", alpha_plot_full, dpi = 600)

################################################################################
# Beta Diversity
################################################################################
comm_matrix <- summary_df_filtered |>
  mutate(full_sample_name = paste0(pipeline, "_", full_sample_name)) |> 
  select(full_sample_name, Name, read_count) |>
  pivot_wider(names_from = Name, values_from = read_count, values_fill = list(read_count = 0)) |>
  column_to_rownames("full_sample_name")

beta <- betadiver(comm_matrix, method = "sor")

# Suppose you have metadata linking samples to pipelines
metadata <- data.frame(full_sample_name = row.names(comm_matrix)) |>
  separate_wider_delim(full_sample_name, delim = "_",cols_remove = FALSE,
                       names = c("pipeline", "healthy", "stool", "total_species", "distribution_score", "seq_depth", "err_model"))
  

# Run PERMANOVA
adonis2(beta ~ pipeline+total_species+distribution_score+seq_depth+err_model, data = metadata)

# PCoA
# Compute beta diversity (e.g., Bray-Curtis or Sørensen)
dist_mat <- vegdist(comm_matrix, method = "bray")  # or use betadiver()

# Run PCoA
pcoa_res <- cmdscale(dist_mat, eig = TRUE, k = 2)  # 2D solution

# Put into a data frame for plotting
pcoa_df <- data.frame(
  Sample = rownames(comm_matrix),
  Axis1 = pcoa_res$points[,1],
  Axis2 = pcoa_res$points[,2]
)


pcoa_df <- pcoa_df %>%
  left_join(metadata, by = c("Sample" = "full_sample_name"))

# Plot
ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = pipeline, label = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCoA (Bray-Curtis)", x = "PCoA1", y = "PCoA2") 
  #facet_grid(vars(seq_depth, distribution_score))
  #ggforce::geom_mark_ellipse(aes(fill = seq_depth), show.legend = FALSE, alpha = 0.1)
##################
#NMDS

# Run NMDS (k = 2 dimensions)
nmds_res <- metaMDS(comm_matrix, distance = "bray", k = 2, trymax = 100)

# Extract scores
# Get only the sample coordinates
nmds_sites <- scores(nmds_res, display = "sites")
nmds_df <- as.data.frame(nmds_sites)

# Add sample names
nmds_df$Sample <- rownames(nmds_sites)

# Join with metadata
nmds_df <- nmds_df %>%
  dplyr::left_join(metadata, by = c("Sample" = "full_sample_name"))

# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = pipeline, label = Sample, shape = distribution_score)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS (Bray-Curtis)")


# Self calculate bray distance fromg round truth
calc_bray <- function(x,y) {
  sum(abs(x-y)) / sum(x+y)
}

bray_dist_df <- metadata |> 
  filter(pipeline != "Ground Truth") |>
  mutate(corres_ground_truth = paste0("Ground Truth_healthy_stool_", total_species, "_", distribution_score, "_", seq_depth, "_", err_model)) |>
  rowwise() |> 
  mutate(bray_dist = calc_bray(comm_matrix[full_sample_name,], comm_matrix[corres_ground_truth,])) |>
  ungroup()

bray_dist_df$pipeline <- factor(bray_dist_df$pipeline, levels = c("Centrifuge",
                                                                  "Kraken2",
                                                                  "Bracken",
                                                                  "mOTUs",
                                                                  "PathoScope2",
                                                                  "MetaScope"))
bray_dist_df$distribution_score <- factor(bray_dist_df$distribution_score, levels = c("1", "10", "50", "100"))
bray_dist_df$seq_depth <- factor(bray_dist_df$seq_depth, levels = c("100k", "1mil", "10mil"))
                              

final_beta_plot <- ggplot(bray_dist_df, aes(x = pipeline, y = bray_dist,
                         color = pipeline)) + 
  geom_boxplot() + 
  facet_grid(distribution_score  ~ seq_depth) +
  ylab("Bray-Curtis Dissimilarity") +
  xlab("Profiler") +
  labs(title = paste("Boxplot of Bray-Curtis Dissimilarity from Ground Truth")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(sec.axis = dup_axis(name = "Species Dominance Score",
                                         breaks = NULL, labels = NULL)) + 
  scale_x_discrete(sec.axis = dup_axis(name = "Total Simulated Reads",
                                         breaks = NULL, labels = NULL))

ggsave("results/figures/final_beta_plot.png", final_beta_plot, dpi = 600)
