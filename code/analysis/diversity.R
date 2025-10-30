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
  summarise(read_count = sum(read_count), .groups = "drop")




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
    unnest_longer(differences, values_to = "difference", indices_to = "difference_id")
  
  p1 <- ggplot(shannon_diff_df, aes(x = pipeline, y = difference, fill = pipeline)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter() +
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
adonis2(beta ~ pipeline, data = metadata)

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
  labs(title = "PCoA (Bray-Curtis)", x = "PCoA1", y = "PCoA2") +
  #facet_grid(vars(seq_depth, distribution_score)) + 
  ggforce::geom_mark_ellipse(aes(fill = seq_depth), show.legend = FALSE)
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

