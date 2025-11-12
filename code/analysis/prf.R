# Setup
library(tidyverse)
library("pracma")
summary_df <- read.csv("code/analysis/summary_df.csv")

summary_df_filtered <- summary_df |>
  select(Name, read_count, pipeline, full_sample_name, rel_read_count, species, is_ground_truth) |>
  mutate(read_count = case_when(
    pipeline %in% c("mOTUs", "Ground Truth") & str_detect(full_sample_name, "100k") ~ round(rel_read_count * 100000),
    pipeline %in% c("mOTUs", "Ground Truth") & str_detect(full_sample_name, "1mil") ~ round(rel_read_count * 1000000),
    pipeline  %in% c("mOTUs", "Ground Truth") & str_detect(full_sample_name, "10mil") ~ round(rel_read_count * 10000000),
    TRUE ~ read_count  # keep original value if no condition matches)
  )) |>
  mutate(Name = case_when(
    pipeline == "Ground Truth" ~ species,
    TRUE ~ Name
  )) 

ground_truth_df <- summary_df_filtered |>
  filter(pipeline == "Ground Truth") |>
  rename("gt_rel_read_count" = rel_read_count) |>
  select(Name, full_sample_name, gt_rel_read_count)


# Helper Function to calculate precision, recall, and F1 and species level 
#  relative abundance threshold
calc_prf <- function(abundance_threshold) {
  prf_df <-  summary_df_filtered |> 
    group_by(full_sample_name , Name, pipeline, is_ground_truth) |>
    summarise(value = sum(rel_read_count), .groups = "drop") |>  # Merge same sample and taxids together
    filter(value > abundance_threshold | is_ground_truth == "Incorrect Call") |> # always keep is ground trouth
    left_join(ground_truth_df, by = c("full_sample_name", "Name"), ) |>
    replace_na(list(gt_rel_read_count = 0)) |>
    mutate(abund_weight = 1 - abs(value - gt_rel_read_count), 
           is_tp = ifelse(is_ground_truth == "Incorrect Call", 0, 1)) |> # weight by distance from ground truth abundance
    group_by(full_sample_name, pipeline) |> 
    summarise(
      TP = sum(is_tp * abund_weight),
      precision = TP / n(),
      recall = TP / ifelse(strsplit(full_sample_name, split = "_")[[1]][3] == "100", 100, 10),
      F1 = 2 * precision * recall / (precision + recall),
      .groups = "drop")
  prf_df$abundance_threshold <- abundance_threshold
  return(prf_df)
}

calc_prf <- function(abundance_threshold) {
  prf_df <-  summary_df_filtered |> 
    group_by(full_sample_name , Name, pipeline, is_ground_truth) |>
    summarise(value = sum(rel_read_count), .groups = "drop") |>  # Merge same sample and taxids together
    filter(value > abundance_threshold | is_ground_truth == "Incorrect Call") |> # always keep is ground trouth
    left_join(ground_truth_df, by = c("full_sample_name", "Name"), ) |>
    replace_na(list(gt_rel_read_count = 0)) |>
    mutate(abund_weight = 1, 
           is_tp = ifelse(is_ground_truth == "Incorrect Call", 0, 1)) |> # weight by distance from ground truth abundance
    group_by(full_sample_name, pipeline) |> 
    summarise(
      TP = sum(is_tp * abund_weight),
      precision = TP / n(),
      recall = TP / ifelse(strsplit(full_sample_name, split = "_")[[1]][3] == "100", 100, 10),
      F1 = 2 * precision * recall / (precision + recall),
      .groups = "drop")
  prf_df$abundance_threshold <- abundance_threshold
  return(prf_df)
}

threshold_vals <- c(0,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-1,5e-1, 1)
threshold_vals <- seq(0,0.1, 0.001)
pr_curve_df <- map_dfr(threshold_vals,calc_prf)
pr_curve_df <- pr_curve_df |>
  mutate(across(F1, ~replace(.,is.nan(.),0)))

write.csv(pr_curve_df, "code/analysis/prf_df.csv", row.names = FALSE)
#pr_curve_df <- read.csv("code/analysis/prf_df.csv")

# Calculate AUC for precision recall curve
aauc <- pr_curve_df |>
  group_by(full_sample_name, pipeline) |>
  arrange(recall) |>
  summarize(cum_area = trapz(recall, precision)) |>
  ungroup() |>
  group_by(pipeline) |>
  summarize(aauc = mean(cum_area))

summary_pr_curve_df <- pr_curve_df |>
  group_by(abundance_threshold, pipeline) |>
  summarise(
    mean_precision = mean(precision), 
    mean_recall = mean(recall),
    sd_precision = sd(precision),
    sd_recall = sd(recall),
    mean_F1 = mean(F1),
    sd_F1 = sd(F1)
  )
p5 <- ggplot(summary_pr_curve_df, aes(x = mean_recall, y = mean_precision, fill = pipeline)) + 
  geom_ribbon(aes(xmin = mean_recall - sd_recall, xmax = mean_recall + sd_recall,
                  ymin = mean_precision - sd_precision, ymax = mean_precision + sd_precision), alpha = 0.5) + 
  geom_path(aes(color = pipeline)) +
  facet_grid(rows = vars(pipeline))+ 
  geom_text(
    data = aauc, 
    aes(x = 0.125, y = ifelse(aauc>0.5, 0.25, 0.75),
        label = paste0("AAUC = ", round(aauc,3))),
    inherit.aes = FALSE,
    size = 5 / (14/5)
  ) + 
  theme(text=element_text(size=5),
        legend.position = "none")
p5

pr_curve_df_0 <- pr_curve_df |> 
  filter(abundance_threshold == 0)
p6 <- ggplot(pr_curve_df_0, aes(x = pipeline, y = precision)) +
  geom_boxplot()
p6

p7 <- ggplot(pr_curve_df_0, aes(x = pipeline, y = recall)) +
  geom_boxplot()
p7

p8 <- ggplot(pr_curve_df_0, aes(x = pipeline, y = F1)) +
  geom_boxplot()
p8

###############################################################################
ggsave(filename="figures/p5_pr_curve_2.svg", plot=p5, dpi=450, 
       width=60,height=75,units="mm",device="svg")

pr_curve_df |>
  group_by(name, pipeline) |>
  arrange(recall) |>
  summarize(cum_area = trapz(recall, precision)) |>
  ggplot(aes(x = pipeline, y = cum_area)) + 
  geom_boxplot()
# F1 Curve
max_f1_thresholds <- pr_curve_df |>
  group_by(pipeline, name) |>
  filter(F1 == max(F1, na.rm = TRUE)) |>
  summarise(
    abundance_threshold = mean(abundance_threshold, na.rm = TRUE),
    F1 = F1
  ) |>
  ungroup() |>
  group_by(pipeline) |>
  summarise(abundance_threshold = mean(abundance_threshold), 
            F1 = mean(F1))

p6 <- ggplot(summary_pr_curve_df, aes(x = abundance_threshold, y = mean_F1, fill = pipeline)) + 
  geom_ribbon(aes(ymin = mean_F1 - sd_F1, ymax = mean_F1 + sd_F1), alpha = 0.5) + 
  geom_path(aes(color = pipeline)) +
  facet_grid(rows = vars(pipeline))+ 
  geom_vline(
    data = max_f1_thresholds,
    aes(xintercept = abundance_threshold),
    linetype = "dashed",
    linewidth = 0.5
  )  + 
  geom_text(
    data = max_f1_thresholds, 
    aes(x = abundance_threshold + 0.015, y = ifelse(F1<0.4, 0.7, 0.5) ,
        label = paste0("F1 = ", round(F1,3), "\nThreshold = ", round(abundance_threshold,3))),
    size = 5 / (14/5),
    inherit.aes = FALSE
  ) + 
  geom_vline(
    aes(xintercept = 1/21), 
    linetype = "dotted",
    linewidth = 0.5
  ) +
  theme(text=element_text(size=5), 
        legend.position = "none")
p6

ggsave(filename="figures/p6_f1_curve.svg", plot=p6, dpi=450, 
       width=60,height=75,units="mm",device="svg")

pr_curve_df |> 
  filter(precision == 0, recall == 0) |>
  group_by(pipeline) |>
  summarise(min_abund = min(abundance_threshold), 
            max_abund = max(abundance_threshold))