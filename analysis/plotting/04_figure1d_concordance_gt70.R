#!/usr/bin/env Rscript

# Plot Figure 1D: PRE-GIA cross-classification above 70% GIA concentration.

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "..", "R", "plot_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "descriptive", "tables")
figure_dir <- file.path(results_dir, "figure1", "figures")
pre_levels <- c(
  "Hispanic", "Black", "White", "SAS", "EAS", "Middle Eastern",
  "Native American", "Other/Unknown"
)

full_counts <- read.csv(file.path(source_dir, "figure1_cross_classification_all.csv"))
counts <- read.csv(file.path(source_dir, "figure1_cross_classification_gt70.csv"))
fill_max <- max(full_counts$n)
plot_data <- expand_grid(
  assigned_pre = pre_levels,
  largest_gia_component = ancestry_levels
) |>
  left_join(counts, by = c("assigned_pre", "largest_gia_component")) |>
  mutate(
    n = replace_na(n, 0L),
    assigned_pre_plot = factor(assigned_pre, levels = rev(pre_levels)),
    largest_gia_component = factor(largest_gia_component, levels = ancestry_levels),
    label = ifelse(n == 0, "", n),
    text_color = ifelse(n >= 0.35 * fill_max, "white", "black")
  )

# Reuse the full-cohort fill limit for a directly comparable color scale.
panel <- ggplot(plot_data, aes(largest_gia_component, assigned_pre_plot, fill = n)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = label, color = text_color), size = 2.8) +
  scale_color_identity() +
  scale_fill_gradient(
    low = "#F7FBFF", high = "#0072B2", limits = c(0, fill_max), name = "Count"
  ) +
  scale_x_discrete(limits = ancestry_levels, drop = FALSE) +
  coord_fixed(ratio = 0.62) +
  labs(x = "Largest GIA component (>70% subset)", y = "Assigned PRE") +
  theme_panel(8.5) +
  theme(panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.45))

save_panel(panel, figure_dir, "Figure1D_concordance_gt70", 5.2, 4.3)
