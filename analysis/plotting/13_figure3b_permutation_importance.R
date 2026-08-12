#!/usr/bin/env Rscript

# Plot Figure 3B: held-out permutation importance in the full RF model.

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "..", "R", "plot_setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "mma_model", "runs", "main_117_top10_metabolites", "tables")
figure_dir <- file.path(results_dir, "figure3", "figures")
importance_colors <- c(
  Clinical = "#7A7A7A", Metabolite = "#56B4E9", PRE = "#0072B2", GIA = "#D55E00"
)

plot_data <- fread(file.path(source_dir, "figure_panel_b_source.csv"))
summary <- fread(file.path(source_dir, "permutation_importance_summary.csv"))
order <- summary[group_id %in% unique(plot_data$group_id)][order(mean_delta_brier)]$display_label
plot_data[, display_label := factor(display_label, levels = order)]
means <- plot_data[, .(mean_importance = mean(importance_delta_brier)), by = .(display_label, feature_group)]

# Positive values indicate worse held-out Brier score after predictor permutation.
panel <- ggplot(plot_data, aes(importance_delta_brier, display_label, fill = feature_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.64, linewidth = 0.45) +
  geom_point(
    data = means,
    aes(mean_importance, display_label),
    inherit.aes = FALSE,
    shape = 23, size = 2.05, fill = "white", color = "black"
  ) +
  scale_fill_manual(values = importance_colors) +
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +
  labs(x = "Increase in held-out Brier score after permutation", y = NULL, fill = NULL) +
  theme_panel(10) +
  theme(legend.position = "bottom")

save_panel(panel, figure_dir, "Figure3B_permutation_importance", 11.5, 6.1)
