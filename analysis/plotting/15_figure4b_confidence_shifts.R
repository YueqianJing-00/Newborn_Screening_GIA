#!/usr/bin/env Rscript

# Plot Figure 4B: mean prediction shifts by GIA concentration and outcome.

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "..", "R", "plot_setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "figure4_analysis", "runs", "main_prediction_shifts", "tables")
figure_dir <- file.path(results_dir, "figure4", "figures")
plot_data <- fread(file.path(source_dir, "figure4_confidence_shift_plot_source_aggregate.csv"))
high_group <- plot_data[grepl(">", confidence_group, fixed = TRUE), unique(confidence_group)]
row_order <- c(
  plot_data[confidence_group == high_group & outcome == "TP", row_label],
  plot_data[confidence_group == high_group & outcome == "FP", row_label],
  plot_data[confidence_group != high_group & outcome == "TP", row_label],
  plot_data[confidence_group != high_group & outcome == "FP", row_label]
)
plot_data[, row_label := factor(row_label, levels = rev(row_order))]

# Points show mean raw score changes; horizontal bars are bootstrap 95% CIs.
panel <- ggplot(plot_data, aes(estimate_pp, row_label, color = outcome, shape = outcome)) +
  geom_vline(xintercept = 0, linewidth = 0.36, linetype = "22", color = "#7A7A7A") +
  geom_hline(yintercept = 2.5, linewidth = 0.36, color = "#D9D9D9") +
  geom_errorbarh(aes(xmin = low_pp, xmax = high_pp), height = 0.13, linewidth = 0.42) +
  geom_point(size = 2.25, stroke = 0.45, fill = "white") +
  geom_text(
    aes(x = 10.55, label = estimate_label),
    hjust = 1, color = "#333333", size = 2.45, show.legend = FALSE
  ) +
  scale_color_manual(values = c(TP = "#0072B2", FP = "#D55E00")) +
  scale_shape_manual(values = c(TP = 21, FP = 24)) +
  scale_x_continuous(limits = c(-8, 11), breaks = c(-5, 0, 5)) +
  labs(x = "Raw change in predicted TP probability\n(percentage points; 95% CI)", y = NULL) +
  theme_panel(8) +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_panel(panel, figure_dir, "Figure4B_GIA_confidence_shifts", 3.8, 3.1)
