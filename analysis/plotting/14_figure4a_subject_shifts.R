#!/usr/bin/env Rscript

# Plot Figure 4A: subject-level prediction shifts toward the observed outcome.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "plot_setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "figure4_analysis", "runs", "main_prediction_shifts", "tables")
figure_dir <- file.path(results_dir, "figure4", "figures")
plot_data <- fread(file.path(source_dir, "figure4_subject_shift_plot_source_restricted_internal.csv"))
label_order <- plot_data[, unique(outcome_label), by = outcome][
  match(c("FP", "TP"), outcome), V1
]
plot_data[, outcome_label := factor(
  outcome_label,
  levels = label_order
)]

# Sign each probability change so positive values always favor the observed outcome.
panel <- ggplot(plot_data, aes(shift_pp, outcome_label, color = outcome)) +
  geom_vline(xintercept = 0, linewidth = 0.36, linetype = "22", color = "#7A7A7A") +
  geom_boxplot(
    aes(group = outcome_label), width = 0.22, outlier.shape = NA,
    fill = NA, color = "#333333", linewidth = 0.38
  ) +
  geom_jitter(
    size = 1.25, alpha = 0.58, stroke = 0,
    position = position_jitter(width = 0, height = 0.105, seed = 20260722L)
  ) +
  scale_color_manual(values = c(TP = "#0072B2", FP = "#D55E00")) +
  scale_x_continuous(limits = c(-18, 18), breaks = seq(-15, 15, 5)) +
  labs(x = "Change toward observed outcome (percentage points)", y = NULL) +
  theme_panel(8) +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_panel(panel, figure_dir, "Figure4A_subject_prediction_shifts", 4.1, 3.1)
