#!/usr/bin/env Rscript

# Plot Figure 1G: PRE-GIA agreement across ancestry-concentration thresholds.

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "..", "R", "plot_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tidyr)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "descriptive", "tables")
figure_dir <- file.path(results_dir, "figure1", "figures")

sensitivity <- read.csv(file.path(source_dir, "maximum_gia_threshold_sensitivity.csv"))
plot_data <- sensitivity |>
  pivot_longer(c(observed_agreement, kappa), names_to = "metric", values_to = "estimate") |>
  mutate(metric = factor(
    metric,
    levels = c("observed_agreement", "kappa"),
    labels = c("Observed agreement", "Cohen's kappa")
  ))
labels <- sensitivity |> filter(round(threshold, 2) %in% c(0, 0.50, 0.60, 0.70, 0.80, 0.90))

# The vertical guide marks the prespecified 70% concentration threshold.
panel <- ggplot(plot_data, aes(threshold, estimate, color = metric, shape = metric)) +
  geom_vline(xintercept = 0.70, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.1) +
  geom_text(
    data = labels,
    aes(threshold, 0.685, label = paste0("n=", n_retained)),
    inherit.aes = FALSE,
    size = 2.35,
    angle = 45,
    hjust = 0
  ) +
  scale_color_manual(values = c("Observed agreement" = "#0072B2", "Cohen's kappa" = "#D55E00")) +
  scale_shape_manual(values = c("Observed agreement" = 16, "Cohen's kappa" = 17)) +
  scale_x_continuous(
    limits = c(0, 0.92),
    breaks = c(0, 0.5, 0.6, 0.7, 0.8, 0.9),
    labels = c("All", ">50%", ">60%", ">70%", ">80%", ">90%")
  ) +
  scale_y_continuous(
    limits = c(0.67, 1.02),
    breaks = seq(0.7, 1.0, 0.1),
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    x = "Threshold on largest-component GIA proportion",
    y = "Agreement estimate",
    color = NULL,
    shape = NULL
  ) +
  theme_panel(8.7) +
  theme(legend.position = "bottom")

save_panel(panel, figure_dir, "Figure1G_threshold_sensitivity", 6, 4.4)
