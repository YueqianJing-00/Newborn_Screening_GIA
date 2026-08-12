#!/usr/bin/env Rscript

# Plot Figure 2A: ancestry entropy by PRE reporting status.

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "..", "R", "plot_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "descriptive", "tables")
figure_dir <- file.path(results_dir, "figure2", "figures")
status_map <- c("Single/no multiple report" = "Single/no multiple PRE", Multiple = "Multiple PREs")
status_colors <- c("Single/no multiple PRE" = "#0072B2", "Multiple PREs" = "#D55E00")

plot_data <- read.csv(
  file.path(source_dir, "figure3_entropy_source_restricted_internal.csv")
) |>
  transmute(
    status = factor(unname(status_map[reporting_status]), levels = unname(status_map)),
    entropy_bits
  )
counts <- plot_data |> count(status)
x_labels <- setNames(paste0(counts$status, "\n(n = ", counts$n, ")"), counts$status)

# Violin, boxplot, and points show complementary views of the same distribution.
panel <- ggplot(plot_data, aes(status, entropy_bits, color = status, fill = status)) +
  geom_violin(width = 0.78, alpha = 0.16, linewidth = 0.55, trim = TRUE) +
  geom_boxplot(width = 0.20, outlier.shape = NA, alpha = 0.34, linewidth = 0.55) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 20260722L),
    alpha = 0.34, size = 0.85, stroke = 0
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, fill = "white", color = "black", size = 2.4) +
  scale_color_manual(values = status_colors, guide = "none") +
  scale_fill_manual(values = status_colors, guide = "none") +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  coord_cartesian(ylim = c(0, log2(5))) +
  labs(x = NULL, y = "Shannon entropy (bits)") +
  theme_panel(9.5)

save_panel(panel, figure_dir, "Figure2A_entropy_overall", 5, 5)
