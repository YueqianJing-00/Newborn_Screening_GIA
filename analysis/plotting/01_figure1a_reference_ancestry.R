#!/usr/bin/env Rscript

# Plot Figure 1A: ancestry profiles of selected 1000 Genomes references.

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

reference <- read.csv(
  file.path(source_dir, "figure1_reference_admixture_source_restricted_internal.csv"),
  check.names = FALSE
) |>
  mutate(
    Population = factor(Population, levels = unique(Population)),
    plot_index = anonymous_plot_index
  )

groups <- group_boundaries(reference, "Population", "plot_index")
plot_data <- reference |>
  pivot_longer(all_of(ancestry_levels), names_to = "ancestry", values_to = "proportion") |>
  mutate(ancestry = factor(ancestry, levels = ancestry_levels))

# Stack ancestry components for each anonymous reference in population order.
panel <- ggplot(plot_data, aes(plot_index, proportion, fill = ancestry)) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(
    data = groups[-nrow(groups), , drop = FALSE],
    aes(xintercept = end + 0.5),
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = ancestry_palette, breaks = ancestry_levels, drop = FALSE) +
  scale_x_continuous(
    breaks = groups$center,
    labels = as.character(groups$Population),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Ancestry component") +
  theme_panel(9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    legend.position = "bottom"
  )

save_panel(panel, figure_dir, "Figure1A_reference_ancestry", 12, 3.7)
