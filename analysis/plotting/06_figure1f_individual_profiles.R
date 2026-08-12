#!/usr/bin/env Rscript

# Plot Figure 1F: individual continuous GIA profiles grouped by assigned PRE.

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
pre_levels <- c(
  "Hispanic", "White", "Middle Eastern", "Black", "SAS", "EAS",
  "Native American", "Other/Unknown"
)

cohort <- read.csv(
  file.path(source_dir, "figure1_cohort_admixture_source_restricted_internal.csv"),
  check.names = FALSE
) |>
  mutate(
    assigned_sre = factor(assigned_sre, levels = pre_levels),
    plot_index = anonymous_plot_index
  )
groups <- group_boundaries(cohort, "assigned_sre", "plot_index")
plot_data <- cohort |>
  pivot_longer(all_of(ancestry_levels), names_to = "component", values_to = "proportion") |>
  mutate(component = factor(component, levels = ancestry_levels))

panel <- ggplot(plot_data, aes(plot_index, proportion, fill = component)) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(
    data = groups[-nrow(groups), , drop = FALSE],
    aes(xintercept = end + 0.5),
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = ancestry_palette, breaks = ancestry_levels, drop = FALSE) +
  scale_x_continuous(
    breaks = groups$center,
    labels = as.character(groups$assigned_sre),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(x = "Assigned PRE", y = "GIA proportion", fill = "GIA component") +
  theme_panel(8.7) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1),
    legend.position = "right",
    legend.direction = "vertical"
  )

save_panel(panel, figure_dir, "Figure1F_individual_profiles", 12, 4.4)
