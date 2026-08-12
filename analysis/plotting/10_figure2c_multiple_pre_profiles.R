#!/usr/bin/env Rscript

# Plot Figure 2C: GIA profiles among participants reporting multiple PREs.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "plot_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tidyr)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "descriptive", "tables")
figure_dir <- file.path(results_dir, "figure2", "figures")
pre_levels <- c(
  "Hispanic", "Black", "EAS", "SAS", "Middle Eastern",
  "Native American", "White", "Other/Unknown"
)

profile <- read.csv(
  file.path(source_dir, "figure4_multisre_admixture_source_restricted_internal.csv"),
  check.names = FALSE
) |>
  mutate(
    assigned_sre = factor(assigned_sre, levels = pre_levels),
    profile_index = anonymous_plot_index
  )
groups <- group_boundaries(profile, "assigned_sre", "profile_index") |>
  mutate(axis_label = sprintf("%s (n = %d)", assigned_sre, end - start + 1))
plot_data <- profile |>
  pivot_longer(all_of(ancestry_levels), names_to = "component", values_to = "proportion") |>
  mutate(component = factor(component, levels = ancestry_levels))

# Each column is one anonymous participant and sums to one across GIA components.
panel <- ggplot(plot_data, aes(profile_index, proportion, fill = component)) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(
    data = groups[-nrow(groups), , drop = FALSE],
    aes(xintercept = end + 0.5),
    color = "white",
    linewidth = 0.65
  ) +
  scale_fill_manual(values = ancestry_palette, breaks = ancestry_levels, name = "GIA component") +
  scale_x_continuous(
    breaks = groups$center,
    labels = groups$axis_label,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25), labels = label_percent(),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(x = "Assigned PRE", y = "GIA proportion") +
  theme_panel(9.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

save_panel(panel, figure_dir, "Figure2C_multiple_PRE_profiles", 13.5, 4.5)
