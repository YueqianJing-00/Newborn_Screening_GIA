#!/usr/bin/env Rscript

# Plot Figure 2D: reported PRE selections for the multiple-PRE cohort.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "plot_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "descriptive", "tables")
figure_dir <- file.path(results_dir, "figure2", "figures")
pre_levels <- c(
  "Hispanic", "White", "Black", "EAS", "SAS", "Middle Eastern",
  "Native American", "Other/Unknown"
)
pre_palette <- c(
  Hispanic = ancestry_palette[["AMR"]], White = ancestry_palette[["EUR"]],
  Black = ancestry_palette[["AFR"]], EAS = ancestry_palette[["EAS"]],
  SAS = ancestry_palette[["SAS"]], `Middle Eastern` = "#56B4E9",
  `Native American` = "#F0E442", `Other/Unknown` = "#8C8C8C"
)

tiles <- read.csv(
  file.path(source_dir, "figure4_multisre_selection_tiles_source_restricted_internal.csv")
) |>
  mutate(
    profile_index = anonymous_plot_index,
    reported_category = factor(reported_category, levels = rev(pre_levels)),
    tile_fill = ifelse(as.logical(present), pre_palette[as.character(reported_category)], "#F2F2F2")
  )

# Colored tiles mark the PRE categories selected by each anonymous participant.
panel <- ggplot(tiles, aes(profile_index, reported_category, fill = tile_fill)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0.5, max(tiles$profile_index) + 0.5),
    breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Participants in Figure 2C order", y = "Reported PRE selection") +
  theme_panel(9.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
  )

save_panel(panel, figure_dir, "Figure2D_multiple_PRE_selections", 13.5, 3.2)
