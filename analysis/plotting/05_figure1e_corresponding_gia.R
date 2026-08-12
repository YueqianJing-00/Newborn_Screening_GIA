#!/usr/bin/env Rscript

# Plot Figure 1E: continuous GIA proportion corresponding to each mapped PRE.

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "..", "R", "plot_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "descriptive", "tables")
figure_dir <- file.path(results_dir, "figure1", "figures")

plot_data <- read.csv(
  file.path(source_dir, "figure1_corresponding_gia_source_restricted_internal.csv")
) |>
  mutate(display_group = factor(
    display_group,
    levels = c("Hispanic → AMR", "Black → AFR", "White → EUR", "SAS → SAS", "EAS → EAS")
  ))
counts <- plot_data |> count(display_group)

# Show individual values with violin and boxplot distribution summaries.
panel <- ggplot(
  plot_data,
  aes(display_group, corresponding_gia_proportion, fill = mapped_component)
) +
  geom_violin(width = 0.82, alpha = 0.18, trim = TRUE, linewidth = 0.45) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.42, linewidth = 0.45) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0, seed = 20260722L),
    alpha = 0.25, size = 0.7, stroke = 0
  ) +
  geom_text(
    data = counts,
    aes(display_group, 1.02, label = paste0("n=", n)),
    inherit.aes = FALSE,
    size = 2.6,
    vjust = 0
  ) +
  scale_fill_manual(values = ancestry_palette, guide = "none") +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(x = NULL, y = "Estimated GIA proportion") +
  theme_panel(8.7) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1))

save_panel(panel, figure_dir, "Figure1E_corresponding_GIA", 6.4, 4.4)
