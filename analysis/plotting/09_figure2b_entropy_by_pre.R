#!/usr/bin/env Rscript

# Plot Figure 2B: ancestry entropy by assigned PRE and reporting status.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "plot_setup.R"))
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
pre_levels <- c("Hispanic", "Black", "EAS", "Middle Eastern", "Other/Unknown")

plot_data <- read.csv(
  file.path(source_dir, "figure3_entropy_source_restricted_internal.csv")
) |>
  transmute(
    assigned_pre = factor(assigned_sre, levels = pre_levels),
    status = factor(unname(status_map[reporting_status]), levels = unname(status_map)),
    entropy_bits
  ) |>
  filter(!is.na(assigned_pre))
counts <- plot_data |> count(assigned_pre, status) |> tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0)
labels <- setNames(
  sprintf("%s\nn = %d/%d", counts$assigned_pre, counts[[2]], counts[[3]]),
  counts$assigned_pre
)

panel <- ggplot(plot_data, aes(assigned_pre, entropy_bits, color = status, fill = status)) +
  geom_boxplot(
    width = 0.62, position = position_dodge(width = 0.72),
    outlier.shape = NA, alpha = 0.16, linewidth = 0.55
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.13, dodge.width = 0.72, seed = 20260722L
    ),
    alpha = 0.48, size = 0.85, stroke = 0
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 23, fill = "white", size = 2.1,
    position = position_dodge(width = 0.72)
  ) +
  scale_color_manual(values = status_colors, guide = "none") +
  scale_fill_manual(values = status_colors, name = "PRE reporting status") +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  coord_cartesian(ylim = c(0, log2(5))) +
  labs(x = "Assigned PRE", y = "Shannon entropy (bits)") +
  theme_panel(9.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_panel(panel, figure_dir, "Figure2B_entropy_by_PRE", 7.5, 5)
