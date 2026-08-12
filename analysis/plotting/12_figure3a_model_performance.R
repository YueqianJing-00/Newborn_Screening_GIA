#!/usr/bin/env Rscript

# Plot Figure 3A: repeated cross-validation performance for four RF models.

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "plot_setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_root <- project_root_from_script()
results_dir <- results_root(project_root)
source_dir <- file.path(results_dir, "mma_model", "runs", "main_117_top10_metabolites", "tables")
figure_dir <- file.path(results_dir, "figure3", "figures")
model_order <- c(
  "Clinical + metabolites",
  "Clinical + metabolites + PRE",
  "Clinical + metabolites + GIA",
  "Clinical + metabolites + PRE + GIA"
)
model_labels <- c(
  "Clinical + metabolites" = "Covariates\nonly",
  "Clinical + metabolites + PRE" = "+ PRE",
  "Clinical + metabolites + GIA" = "+ GIA",
  "Clinical + metabolites + PRE + GIA" = "Both"
)
model_colors <- c(
  "Clinical + metabolites" = "#6F6F6F",
  "Clinical + metabolites + PRE" = "#0072B2",
  "Clinical + metabolites + GIA" = "#D55E00",
  "Clinical + metabolites + PRE + GIA" = "#009E73"
)

plot_data <- fread(file.path(source_dir, "figure_panel_a_source.csv"))
plot_data[, model := factor(model, levels = model_order)]
plot_data[, metric := factor(metric, levels = c("AUC", "Specificity at sensitivity >=0.95"))]

# Boxplots summarize algorithmic stability across dependent repeated-CV runs.
panel <- ggplot(plot_data, aes(model, estimate, fill = model, color = model)) +
  geom_boxplot(width = 0.58, alpha = 0.68, linewidth = 0.55, outlier.alpha = 0.42) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  scale_x_discrete(labels = model_labels) +
  scale_y_continuous(labels = function(y) sprintf("%.2f", y)) +
  labs(x = NULL, y = "Performance estimate") +
  theme_panel(10.5) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

save_panel(panel, figure_dir, "Figure3A_model_performance", 11.5, 4)
