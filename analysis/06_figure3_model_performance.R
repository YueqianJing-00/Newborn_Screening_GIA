#!/usr/bin/env Rscript

# Figure 3: cross-validated performance and held-out permutation importance.

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "..", "R", "project_setup.R"))
require_packages(c("data.table", "ggplot2", "patchwork"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

paths <- project_paths(script_path)
project_root <- paths$root
analysis_dir <- file.path(paths$results, "figure3")
source_analysis_dir <- file.path(
  paths$results, "mma_model", "runs", "main_117_top10_metabolites"
)
source_table_dir <- file.path(source_analysis_dir, "tables")
output_table_dir <- file.path(analysis_dir, "tables")
output_figure_dir <- file.path(analysis_dir, "figures")
make_directories(output_table_dir, output_figure_dir)

# Model source tables ----

rel_path <- function(path) {
  relative_to_project(path, project_root)
}

required_tables <- c(
  panel_a = "figure_panel_a_source.csv",
  panel_b = "figure_panel_b_source.csv",
  primary_performance = "primary_performance.csv",
  paired_effects = "paired_effects.csv",
  importance_summary = "permutation_importance_summary.csv",
  metabolite_selection = "metabolite_selection_summary.csv",
  cohort_summary = "cohort_summary.csv"
)
required_paths <- setNames(file.path(source_table_dir, required_tables), names(required_tables))
if (any(!file.exists(required_paths))) {
  stop(
    "Missing required source tables: ",
    paste(rel_path(required_paths[!file.exists(required_paths)]), collapse = ", ")
  )
}

performance_plot_source <- fread(required_paths[["panel_a"]])
importance_plot_source <- fread(required_paths[["panel_b"]])
primary_performance <- fread(required_paths[["primary_performance"]])
paired_effects <- fread(required_paths[["paired_effects"]])
importance_summary <- fread(required_paths[["importance_summary"]])
metabolite_selection <- fread(required_paths[["metabolite_selection"]])
cohort_summary <- fread(required_paths[["cohort_summary"]])

if (
  nrow(cohort_summary) != 1L ||
    cohort_summary$n != 117L ||
    cohort_summary$tp != 85L ||
    cohort_summary$fp != 32L ||
    cohort_summary$metabolite_candidates != 40L ||
    cohort_summary$metabolites_selected_per_fold != 10L
) {
  stop("Source analysis does not match the 117-subject/top-10 protocol.")
}
if (!all(c("FC", "C3_C2") %in% metabolite_selection$metabolite)) {
  stop("FC and C3/C2 are absent from the metabolite candidate ledger.")
}

model_order <- c(
  "Clinical + metabolites",
  "Clinical + metabolites + PRE",
  "Clinical + metabolites + GIA",
  "Clinical + metabolites + PRE + GIA"
)
if (!setequal(unique(performance_plot_source$model), model_order)) {
  stop("Panel A does not contain the four expected models.")
}
if (uniqueN(performance_plot_source$repeat_id) != 100L) {
  stop("Panel A does not contain 100 cross-validation repeats.")
}

if (!"GIA_grouped" %in% importance_plot_source$group_id) {
  stop("Panel B does not contain the grouped GIA importance row.")
}
if (any(grepl("^GIA_(AFR|AMR|EAS|EUR)$", importance_plot_source$group_id))) {
  stop("Panel B unexpectedly contains individual GIA importance rows.")
}

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
importance_colors <- c(
  Clinical = "#7A7A7A",
  Metabolite = "#56B4E9",
  PRE = "#0072B2",
  GIA = "#D55E00"
)

# Figure panels ----

theme_publication <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 4)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_blank(),
      legend.key.height = unit(0.45, "lines"),
      legend.key.width = unit(1.1, "lines"),
      panel.spacing = unit(1.25, "lines"),
      plot.margin = margin(7, 9, 7, 8)
    )
}

# Panel A: four models on the x-axis; performance on the y-axis.
performance_plot_source[, model := factor(model, levels = model_order)]
performance_plot_source[, model_display := model_labels[as.character(model)]]
performance_plot_source[, model_display := gsub("\\n", " ", model_display)]
performance_plot_source[, metric := factor(
  metric,
  levels = c("AUC", "Specificity at sensitivity >=0.95")
)]
fwrite(
  performance_plot_source,
  file.path(output_table_dir, "figure3_panel_a_vertical_boxplots.csv")
)

panel_a <- ggplot(
  performance_plot_source,
  aes(x = model, y = estimate, fill = model, color = model)
) +
  geom_boxplot(
    width = 0.58,
    alpha = 0.68,
    linewidth = 0.55,
    outlier.alpha = 0.42,
    outlier.size = 1.3,
    outlier.stroke = 0.35
  ) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  scale_x_discrete(labels = model_labels) +
  scale_y_continuous(labels = function(y) sprintf("%.2f", y)) +
  labs(title = "A", x = NULL, y = "Performance estimate") +
  theme_publication(10.5) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(hjust = 0.5, vjust = 1)
  )

# Panel B: joint GIA permutation preserves the internal correlation and
# compositional structure of the four included GIA proportions.
importance_order <- importance_summary[
  group_id %in% unique(importance_plot_source$group_id)
][order(mean_delta_brier)]
importance_order[, display_label_plot := display_label]
importance_plot_source[, display_label_plot := display_label]
importance_plot_source[, display_label_plot := factor(
  display_label_plot,
  levels = importance_order$display_label_plot
)]
fwrite(
  importance_plot_source,
  file.path(output_table_dir, "figure3_panel_b_grouped_gia_importance.csv")
)

importance_means <- importance_plot_source[, .(
  mean_importance = mean(importance_delta_brier)
), by = .(display_label_plot, feature_group)]

panel_b <- ggplot(
  importance_plot_source,
  aes(x = importance_delta_brier, y = display_label_plot, fill = feature_group)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.64,
    linewidth = 0.45
  ) +
  geom_point(
    data = importance_means,
    aes(x = mean_importance, y = display_label_plot),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.05,
    stroke = 0.45,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = importance_colors) +
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +
  labs(title = "B", x = "Increase in held-out Brier score after permutation", y = NULL) +
  theme_publication(10) +
  theme(legend.position = "bottom")

figure3 <- panel_a / panel_b + plot_layout(heights = c(0.82, 1.18))

figure_pdf <- file.path(output_figure_dir, "Figure3_RF_model_comparison.pdf")
figure_png <- file.path(output_figure_dir, "Figure3_RF_model_comparison.png")
panel_b_pdf <- file.path(output_figure_dir, "Figure3B_grouped_GIA_importance.pdf")
panel_b_png <- file.path(output_figure_dir, "Figure3B_grouped_GIA_importance.png")
ggsave(figure_pdf, figure3, width = 11.5, height = 10.2, device = cairo_pdf)
ggsave(figure_png, figure3, width = 11.5, height = 10.2, dpi = 300, bg = "white")
ggsave(panel_b_pdf, panel_b, width = 11.5, height = 6.1, device = cairo_pdf)
ggsave(panel_b_png, panel_b, width = 11.5, height = 6.1, dpi = 300, bg = "white")

fwrite(primary_performance, file.path(output_table_dir, "primary_performance.csv"))
fwrite(paired_effects, file.path(output_table_dir, "paired_effects.csv"))
fwrite(importance_summary, file.path(output_table_dir, "permutation_importance_summary.csv"))

cat("Figure 3 PDF:", rel_path(figure_pdf), "\n")
cat("Figure 3 PNG:", rel_path(figure_png), "\n")
cat("Figure 3B grouped-GIA PDF:", rel_path(panel_b_pdf), "\n")
cat("Figure 3B grouped-GIA PNG:", rel_path(panel_b_png), "\n")
