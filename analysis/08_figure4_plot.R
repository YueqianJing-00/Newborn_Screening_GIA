#!/usr/bin/env Rscript

# Two-panel Figure 4 generated from the prediction-shift analysis.

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

parse_args <- function(args) {
  defaults <- list(
    source_run = "main_prediction_shifts",
    run_name = "final"
  )
  for (arg in args) {
    if (!grepl("^--[A-Za-z0-9_-]+=", arg)) stop("Malformed argument: ", arg)
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", parts[[1]])
    if (!key %in% names(defaults)) stop("Unknown argument: --", parts[[1]])
    if (!is_run_name(parts[[2]])) {
      stop(key, " contains unsupported characters.")
    }
    defaults[[key]] <- parts[[2]]
  }
  defaults
}

paths <- project_paths(script_path)
project_root <- paths$root
analysis_dir <- file.path(paths$results, "figure4_analysis")
args <- parse_args(commandArgs(trailingOnly = TRUE))
source_run <- file.path(analysis_dir, "runs", args$source_run)
source_tables <- file.path(source_run, "tables")
run_dir <- file.path(paths$results, "figure4", args$run_name)
table_dir <- file.path(run_dir, "tables")
figure_dir <- file.path(run_dir, "figures")
make_directories(table_dir, figure_dir)

# Prediction-shift sources ----

rel_path <- function(path) {
  relative_to_project(path, project_root)
}

source_files <- c(
  subject = file.path(source_tables, "subject_level_case_patterns_source_restricted_internal.csv"),
  confidence = file.path(source_tables, "gia_confidence_score_shift_summary.csv")
)
require_files(source_files, "source file")

subject <- fread(source_files[["subject"]])
confidence <- fread(source_files[["confidence"]])

if (nrow(subject) != 117L || sum(subject$outcome == "TP") != 85L || sum(subject$outcome == "FP") != 32L) {
  stop("Subject source did not reproduce the 117-newborn cohort.")
}
if (nrow(confidence) != 4L || !setequal(confidence$outcome, c("TP", "FP"))) {
  stop("GIA-confidence summary did not reproduce the four expected outcome-stratum cells.")
}

tp_color <- "#0072B2"
fp_color <- "#D55E00"
neutral_dark <- "#333333"
neutral_mid <- "#7A7A7A"
neutral_light <- "#D9D9D9"

# Figure panels ----

theme_journal <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      line = element_line(linewidth = 0.36, color = neutral_dark),
      axis.line = element_line(linewidth = 0.36, color = neutral_dark),
      axis.ticks = element_line(linewidth = 0.36, color = neutral_dark),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.text = element_text(size = 7, color = neutral_dark),
      axis.title = element_text(size = 8, color = neutral_dark),
      legend.text = element_text(size = 7, color = neutral_dark),
      legend.title = element_blank(),
      plot.tag = element_text(size = 9, face = "bold", color = neutral_dark),
      plot.tag.position = c(0, 1),
      plot.margin = margin(4, 5, 4, 5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# Panel A: subject-level changes signed toward the observed outcome.
subject[, outcome_label := factor(
  outcome,
  levels = c("FP", "TP"),
  labels = c("False positive (n = 32)", "True positive (n = 85)")
)]
subject[, shift_pp := 100 * correct_direction_shift]

panel_a_source <- subject[, .(
  anonymous_plot_index = seq_len(.N),
  outcome,
  outcome_label,
  shift_pp
)]
fwrite(
  panel_a_source,
  file.path(table_dir, "figure4_panel_A_source_restricted_internal.csv")
)

panel_a <- ggplot(
  panel_a_source,
  aes(x = shift_pp, y = outcome_label, color = outcome)
) +
  geom_vline(xintercept = 0, linewidth = 0.36, linetype = "22", color = neutral_mid) +
  geom_boxplot(
    aes(group = outcome_label),
    width = 0.22,
    outlier.shape = NA,
    fill = NA,
    color = neutral_dark,
    linewidth = 0.38
  ) +
  geom_jitter(
    size = 1.25,
    alpha = 0.58,
    stroke = 0,
    position = position_jitter(width = 0, height = 0.105, seed = 20260722L)
  ) +
  scale_color_manual(values = c(TP = tp_color, FP = fp_color)) +
  scale_x_continuous(
    limits = c(-18, 18),
    breaks = seq(-15, 15, by = 5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    tag = "A",
    x = "Change toward observed outcome (percentage points)",
    y = NULL
  ) +
  theme_journal() +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 8, 5, 6)
  )

# Panel B: mean raw score changes by GIA concentration and observed outcome.
confidence[, confidence_group := fifelse(
  grepl(">70%", majority_GIA_confidence, fixed = TRUE),
  ">70%",
  paste0(intToUtf8(8804), "70%")
)]
confidence[, outcome_name := fifelse(outcome == "TP", "TP", "FP")]
confidence[, row_label := paste0(confidence_group, ", ", outcome_name, " (n = ", n, ")")]
row_order <- c(
  ">70%, TP (n = 71)",
  ">70%, FP (n = 15)",
  paste0(intToUtf8(8804), "70%, TP (n = 14)"),
  paste0(intToUtf8(8804), "70%, FP (n = 17)")
)
confidence[, row_label := factor(row_label, levels = rev(row_order))]
confidence[, `:=`(
  estimate_pp = 100 * mean_probability_change,
  low_pp = 100 * change_ci_low,
  high_pp = 100 * change_ci_high,
  estimate_label = sprintf("%+.2f", 100 * mean_probability_change)
)]
fwrite(
  confidence[, .(
    confidence_group, outcome, n, estimate_pp, low_pp, high_pp, estimate_label
  )],
  file.path(table_dir, "figure4_panel_B_source_aggregate.csv")
)

panel_b <- ggplot(
  confidence,
  aes(x = estimate_pp, y = row_label, color = outcome, shape = outcome)
) +
  geom_vline(xintercept = 0, linewidth = 0.36, linetype = "22", color = neutral_mid) +
  geom_hline(yintercept = 2.5, linewidth = 0.36, color = neutral_light) +
  geom_errorbarh(
    aes(xmin = low_pp, xmax = high_pp),
    height = 0.13,
    linewidth = 0.42
  ) +
  geom_point(size = 2.25, stroke = 0.45, fill = "white") +
  geom_text(
    aes(x = 10.55, label = estimate_label),
    hjust = 1,
    color = neutral_dark,
    size = 2.45,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(TP = tp_color, FP = fp_color)) +
  scale_shape_manual(values = c(TP = 21, FP = 24)) +
  scale_x_continuous(
    limits = c(-8.0, 11.0),
    breaks = c(-5, 0, 5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "B",
    x = "Raw change in predicted TP probability\n(percentage points; 95% CI)",
    y = NULL
  ) +
  theme_journal() +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 6, 5, 8)
  )

figure <- (panel_a | panel_b) +
  plot_layout(widths = c(1.05, 0.95)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

width_in <- 180 / 25.4
height_in <- 78 / 25.4
pdf_path <- file.path(figure_dir, "Figure4_GIA_case_patterns_two_panel.pdf")
png_path <- file.path(figure_dir, "Figure4_GIA_case_patterns_two_panel_300dpi.png")
tiff_path <- file.path(figure_dir, "Figure4_GIA_case_patterns_two_panel_600dpi.tiff")

ggsave(pdf_path, figure, width = width_in, height = height_in, device = cairo_pdf, bg = "white")
ggsave(png_path, figure, width = width_in, height = height_in, dpi = 300, bg = "white")
ggsave(
  tiff_path,
  figure,
  width = width_in,
  height = height_in,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

cat("Two-panel Figure 4 redraw complete\n")
cat("PDF:", rel_path(pdf_path), "\n")
cat("PNG:", rel_path(png_path), "\n")
cat("TIFF:", rel_path(tiff_path), "\n")
