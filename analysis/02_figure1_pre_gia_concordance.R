#!/usr/bin/env Rscript

# Figure 1: concordance between parent-reported ethnicity (PRE) and genetically
# inferred ancestry (GIA).

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
helper_dir <- file.path(dirname(script_path), "..", "R")
source(file.path(helper_dir, "project_setup.R"))
source(file.path(helper_dir, "statistical_helpers.R"))

require_packages(c(
  "dplyr", "ggplot2", "patchwork", "ragg", "scales", "tidyr"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

profiles_only <- "--profiles-only" %in% commandArgs(trailingOnly = TRUE)

paths <- project_paths(script_path)
project_root <- paths$root
analysis_dir <- file.path(paths$results, "figure1")
source_dir <- file.path(paths$results, "descriptive", "tables")
figure_dir <- file.path(analysis_dir, "figures")
table_dir <- file.path(analysis_dir, "tables")
make_directories(figure_dir, table_dir)

# Source tables ----

rel_path <- function(path) {
  relative_to_project(path, project_root)
}

cross_file <- file.path(source_dir, "figure2_sre_majority_ga_source.csv")
cohort_file <- file.path(
  source_dir,
  "figure1_cohort_admixture_source_restricted_internal.csv"
)
required_inputs <- c(cross_file, cohort_file)
require_files(required_inputs, "source file")

cross_source_reference <- read.csv(cross_file, check.names = FALSE, stringsAsFactors = FALSE)
cohort <- read.csv(cohort_file, check.names = FALSE, stringsAsFactors = FALSE)
if (sum(cross_source_reference$n) != nrow(cohort)) {
  stop("The aggregate and individual-level source tables have inconsistent totals.")
}

pre_levels <- c(
  "Hispanic", "White", "Middle Eastern", "Black", "SAS", "EAS",
  "Native American", "Other/Unknown"
)
heatmap_pre_levels <- c(
  "Hispanic", "Black", "White", "SAS", "EAS", "Middle Eastern",
  "Native American", "Other/Unknown"
)
gia_levels <- c("AMR", "AFR", "EUR", "SAS", "EAS")
pre_to_gia <- c(
  Hispanic = "AMR",
  White = "EUR",
  Black = "AFR",
  SAS = "SAS",
  EAS = "EAS"
)

ancestry_palette <- c(
  AMR = "#E69F00",
  AFR = "#D55E00",
  EUR = "#0072B2",
  SAS = "#009E73",
  EAS = "#CC79A7"
)

# Cross-classification and profile panels ----

cohort <- cohort %>%
  rename(
    assigned_pre = assigned_sre,
    largest_gia_component = majority_ga,
    maximum_gia_proportion = majority_ga_proportion
  ) %>%
  mutate(
    assigned_pre = factor(assigned_pre, levels = pre_levels),
    largest_gia_component = factor(largest_gia_component, levels = gia_levels)
  )

make_cross_table <- function(data) {
  counts <- data %>%
    transmute(
      assigned_pre = as.character(assigned_pre),
      largest_gia_component = as.character(largest_gia_component)
    ) %>%
    count(assigned_pre, largest_gia_component, name = "n")

  tidyr::expand_grid(
    assigned_pre = heatmap_pre_levels,
    largest_gia_component = gia_levels
  ) %>%
    left_join(counts, by = c("assigned_pre", "largest_gia_component")) %>%
    mutate(
      n = tidyr::replace_na(n, 0L),
      assigned_pre = factor(assigned_pre, levels = heatmap_pre_levels),
      largest_gia_component = factor(largest_gia_component, levels = gia_levels),
      assigned_pre_plot = factor(assigned_pre, levels = rev(heatmap_pre_levels))
    )
}

cross_all <- make_cross_table(cohort)
cross_gt70 <- make_cross_table(
  cohort %>% filter(maximum_gia_proportion > 0.70)
)

write.csv(cross_all, file.path(table_dir, "cross_classification_all.csv"), row.names = FALSE)
write.csv(cross_gt70, file.path(table_dir, "cross_classification_gt70.csv"), row.names = FALSE)

theme_manuscript <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(0.38, "cm"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

count_fill_max <- max(cross_all$n)

make_count_heatmap <- function(data, x_label, show_y = TRUE) {
  plot_data <- data %>%
    mutate(
      label = ifelse(n == 0, "", as.character(n)),
      text_color = ifelse(n >= 0.35 * count_fill_max, "white", "black")
    )

  ggplot(
    plot_data,
    aes(largest_gia_component, assigned_pre_plot, fill = n)
  ) +
    geom_tile(color = "white", linewidth = 0.55) +
    geom_text(aes(label = label, color = text_color), size = 2.8) +
    scale_color_identity() +
    scale_fill_gradient(
      low = "#F7FBFF",
      high = "#0072B2",
      limits = c(0, count_fill_max),
      name = "Count"
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = x_label,
      y = if (show_y) "Assigned PRE" else NULL
    ) +
    scale_x_discrete(limits = gia_levels, drop = FALSE) +
    coord_fixed(ratio = 0.62) +
    theme_manuscript(8.5) +
    theme(panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.45))
}

panel_all <- make_count_heatmap(
  cross_all,
  "Largest GIA component (all individuals)",
  show_y = TRUE
)
panel_gt70 <- make_count_heatmap(
  cross_gt70,
  "Largest GIA component (>70% subset)",
  show_y = FALSE
)

# Panel D: retain continuous information using the analytically corresponding
# GIA proportion in the five directly mapped PRE categories.
mapped_continuous <- cohort %>%
  filter(as.character(assigned_pre) %in% names(pre_to_gia)) %>%
  mutate(
    mapped_component = unname(pre_to_gia[as.character(assigned_pre)]),
    corresponding_gia_proportion = case_when(
      mapped_component == "AMR" ~ AMR,
      mapped_component == "AFR" ~ AFR,
      mapped_component == "EUR" ~ EUR,
      mapped_component == "SAS" ~ SAS,
      mapped_component == "EAS" ~ EAS
    ),
    display_group = factor(
      paste0(as.character(assigned_pre), " → ", mapped_component),
      levels = c(
        "Hispanic → AMR", "Black → AFR", "White → EUR",
        "SAS → SAS", "EAS → EAS"
      )
    )
  )

mapped_summary <- mapped_continuous %>%
  group_by(assigned_pre, mapped_component, display_group) %>%
  summarise(
    n = n(),
    mean = mean(corresponding_gia_proportion),
    median = median(corresponding_gia_proportion),
    q1 = quantile(corresponding_gia_proportion, 0.25),
    q3 = quantile(corresponding_gia_proportion, 0.75),
    .groups = "drop"
  )
write.csv(
  mapped_summary,
  file.path(table_dir, "continuous_corresponding_gia_summary.csv"),
  row.names = FALSE
)

panel_d <- ggplot(
  mapped_continuous,
  aes(display_group, corresponding_gia_proportion, fill = mapped_component)
) +
  geom_violin(width = 0.82, alpha = 0.18, trim = TRUE, linewidth = 0.45) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.42, linewidth = 0.45) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0, seed = 20260722L),
    alpha = 0.25, size = 0.7, stroke = 0
  ) +
  geom_text(
    data = mapped_summary,
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
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Estimated GIA proportion"
  ) +
  theme_manuscript(8.7) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1))

# Panel E: individual-level continuous GIA profiles for the full cohort.
expected_component <- c(
  Hispanic = "AMR", White = "EUR", Black = "AFR", SAS = "SAS", EAS = "EAS",
  `Native American` = "AMR"
)

cohort_ordered <- cohort %>%
  mutate(
    expected_component = unname(expected_component[as.character(assigned_pre)]),
    expected_proportion = case_when(
      expected_component == "AMR" ~ AMR,
      expected_component == "AFR" ~ AFR,
      expected_component == "EUR" ~ EUR,
      expected_component == "SAS" ~ SAS,
      expected_component == "EAS" ~ EAS,
      TRUE ~ maximum_gia_proportion
    )
  ) %>%
  arrange(assigned_pre, largest_gia_component, desc(expected_proportion), desc(maximum_gia_proportion)) %>%
  mutate(plot_index = row_number())

profile_boundaries <- cohort_ordered %>%
  group_by(assigned_pre) %>%
  summarise(
    start = min(plot_index),
    end = max(plot_index),
    center = (start + end) / 2,
    .groups = "drop"
  )

cohort_long <- cohort_ordered %>%
  select(plot_index, assigned_pre, all_of(gia_levels)) %>%
  pivot_longer(all_of(gia_levels), names_to = "GIA_component", values_to = "proportion") %>%
  mutate(GIA_component = factor(GIA_component, levels = gia_levels))

panel_e <- ggplot(cohort_long, aes(plot_index, proportion, fill = GIA_component)) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(
    data = profile_boundaries[-nrow(profile_boundaries), , drop = FALSE],
    aes(xintercept = end + 0.5),
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = ancestry_palette, breaks = gia_levels, drop = FALSE) +
  scale_x_continuous(
    breaks = profile_boundaries$center,
    labels = as.character(profile_boundaries$assigned_pre),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Assigned PRE",
    y = "GIA proportion",
    fill = "GIA component"
  ) +
  theme_manuscript(8.7) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1),
    legend.position = "right",
    legend.direction = "vertical"
  )

# Agreement across ancestry-concentration thresholds.
thresholds <- seq(0, 0.90, 0.05)
threshold_sensitivity <- bind_rows(lapply(thresholds, function(threshold) {
  subset <- cohort %>%
    filter(
      as.character(assigned_pre) %in% names(pre_to_gia),
      maximum_gia_proportion > threshold
    )
  mapped_pre <- unname(pre_to_gia[as.character(subset$assigned_pre)])
  metrics <- cohen_kappa(
    mapped_pre,
    as.character(subset$largest_gia_component),
    gia_levels
  )
  data.frame(
    threshold = threshold,
    n_retained = nrow(subset),
    observed_agreement = unname(metrics$observed),
    kappa = unname(metrics$kappa)
  )
}))
write.csv(
  threshold_sensitivity,
  file.path(table_dir, "maximum_gia_threshold_sensitivity.csv"),
  row.names = FALSE
)

threshold_long <- threshold_sensitivity %>%
  pivot_longer(
    c(observed_agreement, kappa),
    names_to = "metric",
    values_to = "estimate"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("observed_agreement", "kappa"),
      labels = c("Observed agreement", "Cohen's kappa")
    )
  )

label_thresholds <- c(0, 0.50, 0.60, 0.70, 0.80, 0.90)
threshold_labels <- threshold_sensitivity %>%
  filter(round(threshold, 2) %in% label_thresholds)

panel_f <- ggplot(
  threshold_long,
  aes(threshold, estimate, color = metric, shape = metric)
) +
  geom_vline(xintercept = 0.70, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.1) +
  geom_text(
    data = threshold_labels,
    aes(threshold, 0.685, label = paste0("n=", n_retained)),
    inherit.aes = FALSE,
    size = 2.35,
    angle = 45,
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0.70,
    y = 1.005,
    label = "largest GIA component >70%",
    size = 2.45,
    hjust = -0.04,
    vjust = 1
  ) +
  scale_color_manual(
    values = c("Observed agreement" = "#0072B2", "Cohen's kappa" = "#D55E00")
  ) +
  scale_shape_manual(values = c("Observed agreement" = 16, "Cohen's kappa" = 17)) +
  scale_x_continuous(
    limits = c(0, 0.92),
    breaks = c(0, 0.5, 0.6, 0.7, 0.8, 0.9),
    labels = c("All", ">50%", ">60%", ">70%", ">80%", ">90%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0.67, 1.02),
    breaks = seq(0.7, 1.0, 0.1),
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Threshold on largest-component GIA proportion",
    y = "Agreement estimate",
    color = NULL,
    shape = NULL
  ) +
  theme_manuscript(8.7) +
  theme(legend.position = "bottom")

# Save individual plots for flexible manuscript assembly.
panel_list <- list(
  all_cross = panel_all,
  gt70_cross = panel_gt70,
  continuous = panel_d,
  profiles = panel_e,
  threshold = panel_f
)
panel_stems <- c(
  all_cross = "PRE_GIA_cross_classification_all",
  gt70_cross = "PRE_GIA_cross_classification_gt70",
  continuous = "PRE_GIA_continuous_corresponding_proportion",
  profiles = "PRE_GIA_individual_profiles",
  threshold = "PRE_GIA_threshold_sensitivity"
)
panel_dimensions <- list(
  all_cross = c(5.2, 4.3),
  gt70_cross = c(5.2, 4.3),
  continuous = c(6.4, 4.4),
  profiles = c(12.0, 4.4),
  threshold = c(6.0, 4.4)
)
panel_names_to_save <- if (profiles_only) "profiles" else names(panel_list)
for (panel_name in panel_names_to_save) {
  dimensions <- panel_dimensions[[panel_name]]
  stem <- panel_stems[[panel_name]]
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    panel_list[[panel_name]],
    width = dimensions[1], height = dimensions[2], units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    panel_list[[panel_name]],
    width = dimensions[1], height = dimensions[2], units = "in", dpi = 300,
    device = ragg::agg_png, bg = "white"
  )
}

if (profiles_only) {
  message("Wrote individual-profile figure with right-side legend: ", figure_dir)
  quit(save = "no", status = 0)
}

# Assemble the manuscript layout ----

# Panel positions follow the manuscript figure layout.
# PowerPoint slide coordinates are English Metric Units (EMU) on a
# 12,192,000 x 6,858,000 canvas. Converting them to normalized coordinates
# preserves the aspect ratio of every original vector panel.
slide_width_emu <- 12192000
slide_height_emu <- 6858000
layout_boxes <- data.frame(
  panel = c("A", "B", "C", "D", "E"),
  content = c(
    "Individual GIA profiles",
    "PRE-by-largest-GIA counts, all individuals",
    "PRE-by-largest-GIA counts, >70% subset",
    "Corresponding continuous GIA proportions",
    "Agreement across largest-component thresholds"
  ),
  x_emu = c(0, 8681014, 0, 3671888, 8044628),
  y_from_top_emu = c(307118, 307118, 3723833, 3723833, 3723833),
  width_emu = c(8681014, 3456544, 3671888, 4416183, 4147372),
  height_emu = c(3175000, 2855807, 3033726, 3033726, 3033726),
  stringsAsFactors = FALSE
) %>%
  mutate(
    x = x_emu / slide_width_emu,
    y = 1 - (y_from_top_emu + height_emu) / slide_height_emu,
    width = width_emu / slide_width_emu,
    height = height_emu / slide_height_emu
  )
write.csv(
  layout_boxes,
  file.path(table_dir, "figure1_reference_layout_boxes.csv"),
  row.names = FALSE
)

combined_pdf <- file.path(
  figure_dir,
  "Figure1_PRE_GIA_concordance.pdf"
)
combined_png <- file.path(
  figure_dir,
  "Figure1_PRE_GIA_concordance.png"
)
composition_tex <- file.path(
  analysis_dir,
  "Figure1_PRE_GIA_concordance.tex"
)
composition_template <- file.path(project_root, "resources", "figure1_layout.tex.in")
if (!file.exists(composition_template)) {
  stop("Missing vector-composition template: ", rel_path(composition_template))
}
composition_lines <- readLines(composition_template, warn = FALSE)
composition_lines <- gsub(
  "@@FIGURE_DIR@@",
  normalizePath(figure_dir, mustWork = TRUE),
  composition_lines,
  fixed = TRUE
)
writeLines(composition_lines, composition_tex)

pdflatex_command <- Sys.which("pdflatex")
pdftoppm_command <- Sys.which("pdftoppm")
if (!nzchar(pdflatex_command) || !nzchar(pdftoppm_command)) {
  stop("The redraw requires pdflatex and pdftoppm on PATH.")
}

old_working_directory <- getwd()
setwd(project_root)
latex_output <- system2(
  pdflatex_command,
  args = c(
    "-interaction=nonstopmode",
    "-halt-on-error",
    paste0("-output-directory=", shQuote(figure_dir)),
    shQuote(composition_tex)
  ),
  stdout = TRUE,
  stderr = TRUE
)
setwd(old_working_directory)
latex_status <- attr(latex_output, "status")
if (!is.null(latex_status) && latex_status != 0L) {
  stop("Vector composition failed:\n", paste(latex_output, collapse = "\n"))
}
if (!file.exists(combined_pdf)) {
  stop("Vector composition completed without producing the expected PDF.")
}

png_prefix <- sub("\\.png$", "", combined_png)
render_output <- system2(
  pdftoppm_command,
  args = c(
    "-png", "-r", "300", "-singlefile",
    shQuote(combined_pdf),
    shQuote(png_prefix)
  ),
  stdout = TRUE,
  stderr = TRUE
)
render_status <- attr(render_output, "status")
if (!is.null(render_status) && render_status != 0L) {
  stop("PNG rendering failed:\n", paste(render_output, collapse = "\n"))
}
if (!file.exists(combined_png)) {
  stop("PDF rendering completed without producing the expected PNG.")
}

composition_stem <- tools::file_path_sans_ext(basename(composition_tex))
unlink(file.path(figure_dir, paste0(composition_stem, c(".aux", ".log"))))

message("Wrote combined figure: ", combined_png)
