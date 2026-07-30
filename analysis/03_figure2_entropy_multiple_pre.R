#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
helper_dir <- file.path(dirname(script_path), "..", "R")
source(file.path(helper_dir, "project_setup.R"))
source(file.path(helper_dir, "plot_helpers.R"))

required_packages <- c(
  "dplyr", "ggplot2", "patchwork", "ragg", "scales", "tidyr", "digest"
)
require_packages(required_packages)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

paths <- project_paths(script_path)
project_root <- paths$root
run_dir <- file.path(paths$results, "figure2")
figure_dir <- file.path(run_dir, "figures")
table_dir <- file.path(run_dir, "tables")
make_directories(figure_dir, table_dir)

# Source tables ----

source_table_dir <- file.path(paths$results, "descriptive", "tables")
input_files <- c(
  entropy_individual = "figure3_entropy_source_restricted_internal.csv",
  entropy_summary = "figure3_entropy_summary.csv",
  entropy_overall = "entropy_overall_single_vs_multiple_bootstrap.csv",
  entropy_within_pre = "entropy_by_sre_single_vs_multiple_tests.csv",
  profile = "figure4_multisre_admixture_source_restricted_internal.csv",
  selection_tiles = "figure4_multisre_selection_tiles_source_restricted_internal.csv",
  component_mapping = "ancestry_component_mapping.csv"
)
input_paths <- file.path(source_table_dir, unname(input_files))
names(input_paths) <- names(input_files)
require_files(input_paths)

read_source_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

entropy_individual <- read_source_csv(input_paths[["entropy_individual"]])
entropy_summary <- read_source_csv(input_paths[["entropy_summary"]])
entropy_overall <- read_source_csv(input_paths[["entropy_overall"]])
entropy_within_pre <- read_source_csv(input_paths[["entropy_within_pre"]])
profile_source <- read_source_csv(input_paths[["profile"]])
tile_source <- read_source_csv(input_paths[["selection_tiles"]])
component_mapping <- read_source_csv(input_paths[["component_mapping"]])

ancestry_levels <- c("AMR", "AFR", "EUR", "SAS", "EAS")
study_mapping <- component_mapping %>%
  filter(source_file == "1000G_378.5.Q") %>%
  arrange(match(raw_column, paste0("V", seq_along(ancestry_levels))))
if (!identical(study_mapping$ancestry_component, ancestry_levels)) {
  stop("The saved ancestry-component mapping does not match AMR/AFR/EUR/SAS/EAS.")
}

required_entropy_columns <- c(
  "anonymous_plot_index", "assigned_sre", "reporting_status", "entropy_bits"
)
required_profile_columns <- c(
  "anonymous_plot_index", "assigned_sre", "reported_sre_combination",
  "majority_ga", "majority_ga_proportion", ancestry_levels
)
required_tile_columns <- c(
  "anonymous_plot_index", "reported_category", "present"
)
if (!all(required_entropy_columns %in% names(entropy_individual))) {
  stop("Unexpected columns in the entropy individual-level source table.")
}
if (!all(required_profile_columns %in% names(profile_source))) {
  stop("Unexpected columns in the multiple-PRE GIA source table.")
}
if (!all(required_tile_columns %in% names(tile_source))) {
  stop("Unexpected columns in the multiple-PRE selection-tile source table.")
}
if (nrow(entropy_individual) != 378L) {
  stop("Expected 378 entropy rows; found ", nrow(entropy_individual), ".")
}
if (nrow(profile_source) != 52L) {
  stop("Expected 52 multiple-PRE profiles; found ", nrow(profile_source), ".")
}
if (anyDuplicated(profile_source$anonymous_plot_index)) {
  stop("Anonymous profile indices are not unique.")
}
if (nrow(tile_source) != 52L * 8L) {
  stop("Expected 416 selection-tile rows; found ", nrow(tile_source), ".")
}
if (any(abs(rowSums(profile_source[, ancestry_levels]) - 1) > 2e-6)) {
  stop("At least one GIA profile does not sum to 1 within rounding tolerance.")
}

assigned_pre_levels <- c(
  "Hispanic", "Black", "EAS", "SAS", "Middle Eastern",
  "Native American", "White", "Other/Unknown"
)
selection_pre_levels <- c(
  "Hispanic", "White", "Black", "EAS", "SAS", "Middle Eastern",
  "Native American", "Other/Unknown"
)

status_levels_source <- c("Single/no multiple report", "Multiple")
status_levels_display <- c("Single/no multiple PRE", "Multiple PREs")
status_label_map <- setNames(status_levels_display, status_levels_source)

# Prepare entropy and profile data ----

entropy_plot_data <- entropy_individual %>%
  transmute(
    anonymous_plot_index,
    assigned_pre = assigned_sre,
    pre_reporting_status = unname(status_label_map[reporting_status]),
    entropy_bits
  ) %>%
  mutate(
    assigned_pre = factor(assigned_pre, levels = assigned_pre_levels),
    pre_reporting_status = factor(
      pre_reporting_status,
      levels = status_levels_display
    )
  )
if (any(is.na(entropy_plot_data$pre_reporting_status))) {
  stop("Unrecognized PRE reporting status in entropy source table.")
}

tile_source <- tile_source %>%
  mutate(present = as.logical(present))
if (any(is.na(tile_source$present))) {
  stop("Selection-tile presence values could not be parsed as logical values.")
}

# Canonicalize each reported PRE combination as an unordered set. This merges
# identical selection sets that appeared in different source-field orders.
canonical_combinations <- tile_source %>%
  filter(present) %>%
  group_by(anonymous_plot_index) %>%
  summarise(
    reported_pre_combination = paste(
      sort(unique(reported_category)),
      collapse = " + "
    ),
    .groups = "drop"
  )

profile_matrix <- as.matrix(profile_source[, ancestry_levels])
profile_entropy <- -rowSums(
  ifelse(
    profile_matrix > 0,
    profile_matrix * log2(profile_matrix),
    0
  )
)

profile <- profile_source %>%
  mutate(entropy_bits = profile_entropy) %>%
  left_join(canonical_combinations, by = "anonymous_plot_index") %>%
  mutate(
    assigned_pre = factor(assigned_sre, levels = assigned_pre_levels)
  ) %>%
  arrange(
    assigned_pre,
    reported_pre_combination,
    desc(entropy_bits),
    anonymous_plot_index
  ) %>%
  mutate(profile_index = row_number())

if (any(is.na(profile$reported_pre_combination))) {
  stop("At least one multiple-PRE profile has no reconstructed PRE combination.")
}

profile_index_map <- profile %>%
  select(anonymous_plot_index, profile_index)

profile_tiles <- tile_source %>%
  left_join(profile_index_map, by = "anonymous_plot_index") %>%
  transmute(
    profile_index,
    reported_pre = factor(reported_category, levels = rev(selection_pre_levels)),
    present
  )
if (any(is.na(profile_tiles$profile_index))) {
  stop("Selection tiles did not map completely to the profile indices.")
}

profile_long <- profile %>%
  select(profile_index, assigned_pre, all_of(ancestry_levels)) %>%
  pivot_longer(
    cols = all_of(ancestry_levels),
    names_to = "gia_component",
    values_to = "proportion"
  ) %>%
  mutate(gia_component = factor(gia_component, levels = ancestry_levels))

group_boundaries <- profile %>%
  group_by(assigned_pre, .drop = TRUE) %>%
  summarise(
    start = min(profile_index),
    end = max(profile_index),
    center = (start + end) / 2,
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(assigned_pre) %>%
  mutate(axis_label = sprintf("%s (n = %d)", assigned_pre, n))

entropy_summary_plot <- entropy_summary %>%
  transmute(
    assigned_pre = assigned_sre,
    pre_reporting_status = unname(status_label_map[reporting_status]),
    n,
    mean_entropy_bits,
    sd_entropy_bits,
    median_entropy_bits,
    iqr_entropy_bits
  )

entropy_tests <- entropy_within_pre %>%
  transmute(
    assigned_pre = assigned_sre,
    multiple_pre_n = multiple_n,
    single_or_no_multiple_pre_n = single_n,
    mean_difference_multiple_minus_single_bits,
    wilcoxon_p,
    wilcoxon_bh_adjusted_p
  )

entropy_overall_summary <- entropy_overall %>%
  transmute(
    single_or_no_multiple_pre_n = single_n,
    single_or_no_multiple_pre_mean_bits = single_mean_bits,
    single_or_no_multiple_pre_sd_bits = single_sd_bits,
    multiple_pre_n = multiple_n,
    multiple_pre_mean_bits = multiple_mean_bits,
    multiple_pre_sd_bits = multiple_sd_bits,
    mean_difference_multiple_minus_single_bits,
    mean_difference_bootstrap_ci_low,
    mean_difference_bootstrap_ci_high,
    cliffs_delta,
    cliffs_delta_bootstrap_ci_low,
    cliffs_delta_bootstrap_ci_high,
    wilcoxon_w,
    wilcoxon_p,
    bootstrap_replicates,
    seed
  )

largest_gia_counts <- profile %>%
  count(largest_gia_component = majority_ga, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(match(largest_gia_component, ancestry_levels))

assigned_pre_counts <- profile %>%
  count(assigned_pre, name = "n", .drop = TRUE) %>%
  mutate(percent = 100 * n / sum(n))

combination_counts <- profile %>%
  count(reported_pre_combination, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n), reported_pre_combination)

descriptive_patterns <- data.frame(
  metric = c(
    "EAS-assigned multiple-PRE participants with EUR >10%",
    "Hispanic-assigned multiple-PRE participants with AFR >10%",
    "Among Hispanic-assigned participants with AFR >10%, reporting Black",
    "Hispanic-assigned multiple-PRE participants with EAS >10%",
    "Among Hispanic-assigned participants with EAS >10%, reporting EAS or SAS",
    "Hispanic-assigned multiple-PRE participants with EUR >50%",
    "Among Hispanic-assigned participants with EUR >50%, reporting White"
  ),
  numerator = c(
    sum(profile$assigned_sre == "EAS" & profile$EUR > 0.10),
    sum(profile$assigned_sre == "Hispanic" & profile$AFR > 0.10),
    sum(
      profile$assigned_sre == "Hispanic" & profile$AFR > 0.10 &
        grepl("Black", profile$reported_pre_combination, fixed = TRUE)
    ),
    sum(profile$assigned_sre == "Hispanic" & profile$EAS > 0.10),
    sum(
      profile$assigned_sre == "Hispanic" & profile$EAS > 0.10 &
        grepl("EAS|SAS", profile$reported_pre_combination)
    ),
    sum(profile$assigned_sre == "Hispanic" & profile$EUR > 0.50),
    sum(
      profile$assigned_sre == "Hispanic" & profile$EUR > 0.50 &
        grepl("White", profile$reported_pre_combination, fixed = TRUE)
    )
  ),
  denominator = c(
    sum(profile$assigned_sre == "EAS"),
    sum(profile$assigned_sre == "Hispanic"),
    sum(profile$assigned_sre == "Hispanic" & profile$AFR > 0.10),
    sum(profile$assigned_sre == "Hispanic"),
    sum(profile$assigned_sre == "Hispanic" & profile$EAS > 0.10),
    sum(profile$assigned_sre == "Hispanic"),
    sum(profile$assigned_sre == "Hispanic" & profile$EUR > 0.50)
  )
) %>%
  mutate(percent = 100 * numerator / denominator)

write.csv(
  entropy_plot_data,
  file.path(table_dir, "figure_entropy_source_restricted_internal.csv"),
  row.names = FALSE
)
write.csv(
  entropy_summary_plot,
  file.path(table_dir, "entropy_by_pre_summary.csv"),
  row.names = FALSE
)
write.csv(
  entropy_overall_summary,
  file.path(table_dir, "entropy_overall_comparison.csv"),
  row.names = FALSE
)
write.csv(
  entropy_tests,
  file.path(table_dir, "entropy_within_assigned_pre_tests.csv"),
  row.names = FALSE
)
write.csv(
  profile %>%
    transmute(
      anonymous_profile_index = profile_index,
      assigned_pre,
      reported_pre_combination,
      largest_gia_component = majority_ga,
      largest_gia_proportion = majority_ga_proportion,
      entropy_bits,
      AMR, AFR, EUR, SAS, EAS
    ),
  file.path(table_dir, "figure_multiple_pre_gia_source_restricted_internal.csv"),
  row.names = FALSE
)
write.csv(
  profile_tiles %>%
    transmute(
      anonymous_profile_index = profile_index,
      reported_pre,
      present
    ),
  file.path(
    table_dir,
    "figure_multiple_pre_selection_tiles_source_restricted_internal.csv"
  ),
  row.names = FALSE
)
write.csv(
  largest_gia_counts,
  file.path(table_dir, "multiple_pre_largest_gia_counts.csv"),
  row.names = FALSE
)
write.csv(
  assigned_pre_counts,
  file.path(table_dir, "multiple_pre_assigned_pre_counts.csv"),
  row.names = FALSE
)
write.csv(
  combination_counts,
  file.path(table_dir, "multiple_pre_combination_counts.csv"),
  row.names = FALSE
)
write.csv(
  descriptive_patterns,
  file.path(table_dir, "multiple_pre_descriptive_patterns.csv"),
  row.names = FALSE
)

ancestry_palette <- c(
  AMR = "#E69F00",
  AFR = "#D55E00",
  EUR = "#0072B2",
  SAS = "#009E73",
  EAS = "#CC79A7"
)
pre_palette <- c(
  Hispanic = ancestry_palette[["AMR"]],
  White = ancestry_palette[["EUR"]],
  Black = ancestry_palette[["AFR"]],
  EAS = ancestry_palette[["EAS"]],
  SAS = ancestry_palette[["SAS"]],
  `Middle Eastern` = "#56B4E9",
  `Native American` = "#F0E442",
  `Other/Unknown` = "#8C8C8C"
)
status_palette <- c(
  `Single/no multiple PRE` = "#0072B2",
  `Multiple PREs` = "#D55E00"
)

# Figure panels ----

theme_publication <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.04), hjust = 0),
      plot.subtitle = element_text(color = "grey30", size = rel(0.90), hjust = 0),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(0.42, "cm"),
      plot.margin = margin(7, 9, 7, 7)
    )
}

overall_row <- entropy_overall_summary[1, ]
overall_x_labels <- c(
  `Single/no multiple PRE` = sprintf(
    "Single/no multiple PRE\n(n = %d)",
    overall_row$single_or_no_multiple_pre_n
  ),
  `Multiple PREs` = sprintf(
    "Multiple PREs\n(n = %d)",
    overall_row$multiple_pre_n
  )
)

panel_a <- ggplot(
  entropy_plot_data,
  aes(
    x = pre_reporting_status,
    y = entropy_bits,
    color = pre_reporting_status,
    fill = pre_reporting_status
  )
) +
  geom_violin(width = 0.78, alpha = 0.16, linewidth = 0.55, trim = TRUE) +
  geom_boxplot(
    width = 0.20,
    outlier.shape = NA,
    alpha = 0.34,
    linewidth = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 20260722),
    alpha = 0.34,
    size = 0.85,
    stroke = 0
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    fill = "white",
    color = "black",
    size = 2.4,
    stroke = 0.6
  ) +
  scale_color_manual(values = status_palette, guide = "none") +
  scale_fill_manual(values = status_palette, guide = "none") +
  scale_x_discrete(labels = overall_x_labels) +
  scale_y_continuous(
    breaks = seq(0, 2, 0.5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, log2(5))) +
  labs(
    title = "A",
    x = NULL,
    y = "Shannon entropy (bits)"
  ) +
  theme_publication(9.5)

subgroup_levels <- c(
  "Hispanic", "Black", "EAS", "Middle Eastern", "Other/Unknown"
)
subgroup_data <- entropy_plot_data %>%
  filter(as.character(assigned_pre) %in% subgroup_levels) %>%
  mutate(assigned_pre = factor(as.character(assigned_pre), levels = subgroup_levels))

subgroup_count_labels <- vapply(subgroup_levels, function(pre_name) {
  single_n <- sum(
    subgroup_data$assigned_pre == pre_name &
      subgroup_data$pre_reporting_status == "Single/no multiple PRE"
  )
  multiple_n <- sum(
    subgroup_data$assigned_pre == pre_name &
      subgroup_data$pre_reporting_status == "Multiple PREs"
  )
  sprintf("%s\nn = %d/%d", pre_name, single_n, multiple_n)
}, character(1))
names(subgroup_count_labels) <- subgroup_levels

panel_b <- ggplot(
  subgroup_data,
  aes(
    x = assigned_pre,
    y = entropy_bits,
    color = pre_reporting_status,
    fill = pre_reporting_status
  )
) +
  geom_boxplot(
    width = 0.62,
    position = position_dodge(width = 0.72),
    outlier.shape = NA,
    alpha = 0.16,
    linewidth = 0.55
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.13,
      jitter.height = 0,
      dodge.width = 0.72,
      seed = 20260722
    ),
    alpha = 0.48,
    size = 0.85,
    stroke = 0
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    fill = "white",
    size = 2.1,
    stroke = 0.55,
    position = position_dodge(width = 0.72)
  ) +
  scale_color_manual(values = status_palette, guide = "none") +
  scale_fill_manual(
    values = status_palette,
    name = "PRE reporting status",
    drop = FALSE
  ) +
  scale_x_discrete(labels = subgroup_count_labels) +
  scale_y_continuous(
    breaks = seq(0, 2, 0.5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, log2(5))) +
  labs(
    title = "B",
    x = "Assigned PRE",
    y = "Shannon entropy (bits)"
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.35))) +
  theme_publication(9.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

separator_data <- group_boundaries %>%
  slice_head(n = max(0L, nrow(group_boundaries) - 1L))

panel_c <- ggplot(
  profile_long,
  aes(x = profile_index, y = proportion, fill = gia_component)
) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(
    data = separator_data,
    aes(xintercept = end + 0.5),
    color = "white",
    linewidth = 0.65
  ) +
  scale_fill_manual(
    values = ancestry_palette,
    breaks = ancestry_levels,
    name = "GIA component",
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "C",
    x = NULL,
    y = "GIA proportion"
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme_publication(9.5) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

panel_d <- ggplot(
  profile_tiles,
  aes(
    x = profile_index,
    y = reported_pre,
    fill = ifelse(present, as.character(reported_pre), "Not selected")
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_vline(
    data = separator_data,
    aes(xintercept = end + 0.5),
    color = "grey35",
    linewidth = 0.42
  ) +
  scale_fill_manual(
    values = c(pre_palette, `Not selected` = "#F2F2F2"),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0.5, nrow(profile) + 0.5),
    breaks = group_boundaries$center,
    labels = group_boundaries$axis_label,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = NULL,
    x = "Assigned PRE",
    y = "Reported PRE selection"
  ) +
  theme_publication(9.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
  )

entropy_figure <- (panel_a | panel_b) +
  plot_layout(widths = c(0.40, 0.60), guides = "collect") &
  theme(legend.position = "right")

profile_figure <- (panel_c / panel_d) +
  plot_layout(heights = c(2.7, 1.75), guides = "collect") &
  theme(legend.position = "right")

combined_figure <- (entropy_figure / profile_figure) +
  plot_layout(heights = c(1.0, 1.55), guides = "collect") &
  theme(legend.position = "right")

entropy_paths <- save_figure_pair(
  entropy_figure,
  figure_dir,
  "figure2_entropy",
  width = 12.5,
  height = 5.0
)
profile_paths <- save_figure_pair(
  profile_figure,
  figure_dir,
  "figure2_multiple_pre_profiles",
  width = 13.5,
  height = 7.2
)
combined_paths <- save_figure_pair(
  combined_figure,
  figure_dir,
  "figure2",
  width = 13.5,
  height = 11.2
)

# Output records ----

input_checksums <- data.frame(
  input_role = names(input_paths),
  source_path = file.path("results", "descriptive", "tables", basename(input_paths)),
  sha256 = sha256_files(input_paths)
)
write.csv(
  input_checksums,
  file.path(table_dir, "input_checksums_sha256.csv"),
  row.names = FALSE
)

output_paths <- c(
  entropy_paths,
  profile_paths,
  combined_paths,
  list.files(table_dir, full.names = TRUE)
)
output_manifest <- data.frame(
  file = vapply(
    output_paths,
    function(path) relative_to_project(path, project_root),
    character(1)
  ),
  bytes = unname(file.info(output_paths)$size),
  sha256 = sha256_files(output_paths)
)
write.csv(
  output_manifest,
  file.path(run_dir, "output_manifest_sha256.csv"),
  row.names = FALSE
)

capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))

caption_text <- c(
  "**Figure X. Global genetically inferred ancestry complexity and profiles among participants with multiple parent-reported ethnicity (PRE) selections.**",
  "",
  "(A) Shannon entropy of the five global GIA proportions among participants with zero or one nonmissing PRE selection (Single/no multiple PRE) and those with at least two nonmissing PRE selections (Multiple PREs). Higher entropy indicates a more even distribution across GIA components; the theoretical maximum for five equally represented components is log2(5) = 2.322 bits. (B) Entropy distributions within assigned PRE categories represented among participants with multiple PRE selections; x-axis counts are shown as Single/no multiple PRE followed by Multiple PREs. Boxplots show the median and interquartile range, whiskers extend to 1.5 times the interquartile range, points represent participants, and white diamonds indicate means. (C) Individual global GIA profiles for the 52 participants with multiple PRE selections, with the harmonized PRE-selection matrix aligned beneath the ancestry profiles in the same anonymous order. Participants are ordered by assigned PRE, canonicalized reported-PRE combination, and entropy. PRE was assigned with the hierarchy described in the Methods, while all original nonmissing selections were retained for this analysis. Multiple-PRE status was defined before harmonization; therefore, two or more source selections that map to one displayed category may appear as a single tile. Panel C is descriptive and contains no sample identifiers."
)
writeLines(caption_text, file.path(run_dir, "figure_caption.md"))

get_summary_row <- function(pre_name, status_name) {
  entropy_summary_plot %>%
    filter(
      assigned_pre == pre_name,
      pre_reporting_status == status_name
    ) %>%
    slice(1)
}
hispanic_single <- get_summary_row("Hispanic", "Single/no multiple PRE")
hispanic_multiple <- get_summary_row("Hispanic", "Multiple PREs")
eas_single <- get_summary_row("EAS", "Single/no multiple PRE")
eas_multiple <- get_summary_row("EAS", "Multiple PREs")
black_single <- get_summary_row("Black", "Single/no multiple PRE")
black_multiple <- get_summary_row("Black", "Multiple PREs")
other_single <- get_summary_row("Other/Unknown", "Single/no multiple PRE")
other_multiple <- get_summary_row("Other/Unknown", "Multiple PREs")

test_row <- function(pre_name) {
  entropy_tests %>% filter(assigned_pre == pre_name) %>% slice(1)
}
hispanic_test <- test_row("Hispanic")
eas_test <- test_row("EAS")
black_test <- test_row("Black")
other_test <- test_row("Other/Unknown")

largest_lookup <- setNames(largest_gia_counts$n, largest_gia_counts$largest_gia_component)
combination_top <- combination_counts %>% slice_head(n = 5)

results_text <- c(
  "### GIA entropy and profiles among participants with multiple PRE selections",
  "",
  sprintf(
    paste0(
      "Participants with multiple PRE selections had substantially greater GIA entropy than those with zero or one selection (Figure XA). ",
      "Mean entropy was %.3f ± %.3f bits among the %d participants with multiple PRE selections and %.3f ± %.3f bits among the %d participants in the Single/no multiple PRE group. ",
      "The mean difference was %.3f bits (bootstrap 95%% CI, %.3f–%.3f), and Cliff's delta was %.3f (95%% CI, %.3f–%.3f; Wilcoxon P = %.3g)."
    ),
    overall_row$multiple_pre_mean_bits,
    overall_row$multiple_pre_sd_bits,
    overall_row$multiple_pre_n,
    overall_row$single_or_no_multiple_pre_mean_bits,
    overall_row$single_or_no_multiple_pre_sd_bits,
    overall_row$single_or_no_multiple_pre_n,
    overall_row$mean_difference_multiple_minus_single_bits,
    overall_row$mean_difference_bootstrap_ci_low,
    overall_row$mean_difference_bootstrap_ci_high,
    overall_row$cliffs_delta,
    overall_row$cliffs_delta_bootstrap_ci_low,
    overall_row$cliffs_delta_bootstrap_ci_high,
    overall_row$wilcoxon_p
  ),
  "",
  sprintf(
    paste0(
      "Within assigned PRE groups represented in both reporting-status strata (Figure XB), the largest difference was observed for EAS: mean entropy was %.3f ± %.3f bits for multiple PREs (n = %d) and %.3f ± %.3f bits for Single/no multiple PRE (n = %d), a difference of %.3f bits (Benjamini–Hochberg-adjusted P = %.4g). ",
      "Hispanic participants with multiple PREs also had higher entropy (%.3f ± %.3f vs %.3f ± %.3f bits; difference, %.3f bits; adjusted P = %.3g). ",
      "For Black participants, the mean difference was %.3f bits (unadjusted P = %.4f; adjusted P = %.4f), whereas Other/Unknown showed little difference (%.3f bits; adjusted P = %.3f). The Middle Eastern group included only one multiple-PRE participant and was not tested inferentially."
    ),
    eas_multiple$mean_entropy_bits,
    eas_multiple$sd_entropy_bits,
    eas_multiple$n,
    eas_single$mean_entropy_bits,
    eas_single$sd_entropy_bits,
    eas_single$n,
    eas_test$mean_difference_multiple_minus_single_bits,
    eas_test$wilcoxon_bh_adjusted_p,
    hispanic_multiple$mean_entropy_bits,
    hispanic_multiple$sd_entropy_bits,
    hispanic_single$mean_entropy_bits,
    hispanic_single$sd_entropy_bits,
    hispanic_test$mean_difference_multiple_minus_single_bits,
    hispanic_test$wilcoxon_bh_adjusted_p,
    black_test$mean_difference_multiple_minus_single_bits,
    black_test$wilcoxon_p,
    black_test$wilcoxon_bh_adjusted_p,
    other_test$mean_difference_multiple_minus_single_bits,
    other_test$wilcoxon_bh_adjusted_p
  ),
  "",
  sprintf(
    paste0(
      "The individual profiles showed substantial heterogeneity among the 52 participants with multiple PRE selections (Figure XC). ",
      "The largest GIA component was EUR for %d participants (%.1f%%), EAS for %d (%.1f%%), AMR for %d (%.1f%%), and AFR for %d (%.1f%%); none had SAS as the largest component. ",
      "The most common canonical PRE combinations were %s. ",
      "Among the nine EAS-assigned participants, seven had more than 10%% EUR ancestry. Among Hispanic-assigned participants, seven of the eight with more than 10%% AFR ancestry also reported Black, and all 12 with more than 50%% EUR ancestry also reported White. These individual-level patterns are descriptive and illustrate that multiple PRE reporting corresponds to heterogeneous, rather than uniform, GIA profiles."
    ),
    largest_lookup[["EUR"]], 100 * largest_lookup[["EUR"]] / nrow(profile),
    largest_lookup[["EAS"]], 100 * largest_lookup[["EAS"]] / nrow(profile),
    largest_lookup[["AMR"]], 100 * largest_lookup[["AMR"]] / nrow(profile),
    largest_lookup[["AFR"]], 100 * largest_lookup[["AFR"]] / nrow(profile),
    paste(
      sprintf(
        "%s (n = %d, %.1f%%)",
        combination_top$reported_pre_combination,
        combination_top$n,
        combination_top$percent
      ),
      collapse = "; "
    )
  )
)
writeLines(results_text, file.path(run_dir, "results_draft.md"))

run_manifest <- c(
  "# Figure 2 entropy and multiple-PRE GIA run",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0(
    "- Command: `Rscript ",
    "analysis/03_figure2_entropy_multiple_pre.R",
    "`"
  ),
  "- Source universe: saved descriptive source tables only; no raw workbook or sample identifiers were read.",
  "- Cohort validation: 378 entropy records; 52 participants with multiple PRE selections.",
  "- Multiple PRE definition: at least two nonmissing PRE selections.",
  "- Multiple-PRE status was defined before category harmonization; multiple source selections may collapse to one displayed category.",
  "- Entropy definition: -sum(p_i log2(p_i)) across AMR, AFR, EUR, SAS, and EAS.",
  "- Component order validated against `ancestry_component_mapping.csv` for `1000G_378.5.Q`.",
  "- Combination labels were canonicalized as unordered selection sets before aggregation and plotting.",
  "- Figure styling: only panel labels A, B, and C are shown; panel C combines the GIA profiles and PRE-selection matrix; inferential annotations were removed from panel A.",
  "- Versioning: Figure 2 outputs are written to a dedicated results directory.",
  "- Sensitivity: anonymous individual-level ancestry profiles remain restricted pending governance review."
)
writeLines(run_manifest, file.path(run_dir, "run_manifest.md"))

message("Created Figure 2 files in: ", figure_dir)
