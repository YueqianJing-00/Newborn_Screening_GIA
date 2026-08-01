#!/usr/bin/env Rscript

# PRE-GIA descriptive analysis for the full sequenced cohort.

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
helper_dir <- file.path(dirname(script_path), "..", "R")
source(file.path(helper_dir, "project_setup.R"))
source(file.path(helper_dir, "ancestry_helpers.R"))
source(file.path(helper_dir, "pre_helpers.R"))
source(file.path(helper_dir, "statistical_helpers.R"))
source(file.path(helper_dir, "plot_helpers.R"))

require_packages(c(
  "data.table", "dplyr", "ggplot2", "patchwork", "ragg",
  "readxl", "scales", "tidyr"
))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(scales)
  library(tidyr)
})

paths <- project_paths(script_path)
input_dir <- paths$data
analysis_dir <- file.path(paths$results, "descriptive")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")
make_directories(table_dir, figure_dir)

seed <- 20260710L
bootstrap_replicates <- 10000L
set.seed(seed)

input_files <- c(
  phenotype_workbook = file.path(
    input_dir, "Scharfelab-NBS1474samples-250207.xlsx"
  ),
  joint_q = file.path(input_dir, "1000G_378.5.Q"),
  joint_fam = file.path(input_dir, "1000G_378.fam"),
  reference_selection = file.path(input_dir, "sample_pure.txt"),
  reference_metadata = file.path(input_dir, "all_phase3.psam"),
  reference_selection_q = file.path(input_dir, "gwas_ld_pruned.5.Q")
)
require_files(input_files)

# Read and validate inputs ----

phenotype <- read_excel(
  input_files[["phenotype_workbook"]],
  sheet = "1474 DBS with full NBS data",
  .name_repair = "unique_quiet"
)
fam <- fread(input_files[["joint_fam"]], header = FALSE)
joint_q_raw <- as.matrix(fread(input_files[["joint_q"]], header = FALSE))
reference_selection <- fread(
  input_files[["reference_selection"]],
  header = FALSE,
  col.names = c("IID", "SuperPop")
)
reference_metadata <- fread(input_files[["reference_metadata"]])
reference_q_raw <- as.matrix(
  fread(input_files[["reference_selection_q"]], header = FALSE)
)

stopifnot(
  nrow(fam) == nrow(joint_q_raw),
  ncol(joint_q_raw) == 5L,
  nrow(reference_selection) == 2158L,
  nrow(reference_metadata) == nrow(reference_q_raw),
  ncol(reference_q_raw) == 5L,
  !anyDuplicated(reference_selection$IID),
  !anyDuplicated(phenotype[["NBS-sample-ID"]])
)

# The reference set is identified explicitly from sample_pure.txt, avoiding a
# positional "first 2,158 rows" assumption.
is_reference <- fam$V2 %in% reference_selection$IID
stopifnot(sum(is_reference) == 2158L, sum(!is_reference) == 378L)

reference_ids_in_joint_order <- as.character(fam$V2[is_reference])
stopifnot(setequal(reference_ids_in_joint_order, reference_selection$IID))

# Resolve the historical false-positive prefix mismatch by a validated key,
# not by fixed row positions.
phenotype_ids <- as.character(phenotype[["NBS-sample-ID"]])
stopifnot(!anyDuplicated(canonical_sample_id(phenotype_ids)))

study_ids_raw <- as.character(fam$V2[!is_reference])
study_match <- match(study_ids_raw, phenotype_ids)
unmatched <- is.na(study_match)
study_match[unmatched] <- match(
  canonical_sample_id(study_ids_raw[unmatched]),
  canonical_sample_id(phenotype_ids)
)
stopifnot(!anyNA(study_match), !anyDuplicated(study_match), length(study_match) == 378L)
cohort_raw <- phenotype[study_match, , drop = FALSE]

# Dataset-specific column mappings, validated below against known 1000G
# superpopulations. The supervised joint Q and unsupervised selection Q have
# different raw component orders.
joint_component_order <- c("AMR", "AFR", "EUR", "SAS", "EAS")
reference_q_component_order <- c("EUR", "AMR", "SAS", "AFR", "EAS")

joint_q <- as.data.frame(joint_q_raw)
names(joint_q) <- joint_component_order
reference_q <- as.data.frame(reference_q_raw)
names(reference_q) <- reference_q_component_order

stopifnot(
  max(abs(rowSums(joint_q) - 1)) < 1e-3,
  max(abs(rowSums(reference_q) - 1)) < 1e-3
)

reference_superpop <- reference_selection$SuperPop[
  match(reference_ids_in_joint_order, reference_selection$IID)
]
joint_reference_validation <- joint_q[is_reference, , drop = FALSE] %>%
  mutate(SuperPop = reference_superpop) %>%
  group_by(SuperPop) %>%
  summarise(across(all_of(joint_component_order), mean), .groups = "drop")

reference_q_validation <- reference_q %>%
  mutate(
    IID = reference_metadata[["#IID"]],
    SuperPop = reference_metadata$SuperPop
  ) %>%
  filter(IID %in% reference_selection$IID) %>%
  group_by(SuperPop) %>%
  summarise(across(all_of(joint_component_order), mean), .groups = "drop")

validate_component_mapping(joint_reference_validation, joint_component_order)
validate_component_mapping(reference_q_validation, joint_component_order)

component_mapping <- bind_rows(
  data.frame(
    source_file = basename(input_files[["joint_q"]]),
    raw_column = paste0("V", seq_along(joint_component_order)),
    ancestry_component = joint_component_order
  ),
  data.frame(
    source_file = basename(input_files[["reference_selection_q"]]),
    raw_column = paste0("V", seq_along(reference_q_component_order)),
    ancestry_component = reference_q_component_order
  )
)
write.csv(
  component_mapping,
  file.path(table_dir, "ancestry_component_mapping.csv"),
  row.names = FALSE
)
write.csv(
  joint_reference_validation,
  file.path(table_dir, "joint_q_component_validation_means.csv"),
  row.names = FALSE
)
write.csv(
  reference_q_validation,
  file.path(table_dir, "reference_q_component_validation_means.csv"),
  row.names = FALSE
)

# The saved 2,158-person selection is exactly reproduced by max ancestry >0.80
# in gwas_ld_pruned.5.Q. This resolves the 0.75/0.80 narrative conflict for
# this descriptive release without claiming the missing upstream command log.
selected_from_threshold <- reference_metadata[["#IID"]][
  apply(reference_q, 1, max) > 0.80
]
stopifnot(
  length(selected_from_threshold) == 2158L,
  setequal(selected_from_threshold, reference_selection$IID)
)

reference_selection_audit <- data.frame(
  source_file = basename(input_files[["reference_selection_q"]]),
  criterion = "maximum ancestry proportion > 0.80",
  total_1000g = nrow(reference_q),
  selected_n = length(selected_from_threshold),
  saved_selection_n = nrow(reference_selection),
  exact_set_match = TRUE,
  stringsAsFactors = FALSE
)
write.csv(
  reference_selection_audit,
  file.path(table_dir, "reference_selection_audit.csv"),
  row.names = FALSE
)

# Construct PRE and ancestry variables ----

race_columns <- paste0("RACE_ETH_", 1:4)
stopifnot(all(race_columns %in% names(cohort_raw)))

sre_levels <- c(
  "Hispanic", "White", "Middle Eastern", "Black", "SAS", "EAS",
  "Native American", "Other/Unknown"
)
ancestry_levels <- c("AMR", "AFR", "EUR", "SAS", "EAS")

pre_priority <- c(
  "Hispanic", "Black", "EAS", "SAS", "Middle Eastern",
  "Native American", "White"
)
assigned_sre <- apply(
  cohort_raw[, race_columns, drop = FALSE],
  1,
  assign_pre,
  hierarchy = pre_priority
)
reported_selection_count <- apply(
  cohort_raw[, race_columns, drop = FALSE],
  1,
  function(values) {
    values <- trimws(as.character(values))
    sum(!is.na(values) & nzchar(values))
  }
)
reporting_status <- ifelse(
  reported_selection_count >= 2,
  "Multiple",
  "Single/no multiple report"
)

study_q <- joint_q[!is_reference, , drop = FALSE]
majority_ga <- ancestry_levels[max.col(as.matrix(study_q), ties.method = "first")]
majority_ga_proportion <- apply(study_q, 1, max)
entropy_bits <- -rowSums(
  ifelse(as.matrix(study_q) > 0, as.matrix(study_q) * log2(as.matrix(study_q)), 0)
)

cohort <- bind_cols(
  data.frame(
    assigned_sre = assigned_sre,
    reported_selection_count = reported_selection_count,
    reporting_status = reporting_status,
    majority_ga = majority_ga,
    majority_ga_proportion = majority_ga_proportion,
    entropy_bits = entropy_bits,
    raw_group = as.character(cohort_raw[["Group (patient, controls, falsepos)"]]),
    stringsAsFactors = FALSE
  ),
  study_q
)
cohort$outcome <- ifelse(cohort$raw_group == "patients", "True positive", "False positive")
cohort$assigned_sre <- factor(cohort$assigned_sre, levels = sre_levels)
cohort$reporting_status <- factor(
  cohort$reporting_status,
  levels = c("Single/no multiple report", "Multiple")
)
cohort$majority_ga <- factor(cohort$majority_ga, levels = ancestry_levels)
stopifnot(nrow(cohort) == 378L, sum(cohort$reporting_status == "Multiple") == 52L)

# Cohort-count table: aggregate only; no sample identifiers.
make_count_section <- function(section, values, levels = NULL) {
  values <- if (is.null(levels)) values else factor(values, levels = levels)
  counts <- as.data.frame(table(values), stringsAsFactors = FALSE)
  names(counts) <- c("category", "n")
  counts %>%
    mutate(
      section = section,
      denominator = sum(n),
      percent = n / denominator
    ) %>%
    select(section, category, n, denominator, percent)
}

cohort_counts <- bind_rows(
  data.frame(
    section = "Total cohort",
    category = "All matched study participants",
    n = nrow(cohort),
    denominator = nrow(cohort),
    percent = 1
  ),
  make_count_section("Outcome", cohort$outcome),
  make_count_section("Raw cohort group", cohort$raw_group),
  make_count_section("Assigned PRE", cohort$assigned_sre, sre_levels),
  make_count_section(
    "PRE reporting status",
    cohort$reporting_status,
    c("Single/no multiple report", "Multiple")
  ),
  make_count_section(
    "Number of nonmissing PRE selections",
    as.character(cohort$reported_selection_count),
    as.character(0:3)
  ),
  make_count_section("Majority genetic ancestry", cohort$majority_ga, ancestry_levels)
)
write.csv(
  cohort_counts,
  file.path(table_dir, "cohort_counts.csv"),
  row.names = FALSE
)

reference_counts <- reference_selection %>%
  left_join(
    reference_metadata %>%
      transmute(IID = .data[["#IID"]], Population, metadata_superpop = SuperPop),
    by = "IID"
  )
stopifnot(!anyNA(reference_counts$Population), all(reference_counts$SuperPop == reference_counts$metadata_superpop))
reference_counts_table <- reference_counts %>%
  count(SuperPop, Population, name = "n") %>%
  group_by(SuperPop) %>%
  mutate(superpopulation_total = sum(n)) %>%
  ungroup()
write.csv(
  reference_counts_table,
  file.path(table_dir, "reference_population_counts.csv"),
  row.names = FALSE
)

# PRE-GIA agreement ----

cross_source <- cohort %>%
  count(assigned_sre, majority_ga, .drop = FALSE, name = "n") %>%
  group_by(assigned_sre) %>%
  mutate(row_total = sum(n), row_proportion = n / row_total) %>%
  ungroup() %>%
  group_by(majority_ga) %>%
  mutate(column_total = sum(n), column_proportion = n / column_total) %>%
  ungroup() %>%
  mutate(
    assigned_sre = factor(assigned_sre, levels = sre_levels),
    majority_ga = factor(majority_ga, levels = ancestry_levels)
  ) %>%
  arrange(assigned_sre, majority_ga)
write.csv(
  cross_source,
  file.path(table_dir, "figure2_sre_majority_ga_source.csv"),
  row.names = FALSE
)

sre_to_ga <- c(Hispanic = "AMR", White = "EUR", Black = "AFR", SAS = "SAS", EAS = "EAS")
kappa_keep <- as.character(cohort$assigned_sre) %in% names(sre_to_ga)
kappa_sre <- unname(sre_to_ga[as.character(cohort$assigned_sre[kappa_keep])])
kappa_ga <- as.character(cohort$majority_ga[kappa_keep])
kappa_result <- cohen_kappa(kappa_sre, kappa_ga, ancestry_levels)

set.seed(seed + 1L)
kappa_boot <- replicate(bootstrap_replicates, {
  index <- sample.int(length(kappa_sre), length(kappa_sre), replace = TRUE)
  cohen_kappa(kappa_sre[index], kappa_ga[index], ancestry_levels)$kappa
})
kappa_boot <- kappa_boot[is.finite(kappa_boot)]
kappa_ci <- unname(quantile(kappa_boot, c(0.025, 0.975), names = FALSE))

kappa_summary <- data.frame(
  analysis = "Five-category overall Cohen kappa",
  included_sre = "Hispanic/AMR; Black/AFR; White/EUR; SAS/SAS; EAS/EAS",
  excluded_sre = "Middle Eastern; Native American; Other/Unknown",
  n = length(kappa_sre),
  estimate = kappa_result$kappa,
  bootstrap_ci_low = kappa_ci[1],
  bootstrap_ci_high = kappa_ci[2],
  observed_agreement = kappa_result$observed,
  chance_expected_agreement = kappa_result$expected,
  bootstrap_replicates = bootstrap_replicates,
  seed = seed + 1L,
  stringsAsFactors = FALSE
)
write.csv(
  kappa_summary,
  file.path(table_dir, "cohen_kappa_overall_bootstrap.csv"),
  row.names = FALSE
)

kappa_confusion <- as.data.frame.matrix(kappa_result$table)
kappa_confusion <- cbind(PRE_mapped_component = rownames(kappa_confusion), kappa_confusion)
rownames(kappa_confusion) <- NULL
write.csv(
  kappa_confusion,
  file.path(table_dir, "cohen_kappa_five_category_confusion_matrix.csv"),
  row.names = FALSE
)

# Ancestry entropy ----

entropy_source <- cohort %>%
  transmute(
    anonymous_plot_index = row_number(),
    assigned_sre,
    reporting_status,
    entropy_bits
  )
write.csv(
  entropy_source,
  file.path(table_dir, "figure3_entropy_source_restricted_internal.csv"),
  row.names = FALSE
)

entropy_summary <- cohort %>%
  group_by(assigned_sre, reporting_status, .drop = FALSE) %>%
  summarise(
    n = n(),
    mean_entropy_bits = ifelse(n() > 0, mean(entropy_bits), NA_real_),
    sd_entropy_bits = ifelse(n() > 1, sd(entropy_bits), NA_real_),
    median_entropy_bits = ifelse(n() > 0, median(entropy_bits), NA_real_),
    iqr_entropy_bits = ifelse(n() > 0, IQR(entropy_bits), NA_real_),
    .groups = "drop"
  ) %>%
  filter(n > 0)
write.csv(
  entropy_summary,
  file.path(table_dir, "figure3_entropy_summary.csv"),
  row.names = FALSE
)

entropy_single <- cohort$entropy_bits[
  cohort$reporting_status == "Single/no multiple report"
]
entropy_multiple <- cohort$entropy_bits[cohort$reporting_status == "Multiple"]
entropy_mean_difference <- mean(entropy_multiple) - mean(entropy_single)
entropy_cliffs_delta <- cliffs_delta(entropy_multiple, entropy_single)

set.seed(seed + 2L)
entropy_boot <- replicate(bootstrap_replicates, {
  boot_multiple <- sample(entropy_multiple, length(entropy_multiple), replace = TRUE)
  boot_single <- sample(entropy_single, length(entropy_single), replace = TRUE)
  c(
    mean_difference = mean(boot_multiple) - mean(boot_single),
    cliffs_delta = cliffs_delta(boot_multiple, boot_single)
  )
})
entropy_mean_ci <- unname(quantile(entropy_boot["mean_difference", ], c(0.025, 0.975)))
entropy_delta_ci <- unname(quantile(entropy_boot["cliffs_delta", ], c(0.025, 0.975)))
entropy_wilcox <- wilcox.test(
  entropy_multiple,
  entropy_single,
  alternative = "two.sided",
  exact = FALSE
)

entropy_overall_comparison <- data.frame(
  single_status_label = "Single/no multiple report",
  single_n = length(entropy_single),
  single_mean_bits = mean(entropy_single),
  single_sd_bits = sd(entropy_single),
  multiple_n = length(entropy_multiple),
  multiple_mean_bits = mean(entropy_multiple),
  multiple_sd_bits = sd(entropy_multiple),
  mean_difference_multiple_minus_single_bits = entropy_mean_difference,
  mean_difference_bootstrap_ci_low = entropy_mean_ci[1],
  mean_difference_bootstrap_ci_high = entropy_mean_ci[2],
  cliffs_delta = entropy_cliffs_delta,
  cliffs_delta_bootstrap_ci_low = entropy_delta_ci[1],
  cliffs_delta_bootstrap_ci_high = entropy_delta_ci[2],
  wilcoxon_w = unname(entropy_wilcox$statistic),
  wilcoxon_p = entropy_wilcox$p.value,
  bootstrap_replicates = bootstrap_replicates,
  seed = seed + 2L,
  stringsAsFactors = FALSE
)
write.csv(
  entropy_overall_comparison,
  file.path(table_dir, "entropy_overall_single_vs_multiple_bootstrap.csv"),
  row.names = FALSE
)

entropy_sre_tests <- lapply(sre_levels, function(sre_name) {
  group_data <- cohort %>% filter(assigned_sre == sre_name)
  x <- group_data$entropy_bits[group_data$reporting_status == "Multiple"]
  y <- group_data$entropy_bits[
    group_data$reporting_status == "Single/no multiple report"
  ]
  data.frame(
    assigned_sre = sre_name,
    multiple_n = length(x),
    single_n = length(y),
    mean_difference_multiple_minus_single_bits = if (
      length(x) > 0 && length(y) > 0
    ) mean(x) - mean(y) else NA_real_,
    wilcoxon_p = if (length(x) >= 2 && length(y) >= 2) {
      wilcox.test(x, y, exact = FALSE)$p.value
    } else {
      NA_real_
    }
  )
}) %>% bind_rows() %>%
  mutate(wilcoxon_bh_adjusted_p = p.adjust(wilcoxon_p, method = "BH"))
write.csv(
  entropy_sre_tests,
  file.path(table_dir, "entropy_by_sre_single_vs_multiple_tests.csv"),
  row.names = FALSE
)

# Plot styling ----

ancestry_palette <- c(
  AMR = "#E69F00", # orange
  AFR = "#D55E00", # vermillion
  EUR = "#0072B2", # blue
  SAS = "#009E73", # bluish green
  EAS = "#CC79A7"  # reddish purple
)
sre_palette <- c(
  Hispanic = ancestry_palette[["AMR"]],
  White = ancestry_palette[["EUR"]],
  `Middle Eastern` = "#56B4E9",
  Black = ancestry_palette[["AFR"]],
  SAS = ancestry_palette[["SAS"]],
  EAS = ancestry_palette[["EAS"]],
  `Native American` = "#F0E442",
  `Other/Unknown` = "#8C8C8C"
)
status_palette <- c(
  `Single/no multiple report` = "#0072B2",
  Multiple = "#D55E00"
)

theme_journal <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.05), hjust = 0),
      plot.subtitle = element_text(color = "grey30", size = rel(0.9), hjust = 0),
      axis.title = element_text(face = "plain"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(0.42, "cm"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

group_boundaries <- function(data, group_column) {
  data %>%
    group_by(.data[[group_column]]) %>%
    summarise(
      start = min(plot_index),
      end = max(plot_index),
      center = (start + end) / 2,
      .groups = "drop"
    )
}

# Reference and study ancestry profiles ----

population_order <- c(
  "PEL", "MXL", "CLM", "PUR",
  "GWD", "MSL", "YRI", "ESN", "LWK", "ASW", "ACB",
  "FIN", "CEU", "GBR", "IBS", "TSI",
  "PJL", "GIH", "STU", "ITU", "BEB",
  "KHV", "CDX", "CHS", "CHB", "JPT"
)

reference_plot <- bind_cols(
  reference_metadata %>%
    transmute(IID = .data[["#IID"]], SuperPop, Population),
  reference_q
) %>%
  filter(IID %in% reference_selection$IID) %>%
  mutate(
    SuperPop = factor(SuperPop, levels = ancestry_levels),
    Population = factor(Population, levels = population_order),
    expected_component = as.character(SuperPop),
    expected_proportion = case_when(
      expected_component == "AMR" ~ AMR,
      expected_component == "AFR" ~ AFR,
      expected_component == "EUR" ~ EUR,
      expected_component == "SAS" ~ SAS,
      expected_component == "EAS" ~ EAS
    )
  ) %>%
  arrange(Population, desc(expected_proportion), desc(EUR), desc(AMR), desc(SAS), desc(AFR), desc(EAS)) %>%
  mutate(plot_index = row_number())
stopifnot(nrow(reference_plot) == 2158L)

reference_plot_source <- reference_plot %>%
  transmute(
    anonymous_plot_index = plot_index,
    SuperPop,
    Population,
    AMR,
    AFR,
    EUR,
    SAS,
    EAS
  )
write.csv(
  reference_plot_source,
  file.path(table_dir, "figure1_reference_admixture_source_restricted_internal.csv"),
  row.names = FALSE
)

expected_component_for_sre <- c(
  Hispanic = "AMR", White = "EUR", Black = "AFR", SAS = "SAS", EAS = "EAS",
  `Native American` = "AMR"
)
cohort_plot <- cohort %>%
  mutate(
    assigned_sre = factor(assigned_sre, levels = sre_levels),
    expected_component = unname(expected_component_for_sre[as.character(assigned_sre)]),
    expected_proportion = case_when(
      expected_component == "AMR" ~ AMR,
      expected_component == "AFR" ~ AFR,
      expected_component == "EUR" ~ EUR,
      expected_component == "SAS" ~ SAS,
      expected_component == "EAS" ~ EAS,
      TRUE ~ majority_ga_proportion
    )
  ) %>%
  arrange(assigned_sre, majority_ga, desc(expected_proportion), desc(majority_ga_proportion)) %>%
  mutate(plot_index = row_number())

cohort_plot_source <- cohort_plot %>%
  transmute(
    anonymous_plot_index = plot_index,
    assigned_sre,
    reporting_status,
    majority_ga,
    majority_ga_proportion,
    AMR,
    AFR,
    EUR,
    SAS,
    EAS
  )
write.csv(
  cohort_plot_source,
  file.path(table_dir, "figure1_cohort_admixture_source_restricted_internal.csv"),
  row.names = FALSE
)

reference_long <- reference_plot %>%
  select(plot_index, SuperPop, Population, all_of(ancestry_levels)) %>%
  pivot_longer(all_of(ancestry_levels), names_to = "ancestry", values_to = "proportion") %>%
  mutate(ancestry = factor(ancestry, levels = ancestry_levels))
cohort_long <- cohort_plot %>%
  select(plot_index, assigned_sre, all_of(ancestry_levels)) %>%
  pivot_longer(all_of(ancestry_levels), names_to = "ancestry", values_to = "proportion") %>%
  mutate(ancestry = factor(ancestry, levels = ancestry_levels))

reference_groups <- group_boundaries(reference_plot, "Population")
cohort_groups <- group_boundaries(cohort_plot, "assigned_sre")

figure1a <- ggplot(reference_long, aes(plot_index, proportion, fill = ancestry)) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(
    data = reference_groups[-nrow(reference_groups), , drop = FALSE],
    aes(xintercept = end + 0.5),
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = ancestry_palette, breaks = ancestry_levels, drop = FALSE) +
  scale_x_continuous(
    breaks = reference_groups$center,
    labels = as.character(reference_groups$Population),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "A  1000 Genomes homogeneous references (n = 2,158)",
    subtitle = "Selected by maximum unsupervised K=5 ancestry proportion >0.80",
    x = NULL,
    y = "Ancestry proportion",
    fill = "Ancestry component"
  ) +
  theme_journal(9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    legend.position = "bottom"
  )

figure1b <- ggplot(cohort_long, aes(plot_index, proportion, fill = ancestry)) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(
    data = cohort_groups[-nrow(cohort_groups), , drop = FALSE],
    aes(xintercept = end + 0.5),
    color = "white",
    linewidth = 0.35
  ) +
  scale_fill_manual(values = ancestry_palette, breaks = ancestry_levels, drop = FALSE) +
  scale_x_continuous(
    breaks = cohort_groups$center,
    labels = as.character(cohort_groups$assigned_sre),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "B  Screen-positive newborn cohort (n = 378)",
    subtitle = "Supervised K=5 estimates grouped by assigned PRE",
    x = "Assigned PRE",
    y = "Ancestry proportion",
    fill = "Ancestry component"
  ) +
  theme_journal(9) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 8),
    legend.position = "bottom"
  )

figure1 <- (figure1a / figure1b) +
  plot_layout(heights = c(1.05, 1), guides = "collect") &
  theme(legend.position = "bottom")
save_figure_pair(
  figure1,
  figure_dir,
  "reference_and_study_ancestry",
  width = 12,
  height = 7.4
)

# PRE by largest GIA component ----

figure2_heat <- cross_source %>%
  mutate(
    label = ifelse(n == 0, "", sprintf("%d\n%s", n, percent(row_proportion, accuracy = 1))),
    text_color = ifelse(n >= 40, "white", "black"),
    assigned_sre_plot = factor(assigned_sre, levels = rev(sre_levels))
  )

figure2a <- ggplot(
  figure2_heat,
  aes(majority_ga, assigned_sre_plot, fill = n)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = text_color), size = 3.1, lineheight = 0.9) +
  scale_color_identity() +
  scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", name = "Count") +
  labs(
    title = "A  Exact cross-classification",
    subtitle = "Cell labels show count and row percentage",
    x = "Majority genetic ancestry",
    y = "Assigned PRE"
  ) +
  coord_fixed(ratio = 0.65) +
  theme_journal(9) +
  theme(panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.5))

figure2b <- ggplot(
  cross_source,
  aes(assigned_sre, row_proportion, fill = majority_ga)
) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = ancestry_palette, breaks = ancestry_levels, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "B  Largest GIA component within PRE",
    x = "Assigned PRE",
    y = "Within-PRE proportion",
    fill = "Majority GA"
  ) +
  theme_journal(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

figure2c <- ggplot(
  cross_source,
  aes(majority_ga, column_proportion, fill = assigned_sre)
) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = sre_palette, breaks = sre_levels, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "C  PRE within largest GIA component",
    x = "Majority genetic ancestry",
    y = "Within-GA proportion",
    fill = "Assigned PRE"
  ) +
  theme_journal(9) +
  theme(legend.position = "bottom")

figure2 <- figure2a / (figure2b | figure2c) +
  plot_layout(heights = c(1.15, 1))
save_figure_pair(
  figure2,
  figure_dir,
  "pre_largest_gia_concordance",
  width = 12,
  height = 9.2
)

# Entropy by PRE reporting status ----

figure3a <- ggplot(
  cohort,
  aes(assigned_sre, entropy_bits, color = reporting_status, fill = reporting_status)
) +
  geom_boxplot(
    width = 0.62,
    position = position_dodge(width = 0.72),
    outlier.shape = NA,
    alpha = 0.15,
    linewidth = 0.55
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.72),
    alpha = 0.55,
    size = 1.05,
    stroke = 0
  ) +
  scale_color_manual(values = status_palette, drop = FALSE) +
  scale_fill_manual(values = status_palette, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 2, 0.5),
    expand = expansion(mult = c(0, 0.03))
  ) +
  coord_cartesian(ylim = c(0, log2(5))) +
  labs(
    title = "A  Ancestry entropy by assigned PRE",
    subtitle = "Entropy is calculated in bits across the five ADMIXTURE proportions",
    x = "Assigned PRE",
    y = "Shannon entropy (bits)",
    color = "PRE reporting status",
    fill = "PRE reporting status"
  ) +
  theme_journal(9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

entropy_annotation <- sprintf(
  paste0(
    "Mean difference = %.3f bits\nbootstrap 95%% CI %.3f to %.3f\n",
    "Cliff's delta = %.3f (95%% CI %.3f to %.3f)\nWilcoxon p = %s"
  ),
  entropy_mean_difference,
  entropy_mean_ci[1],
  entropy_mean_ci[2],
  entropy_cliffs_delta,
  entropy_delta_ci[1],
  entropy_delta_ci[2],
  format_p_value(entropy_wilcox$p.value)
)

figure3b <- ggplot(
  cohort,
  aes(reporting_status, entropy_bits, fill = reporting_status, color = reporting_status)
) +
  geom_violin(width = 0.78, alpha = 0.15, linewidth = 0.55, trim = FALSE) +
  geom_boxplot(width = 0.23, outlier.shape = NA, alpha = 0.35, linewidth = 0.55) +
  geom_jitter(width = 0.09, alpha = 0.42, size = 1.05, stroke = 0) +
  annotate(
    "label",
    x = 1.5,
    y = 2.18,
    label = entropy_annotation,
    hjust = 0.5,
    vjust = 1,
    size = 3.0,
    label.size = 0.25,
    fill = "white"
  ) +
  scale_fill_manual(values = status_palette, drop = FALSE) +
  scale_color_manual(values = status_palette, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 2, 0.5),
    expand = expansion(mult = c(0, 0.03))
  ) +
  coord_cartesian(ylim = c(0, log2(5))) +
  scale_x_discrete(labels = c("Single/no\nmultiple report", "Multiple")) +
  labs(
    title = "B  Overall single-vs-multiple comparison",
    x = NULL,
    y = "Shannon entropy (bits)"
  ) +
  theme_journal(9) +
  theme(legend.position = "none")

figure3 <- figure3a | figure3b
save_figure_pair(
  figure3,
  figure_dir,
  "entropy_by_pre_and_reporting_status",
  width = 12,
  height = 5.9
)

# Multiple-PRE ancestry profiles ----

multi_indices <- which(cohort$reporting_status == "Multiple")
multi_plot <- cohort[multi_indices, , drop = FALSE]
multi_raw_race <- cohort_raw[multi_indices, race_columns, drop = FALSE]
multi_combinations <- apply(multi_raw_race, 1, function(values) {
  categories <- unique(na.omit(vapply(values, harmonize_pre_value, character(1))))
  paste(categories, collapse = " + ")
})
multi_plot$reported_sre_combination <- multi_combinations
multi_plot <- multi_plot %>%
  mutate(
    assigned_sre = factor(assigned_sre, levels = sre_levels),
    expected_component = unname(expected_component_for_sre[as.character(assigned_sre)]),
    expected_proportion = case_when(
      expected_component == "AMR" ~ AMR,
      expected_component == "AFR" ~ AFR,
      expected_component == "EUR" ~ EUR,
      expected_component == "SAS" ~ SAS,
      expected_component == "EAS" ~ EAS,
      TRUE ~ majority_ga_proportion
    )
  ) %>%
  arrange(assigned_sre, reported_sre_combination, majority_ga, desc(expected_proportion)) %>%
  mutate(plot_index = row_number())
stopifnot(nrow(multi_plot) == 52L)

write.csv(
  multi_plot %>%
    transmute(
      anonymous_plot_index = plot_index,
      assigned_sre,
      reported_sre_combination,
      majority_ga,
      majority_ga_proportion,
      AMR,
      AFR,
      EUR,
      SAS,
      EAS
    ),
  file.path(table_dir, "figure4_multisre_admixture_source_restricted_internal.csv"),
  row.names = FALSE
)

multi_raw_categories <- lapply(seq_len(nrow(multi_raw_race)), function(i) {
  unique(na.omit(vapply(multi_raw_race[i, ], harmonize_pre_value, character(1))))
})

# A temporary row key preserves repeated PRE combinations while sorting.
multi_with_key <- cohort[multi_indices, , drop = FALSE]
multi_with_key$temp_row_key <- seq_len(nrow(multi_with_key))
multi_with_key$reported_sre_combination <- multi_combinations
multi_with_key <- multi_with_key %>%
  mutate(
    assigned_sre = factor(assigned_sre, levels = sre_levels),
    expected_component = unname(expected_component_for_sre[as.character(assigned_sre)]),
    expected_proportion = case_when(
      expected_component == "AMR" ~ AMR,
      expected_component == "AFR" ~ AFR,
      expected_component == "EUR" ~ EUR,
      expected_component == "SAS" ~ SAS,
      expected_component == "EAS" ~ EAS,
      TRUE ~ majority_ga_proportion
    )
  ) %>%
  arrange(assigned_sre, reported_sre_combination, majority_ga, desc(expected_proportion), temp_row_key) %>%
  mutate(plot_index = row_number())
stopifnot(
  identical(
    as.character(multi_plot$reported_sre_combination),
    as.character(multi_with_key$reported_sre_combination)
  )
)

selection_tile_source <- expand_grid(
  plot_index = multi_with_key$plot_index,
  reported_category = sre_levels
) %>%
  left_join(
    multi_with_key %>% select(plot_index, temp_row_key),
    by = "plot_index"
  ) %>%
  rowwise() %>%
  mutate(present = reported_category %in% multi_raw_categories[[temp_row_key]]) %>%
  ungroup() %>%
  mutate(
    reported_category = factor(reported_category, levels = rev(sre_levels)),
    tile_fill = ifelse(
      present,
      sre_palette[as.character(reported_category)],
      "#F1F1F1"
    )
  )

write.csv(
  selection_tile_source %>%
    transmute(
      anonymous_plot_index = plot_index,
      reported_category,
      present
    ),
  file.path(table_dir, "figure4_multisre_selection_tiles_source_restricted_internal.csv"),
  row.names = FALSE
)

multi_long <- multi_with_key %>%
  select(plot_index, assigned_sre, all_of(ancestry_levels)) %>%
  pivot_longer(all_of(ancestry_levels), names_to = "ancestry", values_to = "proportion") %>%
  mutate(ancestry = factor(ancestry, levels = ancestry_levels))
multi_groups <- group_boundaries(multi_with_key, "assigned_sre")

figure4a <- ggplot(multi_long, aes(plot_index, proportion, fill = ancestry)) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(
    data = multi_groups[-nrow(multi_groups), , drop = FALSE],
    aes(xintercept = end + 0.5),
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = ancestry_palette, breaks = ancestry_levels, drop = FALSE) +
  scale_x_continuous(
    breaks = multi_groups$center,
    labels = as.character(multi_groups$assigned_sre),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "A  Genetic ancestry profiles among participants with multiple PRE selections (n = 52)",
    x = "Assigned PRE",
    y = "Ancestry proportion",
    fill = "Ancestry component"
  ) +
  theme_journal(9) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

figure4b <- ggplot(
  selection_tile_source,
  aes(plot_index, reported_category, fill = tile_fill)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_vline(
    data = multi_groups[-nrow(multi_groups), , drop = FALSE],
    aes(xintercept = end + 0.5),
    color = "grey35",
    linewidth = 0.35
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0.5, nrow(multi_with_key) + 0.5),
    expand = expansion(mult = c(0, 0)),
    breaks = NULL
  ) +
  labs(
    title = "B  Harmonized reported PRE selections in the same anonymous order",
    x = NULL,
    y = "Reported selection"
  ) +
  theme_journal(9) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
  )

figure4 <- figure4a / figure4b + plot_layout(heights = c(2.8, 1.55))
save_figure_pair(
  figure4,
  figure_dir,
  "multiple_pre_ancestry_profiles",
  width = 12,
  height = 7.8
)

# Privacy guard: no output table may contain source sample identifiers.
output_tables <- list.files(table_dir, pattern = "\\.csv$", full.names = TRUE)
for (table_path in output_tables) {
  table_text <- paste(readLines(table_path, warn = FALSE), collapse = "\n")
  if (grepl("NBSfalsepos_|NBSpatients_|HG[0-9]{4,}", table_text)) {
    stop("Privacy guard failed; sample identifier found in output table: ", table_path)
  }
}

message("Descriptive analysis complete.")
message("Figures: ", figure_dir)
