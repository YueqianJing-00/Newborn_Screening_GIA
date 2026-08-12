#!/usr/bin/env Rscript

# PRE-GIA descriptive analysis for the full sequenced cohort.

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
helper_dir <- file.path(project_root, "R")
source(file.path(helper_dir, "ancestry_helpers.R"))
source(file.path(helper_dir, "pre_helpers.R"))
source(file.path(helper_dir, "statistical_helpers.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(tidyr)
})

input_dir <- normalizePath(
  Sys.getenv("HGG_DATA_DIR", file.path(project_root, "data", "raw")),
  mustWork = FALSE
)
results_dir <- normalizePath(
  Sys.getenv("HGG_RESULTS_DIR", file.path(project_root, "results")),
  mustWork = FALSE
)
analysis_dir <- file.path(results_dir, "descriptive")
table_dir <- file.path(analysis_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

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
  reference_selection_q = file.path(input_dir, "1000G_impact.5.Q")
)

# Read inputs ----

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

# Separate reference and study rows using the saved reference IDs.
is_reference <- fam$V2 %in% reference_selection$IID
reference_ids_in_joint_order <- as.character(fam$V2[is_reference])

# Match every study row to the phenotype workbook by sample ID.
phenotype_ids <- as.character(phenotype[["NBS-sample-ID"]])
study_ids_raw <- as.character(fam$V2[!is_reference])
study_match <- match(study_ids_raw, phenotype_ids)
unmatched <- is.na(study_match)
study_match[unmatched] <- match(
  canonical_sample_id(study_ids_raw[unmatched]),
  canonical_sample_id(phenotype_ids)
)
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

# Confirm the hard-coded Q-column labels against known reference populations.
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

# Reproduce the saved reference selection from the MSK-IMPACT ancestry matrix.
reference_threshold <- 0.80
selected_from_threshold <- reference_metadata[["#IID"]][
  apply(reference_q, 1, max) > reference_threshold
]
stopifnot(
  setequal(selected_from_threshold, reference_selection$IID)
)

reference_selection_audit <- data.frame(
  source_file = basename(input_files[["reference_selection_q"]]),
  criterion = sprintf("maximum ancestry proportion > %.2f", reference_threshold),
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

# Save the threshold-specific count tables used by the two concordance plots.
cross_count_source <- cross_source %>%
  transmute(
    assigned_pre = assigned_sre,
    largest_gia_component = majority_ga,
    n
  )
cross_gt70_count_source <- cohort %>%
  filter(majority_ga_proportion > 0.70) %>%
  count(
    assigned_pre = assigned_sre,
    largest_gia_component = majority_ga,
    name = "n"
  )
write.csv(
  cross_count_source,
  file.path(table_dir, "figure1_cross_classification_all.csv"),
  row.names = FALSE
)
write.csv(
  cross_gt70_count_source,
  file.path(table_dir, "figure1_cross_classification_gt70.csv"),
  row.names = FALSE
)

# Save the continuous GIA values corresponding to directly mapped PRE categories.
pre_to_gia <- c(Hispanic = "AMR", White = "EUR", Black = "AFR", SAS = "SAS", EAS = "EAS")
corresponding_gia_source <- cohort %>%
  filter(as.character(assigned_sre) %in% names(pre_to_gia)) %>%
  mutate(
    mapped_component = unname(pre_to_gia[as.character(assigned_sre)]),
    corresponding_gia_proportion = case_when(
      mapped_component == "AMR" ~ AMR,
      mapped_component == "AFR" ~ AFR,
      mapped_component == "EUR" ~ EUR,
      mapped_component == "SAS" ~ SAS,
      mapped_component == "EAS" ~ EAS
    ),
    display_group = paste0(as.character(assigned_sre), " → ", mapped_component)
  ) %>%
  transmute(
    anonymous_plot_index = row_number(),
    display_group,
    mapped_component,
    corresponding_gia_proportion
  )
write.csv(
  corresponding_gia_source,
  file.path(table_dir, "figure1_corresponding_gia_source_restricted_internal.csv"),
  row.names = FALSE
)

sre_to_ga <- c(Hispanic = "AMR", White = "EUR", Black = "AFR", SAS = "SAS", EAS = "EAS")
# Calculate agreement after progressively restricting to concentrated GIA profiles.
threshold_sensitivity <- bind_rows(lapply(seq(0, 0.90, 0.05), function(threshold) {
  subset <- cohort %>%
    filter(
      as.character(assigned_sre) %in% names(sre_to_ga),
      majority_ga_proportion > threshold
    )
  metrics <- cohen_kappa(
    unname(sre_to_ga[as.character(subset$assigned_sre)]),
    as.character(subset$majority_ga),
    ancestry_levels
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

# Plot annotations read these estimates instead of recomputing statistics.
write.csv(
  data.frame(
    mean_difference_bits = entropy_mean_difference,
    mean_difference_ci_low = entropy_mean_ci[1],
    mean_difference_ci_high = entropy_mean_ci[2],
    cliffs_delta = entropy_cliffs_delta,
    cliffs_delta_ci_low = entropy_delta_ci[1],
    cliffs_delta_ci_high = entropy_delta_ci[2],
    wilcoxon_p = entropy_wilcox$p.value
  ),
  file.path(table_dir, "entropy_plot_annotation.csv"),
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

# Figure source tables ----

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

# Multiple-PRE ancestry profiles ----

multi_indices <- which(cohort$reporting_status == "Multiple")
multi_raw_race <- cohort_raw[multi_indices, race_columns, drop = FALSE]
multi_combinations <- apply(multi_raw_race, 1, function(values) {
  categories <- unique(na.omit(vapply(values, harmonize_pre_value, character(1))))
  paste(categories, collapse = " + ")
})
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

# Save ancestry profiles and PRE-selection tiles in one synchronized order.
write.csv(
  multi_with_key %>%
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
  ungroup()

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

# Privacy guard: no output table may contain source sample identifiers.
output_tables <- list.files(table_dir, pattern = "\\.csv$", full.names = TRUE)
for (table_path in output_tables) {
  table_text <- paste(readLines(table_path, warn = FALSE), collapse = "\n")
  if (grepl("NBSfalsepos_|NBSpatients_|HG[0-9]{4,}", table_text)) {
    stop("Privacy guard failed; sample identifier found in output table: ", table_path)
  }
}

message("Descriptive analysis complete.")
message("Tables: ", table_dir)
