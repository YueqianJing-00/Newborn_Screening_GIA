#!/usr/bin/env Rscript

# Figure 4 analysis of prediction changes after adding GIA to the PRE model.

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
helper_dir <- file.path(project_root, "R")
source(file.path(helper_dir, "ancestry_helpers.R"))
source(file.path(helper_dir, "pre_helpers.R"))
source(file.path(helper_dir, "statistical_helpers.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function(args) {
  defaults <- list(
    source_analysis = "mma_model",
    source_run = "main_117_top10_metabolites",
    run_name = "main_prediction_shifts",
    bootstrap = 2000L,
    seed = 20260722L,
    confidence_cutoff = 0.70,
    target_sensitivity = 0.95
  )
  for (arg in args) {
    if (!grepl("^--[A-Za-z0-9_-]+=", arg)) stop("Malformed argument: ", arg)
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", parts[[1]])
    if (!key %in% names(defaults)) stop("Unknown argument: --", parts[[1]])
    if (key %in% c("source_analysis", "source_run", "run_name")) {
      if (!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", parts[[2]])) {
        stop(key, " contains unsupported characters.")
      }
      defaults[[key]] <- parts[[2]]
    } else if (key %in% c("confidence_cutoff", "target_sensitivity")) {
      value <- suppressWarnings(as.numeric(parts[[2]]))
      if (is.na(value) || value <= 0 || value >= 1) stop(key, " must be between 0 and 1.")
      defaults[[key]] <- value
    } else {
      value <- suppressWarnings(as.integer(parts[[2]]))
      if (is.na(value) || value < 1L) stop(key, " must be a positive integer.")
      defaults[[key]] <- value
    }
  }
  defaults
}

input_dir <- normalizePath(
  Sys.getenv("HGG_DATA_DIR", file.path(project_root, "data", "raw")),
  mustWork = FALSE
)
results_dir <- normalizePath(
  Sys.getenv("HGG_RESULTS_DIR", file.path(project_root, "results")),
  mustWork = FALSE
)
analysis_dir <- file.path(results_dir, "figure4_analysis")
args <- parse_args(commandArgs(trailingOnly = TRUE))

source_analysis_dir <- file.path(results_dir, args$source_analysis)
source_run_dir <- file.path(source_analysis_dir, "runs", args$source_run)
source_table_dir <- file.path(source_run_dir, "tables")

run_dir <- file.path(analysis_dir, "runs", args$run_name)
table_dir <- file.path(run_dir, "tables")
figure_dir <- file.path(run_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cat("Figure 4 GIA case-pattern analysis\n")
cat("Start time:", format(Sys.time(), tz = "America/New_York"), "\n")
cat("Source RF run:", source_run_dir, "\n")
cat("Output run:", run_dir, "\n")

# Reconstruct the analysis cohort ----

rel_path <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(project_root, .Platform$file.sep)
  if (startsWith(path, prefix)) {
    substring(path, nchar(prefix) + 1L)
  } else {
    file.path("external", basename(path))
  }
}

input_files <- c(
  phenotype_workbook = file.path(input_dir, "Scharfelab-NBS1474samples-250207.xlsx"),
  admixture_q = file.path(input_dir, "1000G_378.5.Q"),
  admixture_fam = file.path(input_dir, "1000G_378.fam"),
  reference_metadata = file.path(input_dir, "all_phase3.psam"),
  mean_oof_predictions = file.path(source_table_dir, "mean_oof_predictions.csv"),
  repeat_oof_predictions = file.path(source_table_dir, "oof_predictions_by_repeat.csv"),
  source_modeling_dataset = file.path(
    source_table_dir,
    "modeling_dataset_deidentified_restricted_internal.csv"
  )
)

numeric_clean <- function(x) suppressWarnings(as.numeric(as.character(x)))

global_ancestry <- read_global_ancestry(
  input_files[["admixture_q"]],
  input_files[["admixture_fam"]],
  input_files[["reference_metadata"]]
)
fwrite(global_ancestry$mapping, file.path(table_dir, "global_component_mapping.csv"))

phenotype <- as.data.table(read_excel(
  input_files[["phenotype_workbook"]],
  sheet = "1474 DBS with full NBS data"
))
id_col <- "NBS-sample-ID"
if (!id_col %in% names(phenotype)) stop("Phenotype workbook is missing NBS-sample-ID.")
if (anyDuplicated(phenotype[[id_col]])) stop("Phenotype IDs are not unique.")
if (!all(global_ancestry$study$raw_id %in% phenotype[[id_col]])) stop("Study IDs do not all match.")
phenotype <- phenotype[match(global_ancestry$study$raw_id, phenotype[[id_col]])]
if (!identical(as.character(phenotype[[id_col]]), global_ancestry$study$raw_id)) {
  stop("ID-based phenotype ordering failed.")
}

pre_columns <- paste0("RACE_ETH_", 1:4)
for (column in pre_columns) {
  values <- as.character(phenotype[[column]])
  values[values %in% east_asian_pre_labels] <- "EAS"
  values[values %in% south_asian_pre_labels] <- "SAS"
  set(phenotype, j = column, value = values)
}

phenotype[, PRE_category := apply(.SD, 1, assign_pre), .SDcols = pre_columns]
phenotype[
  , n_PRE_selections := apply(.SD, 1, function(x) length(normalize_pre_values(x))),
  .SDcols = pre_columns
]
phenotype[, PRE_reporting := fifelse(n_PRE_selections > 1L, "Multiple PRE", "Single PRE")]

metabolite_features <- c(
  "ALA", "ARG", "C02", "C03", "C03DC", "C04", "C05", "C051", "C05DC", "C05OH",
  "C06", "C08", "C081", "C10", "C101", "C12", "C121", "C14", "C141", "C14OH",
  "C16", "C161", "C16OH", "C18", "C181", "C181OH", "C182", "C18OH", "CIT",
  "GLY", "MET", "ORN", "OXP", "PHE", "PRO", "TYR", "VAL", "XLE"
)
required_columns <- c(
  id_col, "Group (patient, controls, falsepos)", "BIRTH_WT", "TPN_HYPERAL",
  "AGE_AT_COLCTN", "GENDER", pre_columns, metabolite_features
)
missing_columns <- setdiff(required_columns, names(phenotype))
if (length(missing_columns) > 0L) stop("Missing phenotype columns: ", paste(missing_columns, collapse = ", "))

numeric_columns <- c("BIRTH_WT", "TPN_HYPERAL", "AGE_AT_COLCTN", metabolite_features)
phenotype[, (numeric_columns) := lapply(.SD, numeric_clean), .SDcols = numeric_columns]
phenotype[, male := fcase(
  GENDER == "M", 1,
  GENDER == "F", 0,
  default = NA_real_
)]
phenotype[, outcome := fcase(
  `Group (patient, controls, falsepos)` == "patients", "TP",
  `Group (patient, controls, falsepos)` %in% c("falsepos284", "falsepos754"), "FP",
  default = NA_character_
)]
if (anyNA(phenotype$outcome)) stop("Unexpected outcome group labels.")
phenotype[, c3_c2_ratio := fifelse(C02 > 0, C03 / C02, NA_real_)]
phenotype[, mma_evaluable := !is.na(C03) & (!is.na(c3_c2_ratio) | C03 >= 6.3)]
phenotype[, mma_screen_positive := mma_evaluable & (C03 >= 6.3 | c3_c2_ratio >= 0.3)]

joined <- cbind(phenotype, global_ancestry$study[, .(AFR, AMR, EAS, EUR, SAS)])
pre_levels <- c(
  "White", "Black", "Hispanic", "EAS", "SAS", "Middle Eastern",
  "Native American", "Other/Unknown"
)
joined[, PRE_category := factor(PRE_category, levels = pre_levels)]

core_matrix_all <- cbind(
  data.frame(
    BIRTH_WT = joined$BIRTH_WT,
    AGE_AT_COLCTN = joined$AGE_AT_COLCTN,
    male = joined$male,
    check.names = FALSE
  ),
  as.data.frame(joined[, ..metabolite_features], check.names = FALSE)
)
pre_matrix_all <- model.matrix(~ PRE_category, data = joined)[, -1, drop = FALSE]
gia_matrix_all <- as.data.frame(joined[, .(AFR, AMR, EAS, EUR)], check.names = FALSE)
full_matrix_all <- cbind(core_matrix_all, pre_matrix_all, gia_matrix_all)

eligible <- joined$mma_evaluable & joined$mma_screen_positive &
  !is.na(joined$BIRTH_WT) & joined$BIRTH_WT >= 1000 & joined$BIRTH_WT <= 5000 &
  !is.na(joined$AGE_AT_COLCTN) & joined$AGE_AT_COLCTN >= 12 & joined$AGE_AT_COLCTN <= 168 &
  !is.na(joined$TPN_HYPERAL) & joined$TPN_HYPERAL == 0 &
  complete.cases(full_matrix_all)
model_data <- joined[eligible]
if (nrow(model_data) != 117L || sum(model_data$outcome == "TP") != 85L || sum(model_data$outcome == "FP") != 32L) {
  stop("Cohort did not reproduce 117 newborns (85 TP, 32 FP).")
}

sorted_ids <- sort(model_data[[id_col]])
analysis_id_map <- setNames(sprintf("MMA_TPN0_%03d", seq_along(sorted_ids)), sorted_ids)
model_data[, analysis_id := unname(analysis_id_map[get(id_col)])]
if (anyNA(model_data$analysis_id) || anyDuplicated(model_data$analysis_id)) stop("Analysis ID creation failed.")

ancestry_components <- c("AFR", "AMR", "EAS", "EUR", "SAS")
ancestry_matrix <- as.matrix(model_data[, ..ancestry_components])
majority_index <- max.col(ancestry_matrix, ties.method = "first")
model_data[, majority_GIA := ancestry_components[majority_index]]
model_data[, majority_GIA_proportion := ancestry_matrix[cbind(seq_len(.N), majority_index)]]
model_data[, shannon_entropy_bits := apply(
  ancestry_matrix,
  1,
  function(p) -sum(p[p > 0] * log2(p[p > 0]))
)]
model_data[, normalized_entropy := shannon_entropy_bits / log2(length(ancestry_components))]
model_data[, majority_GIA_confidence := fifelse(
  majority_GIA_proportion > args$confidence_cutoff,
  paste0("Largest GIA component >", round(args$confidence_cutoff * 100), "%"),
  paste0("Largest GIA component <=", round(args$confidence_cutoff * 100), "%")
)]

pre_to_gia <- c(White = "EUR", Black = "AFR", Hispanic = "AMR", EAS = "EAS", SAS = "SAS")
model_data[, mapped_PRE_component := unname(pre_to_gia[as.character(PRE_category)])]
model_data[, PRE_GIA_relationship := fcase(
  is.na(mapped_PRE_component), "No direct five-component mapping",
  mapped_PRE_component == majority_GIA, "PRE-GIA concordant",
  default = "PRE-GIA discordant"
)]
model_data[, PRE_aligned_GIA_proportion := mapply(
  function(component, row_index) {
    if (is.na(component)) return(NA_real_)
    ancestry_matrix[row_index, component]
  },
  mapped_PRE_component,
  seq_len(.N)
)]

mean_predictions <- fread(input_files[["mean_oof_predictions"]])
repeat_predictions <- fread(input_files[["repeat_oof_predictions"]])
source_modeling_dataset <- fread(input_files[["source_modeling_dataset"]])
if (!all(c("analysis_id", "outcome", "FC", "C3_C2") %in% names(source_modeling_dataset))) {
  stop("Source modeling dataset is missing analysis_id, outcome, FC, or C3_C2.")
}
if (
  nrow(source_modeling_dataset) != 117L ||
  !setequal(source_modeling_dataset$analysis_id, model_data$analysis_id) ||
  !identical(
    source_modeling_dataset[order(analysis_id), .(analysis_id, outcome)],
    model_data[order(analysis_id), .(analysis_id, outcome)]
  )
) {
  stop("Reconstructed cohort does not match the source RF modeling dataset by analysis ID and outcome.")
}
model_names <- c(
  core = "Clinical + metabolites",
  pre = "Clinical + metabolites + PRE",
  gia = "Clinical + metabolites + GIA",
  both = "Clinical + metabolites + PRE + GIA"
)
if (!setequal(unique(mean_predictions$model), unname(model_names))) {
  stop("Unexpected model names in mean OOF predictions.")
}

mean_wide <- dcast(mean_predictions, analysis_id + outcome ~ model, value.var = "probability")
setnames(mean_wide, unname(model_names), paste0("p_", names(model_names)))
subject <- merge(
  model_data[, c(
    "analysis_id", "outcome", "PRE_category", "n_PRE_selections", "PRE_reporting",
    "AFR", "AMR", "EAS", "EUR", "SAS", "majority_GIA",
    "majority_GIA_proportion", "majority_GIA_confidence", "shannon_entropy_bits",
    "normalized_entropy", "mapped_PRE_component", "PRE_GIA_relationship",
    "PRE_aligned_GIA_proportion"
  ), with = FALSE],
  mean_wide,
  by = c("analysis_id", "outcome"),
  all = FALSE,
  sort = FALSE
)
if (nrow(subject) != 117L) stop("OOF prediction join did not reproduce 117 subjects.")

# Prediction shifts and operating points ----

subject[, outcome_numeric := as.integer(outcome == "TP")]
subject[, probability_change := p_both - p_pre]
subject[, correct_direction_shift := fifelse(outcome == "TP", probability_change, -probability_change)]
subject[, brier_improvement :=
  (outcome_numeric - p_pre)^2 - (outcome_numeric - p_both)^2]
subject[, benefit_direction := fcase(
  correct_direction_shift > 0, "Moved toward correct outcome",
  correct_direction_shift < 0, "Moved away from correct outcome",
  default = "No change"
)]

repeat_wide <- dcast(
  repeat_predictions[model %in% c(model_names[["pre"]], model_names[["both"]])],
  repeat_id + analysis_id + outcome ~ model,
  value.var = "probability"
)
setnames(
  repeat_wide,
  c(model_names[["pre"]], model_names[["both"]]),
  c("p_pre", "p_both")
)
repeat_wide[, correct_direction_shift := fifelse(
  outcome == "TP", p_both - p_pre, p_pre - p_both
)]
repeat_wide[, brier_improvement :=
  (as.integer(outcome == "TP") - p_pre)^2 - (as.integer(outcome == "TP") - p_both)^2]
subject_stability <- repeat_wide[, .(
  repeats = .N,
  mean_repeat_shift = mean(correct_direction_shift),
  median_repeat_shift = median(correct_direction_shift),
  q25_repeat_shift = quantile(correct_direction_shift, 0.25),
  q75_repeat_shift = quantile(correct_direction_shift, 0.75),
  proportion_repeats_beneficial = mean(correct_direction_shift > 0)
), by = .(analysis_id, outcome)]
subject <- merge(subject, subject_stability, by = c("analysis_id", "outcome"), all.x = TRUE, sort = FALSE)

op_pre <- operating_point(subject$outcome, subject$p_pre, args$target_sensitivity)
op_both <- operating_point(subject$outcome, subject$p_both, args$target_sensitivity)
operating_points <- rbind(
  cbind(model = model_names[["pre"]], op_pre),
  cbind(model = model_names[["both"]], op_both)
)
fwrite(operating_points, file.path(table_dir, "operating_points_mean_oof.csv"))

subject[, disposition_pre := fifelse(p_pre >= op_pre$threshold, "Refer", "Do not refer")]
subject[, disposition_both := fifelse(p_both >= op_both$threshold, "Refer", "Do not refer")]
classify_transition <- function(outcome, before, after) {
  beneficial <-
    (outcome == "TP" & before == "Do not refer" & after == "Refer") |
    (outcome == "FP" & before == "Refer" & after == "Do not refer")
  fcase(
    before == after, "No change",
    beneficial, "Beneficial change",
    default = "Harmful change"
  )
}
subject[, operational_impact := classify_transition(
  outcome, disposition_pre, disposition_both
)]
operational_counts <- subject[, .N, by = .(outcome, disposition_pre, disposition_both)]
operational_grid <- CJ(
  outcome = c("TP", "FP"),
  disposition_pre = c("Do not refer", "Refer"),
  disposition_both = c("Do not refer", "Refer"),
  unique = TRUE
)
operational_transitions <- merge(
  operational_grid,
  operational_counts,
  by = c("outcome", "disposition_pre", "disposition_both"),
  all.x = TRUE,
  sort = FALSE
)
operational_transitions[is.na(N), N := 0L]
operational_transitions[, operational_impact := classify_transition(
  outcome, disposition_pre, disposition_both
)]
operational_transitions[, total_outcome := sum(N), by = outcome]
operational_transitions[, proportion := N / total_outcome]
fwrite(operational_transitions, file.path(table_dir, "operational_reclassification_counts.csv"))

# Individual outcome, PRE, GIA, and prediction fields remain restricted.
restricted_subject_path <- file.path(
  table_dir,
  "subject_level_case_patterns_source_restricted_internal.csv"
)
fwrite(subject, restricted_subject_path)

# Stability and subgroup summaries ----

repeat_operating <- rbindlist(lapply(unique(repeat_wide$repeat_id), function(repeat_id_value) {
  repeat_data <- repeat_wide[repeat_id == repeat_id_value]
  repeat_pre <- operating_point(repeat_data$outcome, repeat_data$p_pre, args$target_sensitivity)
  repeat_both <- operating_point(repeat_data$outcome, repeat_data$p_both, args$target_sensitivity)
  data.table(
    repeat_id = repeat_id_value,
    pre_sensitivity = repeat_pre$sensitivity,
    pre_specificity = repeat_pre$specificity,
    both_sensitivity = repeat_both$sensitivity,
    both_specificity = repeat_both$specificity,
    delta_specificity = repeat_both$specificity - repeat_pre$specificity
  )
}))
fwrite(repeat_operating, file.path(table_dir, "repeat_operating_point_stability.csv"))

make_subgroup_long <- function(data) {
  rbindlist(list(
    data[, .(analysis_id, outcome, correct_direction_shift, brier_improvement,
             proportion_repeats_beneficial, dimension = "Overall", stratum = "All subjects")],
    data[, .(analysis_id, outcome, correct_direction_shift, brier_improvement,
             proportion_repeats_beneficial, dimension = "PRE reporting", stratum = PRE_reporting)],
    data[PRE_GIA_relationship != "No direct five-component mapping",
         .(analysis_id, outcome, correct_direction_shift, brier_improvement,
           proportion_repeats_beneficial, dimension = "PRE-GIA relationship", stratum = PRE_GIA_relationship)],
    data[, .(analysis_id, outcome, correct_direction_shift, brier_improvement,
             proportion_repeats_beneficial, dimension = "GIA mixture", stratum = majority_GIA_confidence)]
  ), use.names = TRUE)
}

subgroup_long <- make_subgroup_long(subject)
bootstrap_indices <- function(n, replicates) {
  index <- replicate(replicates, sample.int(n, n, replace = TRUE))
  if (n == 1L) matrix(index, nrow = 1L) else index
}

set.seed(args$seed + 1000L)
subgroup_summary <- subgroup_long[, {
  n_group <- .N
  bootstrap_index <- bootstrap_indices(n_group, args$bootstrap)
  shift_boot <- colMeans(matrix(correct_direction_shift[bootstrap_index], nrow = n_group))
  benefit_boot <- colMeans(matrix((correct_direction_shift > 0)[bootstrap_index], nrow = n_group))
  .(
    n = n_group,
    mean_correct_direction_shift = mean(correct_direction_shift),
    median_correct_direction_shift = median(correct_direction_shift),
    q25_correct_direction_shift = quantile(correct_direction_shift, 0.25),
    q75_correct_direction_shift = quantile(correct_direction_shift, 0.75),
    shift_ci_low = quantile(shift_boot, 0.025, type = 6),
    shift_ci_high = quantile(shift_boot, 0.975, type = 6),
    proportion_benefiting = mean(correct_direction_shift > 0),
    benefit_ci_low = quantile(benefit_boot, 0.025, type = 6),
    benefit_ci_high = quantile(benefit_boot, 0.975, type = 6),
    mean_brier_improvement = mean(brier_improvement),
    median_repeat_stability = median(proportion_repeats_beneficial),
    publishable_cell = n_group >= 5L
  )
}, by = .(dimension, stratum, outcome)]
fwrite(subgroup_summary, file.path(table_dir, "case_pattern_subgroup_summary.csv"))

# Aggregate mechanism check for the primary case-pattern contrast. The plotted
# quantity is the raw change in predicted TP probability, so a positive value
# means that adding GIA raised the model score regardless of observed outcome.
set.seed(args$seed + 2000L)
gia_confidence_summary <- subject[, {
  n_group <- .N
  bootstrap_index <- bootstrap_indices(n_group, args$bootstrap)
  raw_change_boot <- colMeans(matrix(probability_change[bootstrap_index], nrow = n_group))
  .(
    n = n_group,
    mean_probability_change = mean(probability_change),
    change_ci_low = quantile(raw_change_boot, 0.025, type = 6),
    change_ci_high = quantile(raw_change_boot, 0.975, type = 6),
    proportion_raw_change_positive = mean(probability_change > 0)
  )
}, by = .(majority_GIA_confidence, outcome)]

gia_confidence_composition <- subject[, .(
  n_total = .N,
  n_tp = sum(outcome == "TP"),
  n_fp = sum(outcome == "FP"),
  tp_proportion = mean(outcome == "TP")
), by = majority_GIA_confidence]
gia_confidence_summary <- merge(
  gia_confidence_summary,
  gia_confidence_composition,
  by = "majority_GIA_confidence",
  all.x = TRUE,
  sort = FALSE
)
gia_confidence_summary[, display_label := paste0(
  gsub("<=", intToUtf8(8804), majority_GIA_confidence, fixed = TRUE),
  " (n=", n_total, "; TP=", n_tp, ")"
)]
fwrite(gia_confidence_summary, file.path(table_dir, "gia_confidence_score_shift_summary.csv"))
fwrite(gia_confidence_composition, file.path(table_dir, "gia_confidence_outcome_composition.csv"))

continuous_features <- c(
  "majority_GIA_proportion", "shannon_entropy_bits", "normalized_entropy",
  "PRE_aligned_GIA_proportion"
)
continuous_associations <- rbindlist(lapply(c("All", "TP", "FP"), function(outcome_group) {
  data_group <- if (outcome_group == "All") subject else subject[outcome == outcome_group]
  rbindlist(lapply(continuous_features, function(feature) {
    keep <- !is.na(data_group[[feature]])
    if (sum(keep) < 5L || length(unique(data_group[[feature]][keep])) < 2L) {
      return(data.table(
        outcome_group = outcome_group, feature = feature, n = sum(keep),
        spearman_rho = NA_real_, p_value_exploratory = NA_real_
      ))
    }
    test <- suppressWarnings(cor.test(
      data_group[[feature]][keep],
      data_group$correct_direction_shift[keep],
      method = "spearman",
      exact = FALSE
    ))
    data.table(
      outcome_group = outcome_group,
      feature = feature,
      n = sum(keep),
      spearman_rho = unname(test$estimate),
      p_value_exploratory = test$p.value
    )
  }))
}))
continuous_associations[, note := "Exploratory, unadjusted association; not confirmatory inference"]
fwrite(continuous_associations, file.path(table_dir, "continuous_case_pattern_associations.csv"))

raw_score_associations <- rbindlist(lapply(c("TP", "FP"), function(outcome_group) {
  data_group <- subject[outcome == outcome_group]
  test <- suppressWarnings(cor.test(
    data_group$majority_GIA_proportion,
    data_group$probability_change,
    method = "spearman",
    exact = FALSE
  ))
  data.table(
    outcome = outcome_group,
    n = nrow(data_group),
    feature = "Largest GIA component proportion",
    response = "Change in predicted TP probability after adding GIA",
    spearman_rho = unname(test$estimate),
    p_value_exploratory = test$p.value,
    note = "Exploratory, unadjusted association; not confirmatory inference"
  )
}))
fwrite(raw_score_associations, file.path(table_dir, "raw_score_change_associations.csv"))

# Figure source tables and plot ----

individual_plot_source <- subject[, .(
  analysis_id, outcome, correct_direction_shift, benefit_direction,
  proportion_repeats_beneficial
)]
individual_plot_source[, outcome := factor(outcome, levels = c("TP", "FP"))]
individual_plot_source[, rank_within_outcome := frank(correct_direction_shift, ties.method = "first"), by = outcome]
fwrite(individual_plot_source, file.path(table_dir, "figure4_panel_a_source_restricted_internal.csv"))

stratum_order <- c(
  "All subjects",
  "Single PRE", "Multiple PRE",
  "PRE-GIA concordant", "PRE-GIA discordant",
  paste0("Largest GIA component >", round(args$confidence_cutoff * 100), "%"),
  paste0("Largest GIA component <=", round(args$confidence_cutoff * 100), "%")
)
subgroup_plot_source <- subgroup_summary[publishable_cell == TRUE]
subgroup_plot_source[, stratum := factor(stratum, levels = rev(stratum_order))]
fwrite(subgroup_plot_source, file.path(table_dir, "figure4_panel_b_source_aggregate.csv"))

confidence_order <- c(
  paste0("Largest GIA component >", round(args$confidence_cutoff * 100), "%"),
  paste0("Largest GIA component <=", round(args$confidence_cutoff * 100), "%")
)
confidence_display_order <- gia_confidence_summary[
  match(confidence_order, majority_GIA_confidence),
  unique(display_label)
]
gia_confidence_summary[, display_label := factor(
  display_label,
  levels = rev(confidence_display_order)
)]
fwrite(
  gia_confidence_summary,
  file.path(table_dir, "figure4_panel_b_primary_source_aggregate.csv")
)

disposition_levels <- c("Do not refer", "Refer")
operational_plot_source <- copy(operational_transitions)
operational_plot_source[, disposition_pre := factor(disposition_pre, levels = disposition_levels)]
operational_plot_source[, disposition_both := factor(disposition_both, levels = disposition_levels)]
operational_plot_source[, outcome_label := factor(
  outcome,
  levels = c("TP", "FP"),
  labels = c("True-positive newborns", "False-positive newborns")
)]
fwrite(operational_plot_source, file.path(table_dir, "figure4_panel_c_source_aggregate.csv"))

theme_publication <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_blank(),
      plot.margin = margin(7, 9, 7, 7),
      plot.tag = element_text(face = "bold", size = 13)
    )
}

direction_colors <- c(
  "Moved toward correct outcome" = "#00796B",
  "Moved away from correct outcome" = "#D95F02",
  "No change" = "#7F7F7F"
)
outcome_colors <- c(TP = "#0072B2", FP = "#CC79A7")
impact_colors <- c(
  "Beneficial change" = "#009E73",
  "Harmful change" = "#D55E00",
  "No change" = "#D9D9D9"
)

panel_a <- ggplot(
  individual_plot_source,
  aes(
    x = rank_within_outcome,
    y = 100 * correct_direction_shift,
    fill = benefit_direction
  )
) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.35) +
  geom_col(width = 0.88, linewidth = 0) +
  facet_wrap(
    ~ outcome,
    nrow = 1,
    scales = "free_x",
    labeller = as_labeller(c(TP = "True-positive newborns", FP = "False-positive newborns"))
  ) +
  scale_fill_manual(values = direction_colors) +
  labs(
    x = "Newborns ordered by individual change",
    y = "Probability shift toward correct outcome\n(percentage points)"
  ) +
  theme_publication(10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

panel_b <- ggplot(
  gia_confidence_summary,
  aes(
    x = 100 * mean_probability_change,
    y = display_label,
    color = outcome
  )
) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.45) +
  geom_errorbarh(
    aes(
      xmin = 100 * change_ci_low,
      xmax = 100 * change_ci_high
    ),
    height = 0.13,
    linewidth = 0.62,
    position = position_dodge(width = 0.36)
  ) +
  geom_point(
    size = 3.2,
    position = position_dodge(width = 0.36),
    stroke = 0.45
  ) +
  scale_color_manual(
    values = outcome_colors,
    breaks = c("TP", "FP"),
    labels = c("True positive", "False positive")
  ) +
  labs(
    x = "Mean change in predicted true-positive probability\nafter adding GIA (percentage points; bootstrap 95% CI)",
    y = NULL
  ) +
  theme_publication(10) +
  theme(
    legend.position = "bottom"
  )

panel_c <- ggplot(
  operational_plot_source,
  aes(x = disposition_both, y = disposition_pre, fill = operational_impact)
) +
  geom_tile(color = "white", linewidth = 1.1) +
  geom_text(aes(label = N), size = 5, fontface = "bold") +
  facet_wrap(~ outcome_label, nrow = 1) +
  scale_fill_manual(values = impact_colors) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 0.78) +
  labs(
    x = "Disposition with PRE + GIA model",
    y = "Disposition with PRE model"
  ) +
  theme_publication(10) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = grid::unit(1.1, "lines"),
    legend.position = "bottom"
  )

figure <- panel_a / panel_b / panel_c +
  plot_layout(heights = c(1.05, 0.72, 1.02), guides = "keep") +
  plot_annotation(tag_levels = "A")

pdf_path <- file.path(figure_dir, "Figure4_GIA_incremental_case_patterns_analysis.pdf")
png_path <- file.path(figure_dir, "Figure4_GIA_incremental_case_patterns_analysis.png")
ggsave(pdf_path, figure, width = 10.5, height = 13.0, device = cairo_pdf, bg = "white")
ggsave(png_path, figure, width = 10.5, height = 13.0, dpi = 300, bg = "white")

# Aggregate report ----

high_confidence_label <- paste0(
  "Largest GIA component >", round(args$confidence_cutoff * 100), "%"
)
low_confidence_label <- paste0(
  "Largest GIA component <=", round(args$confidence_cutoff * 100), "%"
)
high_composition <- gia_confidence_composition[majority_GIA_confidence == high_confidence_label]
low_composition <- gia_confidence_composition[majority_GIA_confidence == low_confidence_label]

summary_metrics <- data.table(
  metric = c(
    "subjects", "true_positives", "false_positives",
    "mean_correct_direction_shift", "median_correct_direction_shift",
    "proportion_benefiting", "mean_brier_improvement",
    "high_confidence_n", "high_confidence_tp_proportion",
    "low_confidence_n", "low_confidence_tp_proportion",
    "pre_specificity_at_target_sensitivity", "both_specificity_at_target_sensitivity",
    "specificity_difference", "repeat_fraction_specificity_improved",
    "repeat_fraction_specificity_unchanged", "repeat_fraction_specificity_worsened"
  ),
  value = c(
    nrow(subject), sum(subject$outcome == "TP"), sum(subject$outcome == "FP"),
    mean(subject$correct_direction_shift), median(subject$correct_direction_shift),
    mean(subject$correct_direction_shift > 0), mean(subject$brier_improvement),
    high_composition$n_total, high_composition$tp_proportion,
    low_composition$n_total, low_composition$tp_proportion,
    op_pre$specificity, op_both$specificity, op_both$specificity - op_pre$specificity,
    mean(repeat_operating$delta_specificity > 0),
    mean(repeat_operating$delta_specificity == 0),
    mean(repeat_operating$delta_specificity < 0)
  )
)
fwrite(summary_metrics, file.path(table_dir, "case_pattern_key_metrics.csv"))

cat("\nKey metrics:\n")
print(summary_metrics)
cat("\nSubgroup summary:\n")
print(subgroup_summary[order(dimension, stratum, outcome)])
cat("\nGIA-confidence score-shift summary:\n")
print(gia_confidence_summary[order(majority_GIA_confidence, outcome)])
cat("\nOperational transitions:\n")
print(operational_transitions[order(outcome, disposition_pre, disposition_both)])
cat("\nContinuous exploratory associations:\n")
print(continuous_associations)
cat("\nRaw score-change associations:\n")
print(raw_score_associations)
cat("\nFigure:\n", rel_path(pdf_path), "\n", rel_path(png_path), "\n", sep = "")
cat("End time:", format(Sys.time(), tz = "America/New_York"), "\n")
