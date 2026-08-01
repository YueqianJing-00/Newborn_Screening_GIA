#!/usr/bin/env Rscript

# Independent, read-only checks of the random-forest outputs.

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "..", "R", "project_setup.R"))
require_packages(c("data.table", "pROC"))

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

paths <- project_paths(script_path)
project_root <- paths$root
analysis_dir <- file.path(paths$results, "mma_model")

cli_args <- commandArgs(trailingOnly = TRUE)
run_arg <- grep("^--run-name=", cli_args, value = TRUE)
if (length(run_arg) > 1L || length(setdiff(cli_args, run_arg)) > 0L) {
  stop("Usage: Rscript analysis/05_validate_mma_random_forest.R [--run-name=NAME]")
}
run_name <- if (length(run_arg) == 1L) sub("^--run-name=", "", run_arg) else "main_117_top10_metabolites"
if (!is_run_name(run_name)) stop("Invalid run name.")

run_dir <- file.path(analysis_dir, "runs", run_name)
table_dir <- file.path(run_dir, "tables")
qa_dir <- file.path(run_dir, "qa")
make_directories(qa_dir)

checks <- character()
record_check <- function(label) {
  checks <<- c(checks, paste0("PASS: ", label))
  invisible(TRUE)
}
assert_true <- function(condition, label) {
  if (!isTRUE(condition)) stop("VALIDATION FAILED: ", label, call. = FALSE)
  record_check(label)
}
assert_near <- function(observed, expected, label, tolerance = 1e-12) {
  ok <- length(observed) == length(expected) &&
    all(is.finite(observed)) && all(is.finite(expected)) &&
    max(abs(observed - expected)) <= tolerance
  assert_true(ok, label)
}

required_tables <- c(
  "cohort_flow.csv", "cohort_summary.csv",
  "model_specification.csv", "fold_assignments.csv",
  "metabolite_selection_by_fold.csv", "metabolite_selection_summary.csv",
  "forest_fit_audit.csv", "oof_predictions_by_repeat.csv",
  "mean_oof_predictions.csv", "repeat_metrics.csv",
  "primary_operating_points.csv", "primary_performance.csv",
  "subject_bootstrap_metrics.csv", "subject_bootstrap_weights.csv",
  "paired_effects.csv", "permutation_importance_by_repeat.csv",
  "permutation_importance_summary.csv", "figure_panel_a_source.csv",
  "figure_panel_b_source.csv", "modeling_dataset_deidentified_restricted_internal.csv"
)
required_paths <- file.path(table_dir, required_tables)
assert_true(all(file.exists(required_paths)), "all required output tables exist")

read_table <- function(name) fread(file.path(table_dir, name))

# Load model outputs ----

cohort_flow <- read_table("cohort_flow.csv")
cohort <- read_table("cohort_summary.csv")
model_spec <- read_table("model_specification.csv")
folds <- read_table("fold_assignments.csv")
selection <- read_table("metabolite_selection_by_fold.csv")
selection_summary <- read_table("metabolite_selection_summary.csv")
fit_audit <- read_table("forest_fit_audit.csv")
oof <- read_table("oof_predictions_by_repeat.csv")
mean_oof <- read_table("mean_oof_predictions.csv")
repeat_metrics <- read_table("repeat_metrics.csv")
operating_points <- read_table("primary_operating_points.csv")
primary <- read_table("primary_performance.csv")
bootstrap_metrics <- read_table("subject_bootstrap_metrics.csv")
bootstrap_weights <- read_table("subject_bootstrap_weights.csv")
paired <- read_table("paired_effects.csv")
importance <- read_table("permutation_importance_by_repeat.csv")
importance_summary <- read_table("permutation_importance_summary.csv")
panel_a <- read_table("figure_panel_a_source.csv")
panel_b <- read_table("figure_panel_b_source.csv")
modeling_data <- read_table("modeling_dataset_deidentified_restricted_internal.csv")

model_order <- c(
  "Clinical + metabolites",
  "Clinical + metabolites + PRE",
  "Clinical + metabolites + GIA",
  "Clinical + metabolites + PRE + GIA"
)

# Cohort and specification ----
assert_true(
  nrow(cohort) == 1L && cohort$n == 117L && cohort$tp == 85L &&
    cohort$fp == 32L && cohort$tpn_0 == 117L && cohort$tpn_1 == 0L,
  "cohort is exactly 117 newborns after excluding those receiving TPN (85 TP, 32 FP)"
)
assert_true(
  cohort$metabolite_candidates == 40L &&
    cohort$metabolites_selected_per_fold == 10L &&
    cohort$fc_candidate && cohort$c3_c2_candidate,
  "candidate pool is 40 and includes eligible FC and C3/C2 candidates"
)
assert_true(
  identical(cohort_flow$n, c(378L, 378L, 165L, 162L, 161L, 117L, 117L)) &&
    identical(cohort_flow$tp, c(235L, 235L, 98L, 98L, 98L, 85L, 85L)) &&
    identical(cohort_flow$fp, c(143L, 143L, 67L, 64L, 63L, 32L, 32L)),
  "cohort-flow counts match the analysis filters"
)
candidate_spec <- unique(model_spec[feature_group == "Metabolite candidate", feature])
assert_true(
  length(candidate_spec) == 40L && all(c("FC", "C3_C2") %in% candidate_spec) &&
    all(model_spec[feature %in% c("FC", "C3_C2"), role] == "fold-wise candidate"),
  "FC and C3/C2 are candidates rather than forced predictors"
)

# Fold ledger ----
assert_true(
  nrow(folds) == 117L * 100L && uniqueN(folds$repeat_id) == 100L &&
    uniqueN(folds$analysis_id) == 117L,
  "fold ledger has one row per subject in each of 100 repeats"
)
assert_true(
  all(folds[, .N, by = .(repeat_id, analysis_id)]$N == 1L),
  "every subject occurs in exactly one held-out fold per repeat"
)
fold_class_check <- folds[, .(
  outcome_classes = uniqueN(outcome),
  n = .N
), by = .(repeat_id, fold)]
assert_true(
  nrow(fold_class_check) == 1000L &&
    all(fold_class_check$outcome_classes == 2L) &&
    all(fold_class_check$n %in% c(11L, 12L, 13L)),
  "all 1,000 held-out folds contain both outcomes and 11-13 subjects"
)

# Fold-wise feature selection ----
assert_true(
  nrow(selection) == 40L * 100L * 10L &&
    uniqueN(selection$repeat_id) == 100L && uniqueN(selection$fold) == 10L,
  "selection ledger has all 40 candidates for every repeat/fold"
)
selection_check <- selection[, .(
  candidates = .N,
  selected_n = sum(selected),
  selected_unique = uniqueN(metabolite[selected]),
  ranks_unique = uniqueN(rank),
  ranking_consistent = all(selected == (rank <= 10L)),
  distance_consistent = max(abs(absolute_auc_distance - abs(univariate_auc - 0.5))) < 1e-12
), by = .(repeat_id, fold)]
assert_true(
  nrow(selection_check) == 1000L &&
    all(selection_check$candidates == 40L) &&
    all(selection_check$selected_n == 10L) &&
    all(selection_check$selected_unique == 10L) &&
    all(selection_check$ranks_unique == 40L) &&
    all(selection_check$ranking_consistent) &&
    all(selection_check$distance_consistent),
  "each training fold ranks 40 candidates and selects exactly its top 10"
)
assert_true(
  nrow(selection_summary) == 40L &&
    selection_summary[metabolite == "FC", folds_selected] == 1000L &&
    selection_summary[metabolite == "C3_C2", folds_selected] == 1000L,
  "FC and C3/C2 were naturally selected in all 1,000 folds"
)
assert_true(
  nrow(fit_audit) == 4000L &&
    all(fit_audit$requested_ntree == 1000L) &&
    all(fit_audit$fitted_ntree == 1000L) &&
    all(fit_audit$requested_mtry == 7L) &&
    all(fit_audit$fitted_mtry == 7L),
  "all 4,000 forests used 1,000 trees and mtry=7"
)
fit_panel_check <- fit_audit[, .(
  models = uniqueN(model),
  selected_strings = uniqueN(selected_metabolites)
), by = .(repeat_id, fold)]
assert_true(
  nrow(fit_panel_check) == 1000L && all(fit_panel_check$models == 4L) &&
    all(fit_panel_check$selected_strings == 1L),
  "the identical selected panel was reused across all four models per split"
)
expected_p <- c(
  "Clinical + metabolites" = 13L,
  "Clinical + metabolites + PRE" = 19L,
  "Clinical + metabolites + GIA" = 17L,
  "Clinical + metabolites + PRE + GIA" = 23L
)
assert_true(
  all(fit_audit$predictors == expected_p[fit_audit$model]),
  "all four nested model matrices have the expected predictor counts"
)

# OOF predictions and subject-level averaging ----
assert_true(
  nrow(oof) == 117L * 100L * 4L && setequal(unique(oof$model), model_order) &&
    all(is.finite(oof$probability)) && all(oof$probability >= 0 & oof$probability <= 1),
  "raw OOF ledger has 46,800 valid probabilities"
)
assert_true(
  all(oof[, .N, by = .(analysis_id, model)]$N == 100L),
  "every subject/model has exactly 100 OOF predictions"
)
mean_recomputed <- oof[, .(
  outcome = outcome[1L],
  probability_recomputed = mean(probability)
), by = .(analysis_id, model)]
mean_join <- merge(
  mean_oof, mean_recomputed,
  by = c("analysis_id", "model", "outcome"), all = TRUE
)
assert_true(nrow(mean_join) == 117L * 4L, "mean OOF ledger has 468 subject/model rows")
assert_near(
  mean_join$probability, mean_join$probability_recomputed,
  "saved subject-mean probabilities reproduce exactly from raw OOF predictions"
)

discrete_operating_point <- function(outcome, probability, target = 0.95) {
  thresholds <- sort(unique(probability), decreasing = TRUE)
  candidates <- rbindlist(lapply(thresholds, function(threshold) {
    positive <- probability >= threshold
    data.table(
      threshold = threshold,
      tp = sum(positive & outcome == "TP"),
      fn = sum(!positive & outcome == "TP"),
      tn = sum(!positive & outcome == "FP"),
      fp = sum(positive & outcome == "FP")
    )
  }))
  candidates[, `:=`(
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp)
  )]
  candidates <- candidates[sensitivity >= target]
  setorder(candidates, -specificity, -threshold)
  candidates[1L]
}
metric_bundle <- function(outcome, probability) {
  outcome_factor <- factor(outcome, levels = c("FP", "TP"))
  op <- discrete_operating_point(as.character(outcome_factor), probability)
  list(
    AUC = as.numeric(pROC::auc(pROC::roc(
      outcome_factor, probability, levels = c("FP", "TP"),
      direction = "<", quiet = TRUE
    ))),
    `Specificity at >=95% sensitivity` = op$specificity,
    `Brier score` = mean((probability - as.numeric(outcome_factor == "TP"))^2),
    operating_point = op
  )
}

recomputed_rows <- lapply(model_order, function(model_name) {
  d <- mean_oof[model == model_name]
  x <- metric_bundle(d$outcome, d$probability)
  list(
    performance = data.table(
      model = model_name,
      metric = c("AUC", "Specificity at >=95% sensitivity", "Brier score"),
      recomputed = c(x$AUC, x[["Specificity at >=95% sensitivity"]], x[["Brier score"]])
    ),
    operating = data.table(
      model = model_name,
      threshold = x$operating_point$threshold,
      achieved_sensitivity = x$operating_point$sensitivity,
      specificity = x$operating_point$specificity,
      tp = x$operating_point$tp,
      fn = x$operating_point$fn,
      tn = x$operating_point$tn,
      fp = x$operating_point$fp
    )
  )
})
performance_recomputed <- rbindlist(lapply(recomputed_rows, `[[`, "performance"))
operating_recomputed <- rbindlist(lapply(recomputed_rows, `[[`, "operating"))
primary_join <- merge(primary, performance_recomputed, by = c("model", "metric"))
assert_near(
  primary_join$estimate, primary_join$recomputed,
  "AUC, specificity, and Brier estimates reproduce from subject-mean OOF predictions"
)
op_join <- merge(
  operating_points, operating_recomputed,
  by = "model", suffixes = c("_saved", "_recomputed")
)
op_numeric <- c("threshold", "achieved_sensitivity", "specificity", "tp", "fn", "tn", "fp")
assert_true(
  all(vapply(op_numeric, function(column) {
    max(abs(op_join[[paste0(column, "_saved")]] - op_join[[paste0(column, "_recomputed")]])) < 1e-12
  }, logical(1))) && all(operating_points$achieved_sensitivity >= 0.95),
  "saved empirical operating points reproduce and attain sensitivity >=0.95"
)

# Bootstrap confidence intervals ----
assert_true(
  nrow(bootstrap_metrics) == 2000L * 4L &&
    uniqueN(bootstrap_metrics$bootstrap_id) == 2000L &&
    all(bootstrap_metrics$achieved_sensitivity >= 0.95),
  "bootstrap metric ledger contains 2,000 paired replicates for four models"
)
assert_true(
  nrow(bootstrap_weights) == 2000L * 117L &&
    all(bootstrap_weights[, .N, by = bootstrap_id]$N == 117L),
  "bootstrap weight ledger contains all subjects in every replicate"
)
bootstrap_class_check <- bootstrap_weights[, .(
  sampled_tp = sum(multiplicity[outcome == "TP"]),
  sampled_fp = sum(multiplicity[outcome == "FP"]),
  sampled_total = sum(multiplicity)
), by = bootstrap_id]
assert_true(
  all(bootstrap_class_check$sampled_tp == 85L) &&
    all(bootstrap_class_check$sampled_fp == 32L) &&
    all(bootstrap_class_check$sampled_total == 117L),
  "every bootstrap replicate preserves 85 TP and 32 FP draws"
)
ci <- function(x) unname(quantile(x, c(0.025, 0.975), type = 6))
boot_long <- melt(
  bootstrap_metrics,
  id.vars = c("bootstrap_id", "model"),
  measure.vars = c("auc", "specificity_at_95_sensitivity", "brier_score"),
  variable.name = "metric_code", value.name = "value"
)
boot_long[, metric := c(
  auc = "AUC",
  specificity_at_95_sensitivity = "Specificity at >=95% sensitivity",
  brier_score = "Brier score"
)[as.character(metric_code)]]
boot_ci <- boot_long[, .(
  ci_low_recomputed = ci(value)[1L],
  ci_high_recomputed = ci(value)[2L]
), by = .(model, metric)]
primary_ci_join <- merge(primary, boot_ci, by = c("model", "metric"))
assert_near(
  primary_ci_join$ci_low, primary_ci_join$ci_low_recomputed,
  "all lower confidence limits reproduce from the bootstrap ledger"
)
assert_near(
  primary_ci_join$ci_high, primary_ci_join$ci_high_recomputed,
  "all upper confidence limits reproduce from the bootstrap ledger"
)

point_wide <- dcast(primary[, .(model, metric, estimate)], metric ~ model, value.var = "estimate")
boot_wide <- dcast(boot_long, bootstrap_id + metric ~ model, value.var = "value")
paired_check <- paired[, {
  point_row <- point_wide[metric == .BY$metric]
  boot_rows <- boot_wide[metric == .BY$metric]
  delta <- boot_rows[[augmented_model[1L]]] - boot_rows[[reference_model[1L]]]
  list(
    estimate_recomputed = point_row[[augmented_model[1L]]] - point_row[[reference_model[1L]]],
    ci_low_recomputed = ci(delta)[1L],
    ci_high_recomputed = ci(delta)[2L]
  )
}, by = .(contrast, metric, augmented_model, reference_model)]
paired_join <- merge(
  paired, paired_check,
  by = c("contrast", "metric", "augmented_model", "reference_model")
)
assert_near(paired_join$estimate, paired_join$estimate_recomputed, "paired point differences reproduce")
assert_near(paired_join$ci_low, paired_join$ci_low_recomputed, "paired lower limits reproduce")
assert_near(paired_join$ci_high, paired_join$ci_high_recomputed, "paired upper limits reproduce")

# Held-out permutation importance ----
assert_true(
  nrow(importance) == 45L * 100L && uniqueN(importance$group_id) == 45L &&
    all(importance$model == model_order[4L]),
  "importance ledger has 45 predictor groups across 100 repeats in the full model"
)
assert_near(
  importance$importance_delta_brier,
  importance$permuted_brier - importance$baseline_brier,
  "all importance values equal permuted minus baseline held-out Brier score"
)
importance_recomputed <- importance[, .(
  mean_recomputed = mean(importance_delta_brier),
  sd_recomputed = sd(importance_delta_brier),
  q25_recomputed = quantile(importance_delta_brier, 0.25),
  q75_recomputed = quantile(importance_delta_brier, 0.75)
), by = group_id]
importance_join <- merge(importance_summary, importance_recomputed, by = "group_id")
assert_near(
  importance_join$mean_delta_brier, importance_join$mean_recomputed,
  "mean held-out permutation importance reproduces by predictor group"
)
assert_near(
  importance_join$sd_delta_brier, importance_join$sd_recomputed,
  "importance standard deviations reproduce by predictor group"
)
assert_near(
  importance_join$q25_delta_brier, importance_join$q25_recomputed,
  "importance lower quartiles reproduce by predictor group"
)
assert_near(
  importance_join$q75_delta_brier, importance_join$q75_recomputed,
  "importance upper quartiles reproduce by predictor group"
)

# Figure sources and deidentification ----
assert_true(
  nrow(panel_a) == 100L * 4L * 2L && uniqueN(panel_a$repeat_id) == 100L &&
    setequal(unique(panel_a$model), model_order),
  "Figure 3A source has two metrics for four models across 100 repeats"
)
assert_true(
  nrow(panel_b) == 15L * 100L && uniqueN(panel_b$group_id) == 15L &&
    all(c("FC", "C3_C2", "PRE_grouped", "GIA_grouped") %in% unique(panel_b$group_id)),
  "Figure 3B source has 15 groups and includes FC, C3/C2, PRE, and grouped GIA"
)
panel_b_metabolites <- unique(panel_b[feature_group == "Metabolite", group_id])
assert_true(
  all(selection_summary[metabolite %in% panel_b_metabolites, folds_selected] > 0L),
  "Figure 3B excludes metabolites that were never selected"
)
sample_level_tables <- list(folds, oof, mean_oof, bootstrap_weights, modeling_data)
all_analysis_ids <- unlist(lapply(sample_level_tables, function(x) x$analysis_id), use.names = FALSE)
assert_true(
  all(grepl("^MMA_TPN0_[0-9]{3}$", all_analysis_ids)),
  "all sample-level outputs use deidentified analysis IDs"
)
assert_true(
  !any(grepl("NBS-sample-ID|raw_id", unlist(lapply(sample_level_tables, names)), ignore.case = TRUE)),
  "sample-level output schemas contain no raw-ID fields"
)
assert_true(
  nrow(modeling_data) == 117L && !"TPN_HYPERAL" %in% names(modeling_data),
  "restricted modeling dataset has 117 rows and does not export the TPN source field"
)

report_path <- file.path(qa_dir, "validation_report.txt")
report <- c(
  "RF rerun output validation",
  paste0("Validation time: ", format(Sys.time(), tz = "America/New_York")),
  paste0("Run directory: ", relative_to_project(run_dir, project_root)),
  "",
  checks,
  "",
  paste0("TOTAL: ", length(checks), " checks passed; 0 failed.")
)
writeLines(report, report_path)

cat(paste(report, collapse = "\n"), "\n")
