#!/usr/bin/env Rscript

# Random-forest analysis of the 117 MMA screen-positive newborns. Metabolites
# are selected within each training fold, before the four predictor sets are fit.

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
helper_dir <- file.path(dirname(script_path), "..", "R")
source(file.path(helper_dir, "project_setup.R"))
source(file.path(helper_dir, "ancestry_helpers.R"))
source(file.path(helper_dir, "pre_helpers.R"))
source(file.path(helper_dir, "statistical_helpers.R"))

required_packages <- c(
  "data.table", "readxl", "randomForest", "pROC", "ggplot2",
  "patchwork"
)
require_packages(required_packages)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(randomForest)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function(args) {
  detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (length(detected_cores) != 1L || is.na(detected_cores)) {
    detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  }
  if (length(detected_cores) != 1L || is.na(detected_cores)) detected_cores <- 1L

  defaults <- list(
    repeats = 100L,
    folds = 10L,
    ntree = 1000L,
    mtry = 7L,
    top_metabolites = 10L,
    bootstrap = 2000L,
    cores = min(4L, as.integer(detected_cores)),
    seed = 20260724L,
    top_importance = 15L,
    gia_reference = "SAS",
    gia_importance = "grouped",
    run_name = "main_117_top10_metabolites"
  )

  string_args <- c("run_name", "gia_reference", "gia_importance")
  integer_args <- setdiff(names(defaults), string_args)
  for (arg in args) {
    if (!grepl("^--[A-Za-z0-9_-]+=", arg)) {
      stop("Arguments must use --name=value syntax. Unrecognized argument: ", arg)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(parts) != 2L) stop("Malformed argument: ", arg)
    key <- gsub("-", "_", parts[[1]])
    if (!key %in% names(defaults)) stop("Unknown argument: --", parts[[1]])

    if (key == "run_name") {
      if (!is_run_name(parts[[2]])) {
        stop("run_name must contain only letters, numbers, underscores, and hyphens.")
      }
      defaults[[key]] <- parts[[2]]
    } else if (key == "gia_reference") {
      value <- toupper(parts[[2]])
      if (!value %in% c("AFR", "AMR", "EAS", "EUR", "SAS")) {
        stop("gia_reference must be one of AFR, AMR, EAS, EUR, or SAS.")
      }
      defaults[[key]] <- value
    } else if (key == "gia_importance") {
      value <- tolower(parts[[2]])
      if (!value %in% c("grouped", "individual")) {
        stop("gia_importance must be grouped or individual.")
      }
      defaults[[key]] <- value
    } else {
      value <- suppressWarnings(as.integer(parts[[2]]))
      if (is.na(value) || value < 1L) stop("Argument must be a positive integer: ", arg)
      defaults[[key]] <- value
    }
  }

  if (!all(integer_args %in% names(defaults))) stop("Internal argument parsing error.")
  defaults
}

paths <- project_paths(script_path)
project_root <- paths$root
analysis_dir <- file.path(paths$results, "mma_model")
input_dir <- paths$data
args <- parse_args(commandArgs(trailingOnly = TRUE))

run_dir <- file.path(analysis_dir, "runs", args$run_name)
table_dir <- file.path(run_dir, "tables")
figure_dir <- file.path(run_dir, "figures")
make_directories(table_dir, figure_dir)

cat("MMA four-model analysis after excluding newborns receiving TPN\n")
cat("Start time:", format(Sys.time(), tz = "America/New_York"), "\n")
cat("Project root:", project_root, "\n")
cat("Run directory:", run_dir, "\n")
cat("Parameters:", paste(names(args), unlist(args), sep = "=", collapse = "; "), "\n")

# Inputs and cohort ----

input_files <- c(
  phenotype_workbook = file.path(input_dir, "Scharfelab-NBS1474samples-250207.xlsx"),
  admixture_q = file.path(input_dir, "1000G_378.5.Q"),
  admixture_fam = file.path(input_dir, "1000G_378.fam"),
  reference_metadata = file.path(input_dir, "all_phase3.psam")
)
require_files(input_files)

rel_path <- function(path) {
  relative_to_project(path, project_root)
}

numeric_clean <- function(x) suppressWarnings(as.numeric(as.character(x)))

global_ancestry <- read_global_ancestry(
  input_files[["admixture_q"]],
  input_files[["admixture_fam"]],
  input_files[["reference_metadata"]]
)
fwrite(global_ancestry$mapping, file.path(table_dir, "global_component_mapping.csv"))
fwrite(global_ancestry$reference_means, file.path(table_dir, "reference_component_means.csv"))

cat("Reading phenotype and metabolite workbook\n")
phenotype <- as.data.table(read_excel(
  input_files[["phenotype_workbook"]],
  sheet = "1474 DBS with full NBS data"
))
id_col <- "NBS-sample-ID"
if (!id_col %in% names(phenotype)) stop("Workbook is missing NBS-sample-ID.")
if (anyDuplicated(phenotype[[id_col]])) stop("Workbook NBS-sample-ID values are not unique.")
if (!all(global_ancestry$study$raw_id %in% phenotype[[id_col]])) {
  stop("Not all 378 normalized study IDs are present in the phenotype workbook.")
}
phenotype <- phenotype[match(global_ancestry$study$raw_id, phenotype[[id_col]])]
if (!identical(as.character(phenotype[[id_col]]), global_ancestry$study$raw_id)) {
  stop("ID-based phenotype ordering failed.")
}

pre_columns <- paste0("RACE_ETH_", 1:4)
if (!all(pre_columns %in% names(phenotype))) {
  stop("Workbook is missing one or more parent-reported ethnicity columns.")
}
for (column in pre_columns) {
  values <- as.character(phenotype[[column]])
  values[values %in% east_asian_pre_labels] <- "EAS"
  values[values %in% south_asian_pre_labels] <- "SAS"
  set(phenotype, j = column, value = values)
}

phenotype[, PRE_category := apply(.SD, 1, assign_pre), .SDcols = pre_columns]

base_metabolite_features <- c(
  "ALA", "ARG", "C02", "C03", "C03DC", "C04", "C05", "C051", "C05DC", "C05OH",
  "C06", "C08", "C081", "C10", "C101", "C12", "C121", "C14", "C141", "C14OH",
  "C16", "C161", "C16OH", "C18", "C181", "C181OH", "C182", "C18OH", "CIT",
  "GLY", "MET", "ORN", "OXP", "PHE", "PRO", "TYR", "VAL", "XLE"
)
metabolite_source_features <- c(base_metabolite_features, "FC")
metabolite_candidate_features <- c(base_metabolite_features, "FC", "C3_C2")
clinical_source <- c("BIRTH_WT", "TPN_HYPERAL", "AGE_AT_COLCTN", "GENDER")
required_columns <- c(
  id_col,
  "Group (patient, controls, falsepos)",
  "Condition or Positive For",
  clinical_source,
  pre_columns,
  metabolite_source_features
)
missing_columns <- setdiff(required_columns, names(phenotype))
if (length(missing_columns) > 0L) {
  stop("Workbook is missing required columns: ", paste(missing_columns, collapse = ", "))
}

numeric_columns <- c(
  "BIRTH_WT", "TPN_HYPERAL", "AGE_AT_COLCTN",
  metabolite_source_features
)
phenotype[, (numeric_columns) := lapply(.SD, numeric_clean), .SDcols = numeric_columns]
phenotype[, male := fifelse(GENDER == "M", 1, fifelse(GENDER == "F", 0, NA_real_))]
phenotype[, outcome := fifelse(
  `Group (patient, controls, falsepos)` == "patients",
  "TP",
  fifelse(
    `Group (patient, controls, falsepos)` %in% c("falsepos284", "falsepos754"),
    "FP",
    NA_character_
  )
)]
if (anyNA(phenotype$outcome)) stop("Unexpected group labels prevent TP/FP assignment.")

phenotype[, C3_C2 := fifelse(C02 > 0, C03 / C02, NA_real_)]
phenotype[, mma_evaluable := !is.na(C03) & (!is.na(C3_C2) | C03 >= 6.3)]
phenotype[, mma_screen_positive := mma_evaluable & (C03 >= 6.3 | C3_C2 >= 0.3)]

joined <- cbind(
  phenotype,
  global_ancestry$study[, .(AFR, AMR, EAS, EUR, SAS)]
)

pre_levels <- c(
  "White", "Black", "Hispanic", "EAS", "SAS", "Middle Eastern",
  "Native American", "Other/Unknown"
)
joined[, PRE_category := factor(PRE_category, levels = pre_levels)]

clinical_matrix_all <- data.frame(
  BIRTH_WT = joined$BIRTH_WT,
  AGE_AT_COLCTN = joined$AGE_AT_COLCTN,
  male = joined$male,
  check.names = FALSE
)
metabolite_matrix_all <- as.data.frame(
  joined[, ..metabolite_candidate_features],
  check.names = FALSE
)
pre_matrix_all <- model.matrix(~ PRE_category, data = joined)[, -1, drop = FALSE]
colnames(pre_matrix_all) <- paste0(
  "PRE_",
  make.names(sub("^PRE_category", "", colnames(pre_matrix_all)), unique = TRUE)
)
all_gia_components <- c("AFR", "AMR", "EAS", "EUR", "SAS")
gia_features <- setdiff(all_gia_components, args$gia_reference)
gia_matrix_all <- as.data.frame(joined[, ..gia_features], check.names = FALSE)
full_candidate_matrix_all <- cbind(
  clinical_matrix_all,
  metabolite_matrix_all,
  pre_matrix_all,
  gia_matrix_all
)

# Apply the cohort criteria in the order shown in the flow table.
stage_masks <- list()
stage_masks[["Matched 378 study samples"]] <- rep(TRUE, nrow(joined))
stage_masks[["MMA rule evaluable"]] <- stage_masks[[1]] & joined$mma_evaluable
stage_masks[["MMA screen-positive"]] <- stage_masks[[2]] & joined$mma_screen_positive
stage_masks[["Birth weight 1,000-5,000 g"]] <- stage_masks[[3]] &
  !is.na(joined$BIRTH_WT) & joined$BIRTH_WT >= 1000 & joined$BIRTH_WT <= 5000
stage_masks[["Age at collection 12-168 h"]] <- stage_masks[[4]] &
  !is.na(joined$AGE_AT_COLCTN) & joined$AGE_AT_COLCTN >= 12 & joined$AGE_AT_COLCTN <= 168
stage_masks[["Exclude newborns receiving TPN"]] <- stage_masks[[5]] &
  !is.na(joined$TPN_HYPERAL) & joined$TPN_HYPERAL == 0
stage_masks[["Complete data for predictors and metabolite candidates"]] <- stage_masks[[6]] &
  complete.cases(full_candidate_matrix_all)

cohort_flow <- rbindlist(lapply(seq_along(stage_masks), function(index) {
  mask <- stage_masks[[index]]
  data.table(
    stage_order = index,
    stage = names(stage_masks)[index],
    n = sum(mask),
    tp = sum(mask & joined$outcome == "TP"),
    fp = sum(mask & joined$outcome == "FP"),
    excluded_from_previous = if (index == 1L) 0L else sum(stage_masks[[index - 1L]]) - sum(mask)
  )
}))
fwrite(cohort_flow, file.path(table_dir, "cohort_flow.csv"))

eligible <- stage_masks[[length(stage_masks)]]
model_data <- joined[eligible]
y <- factor(model_data$outcome, levels = c("FP", "TP"))

clinical_x <- as.matrix(clinical_matrix_all[eligible, , drop = FALSE])
metabolite_x <- as.matrix(metabolite_matrix_all[eligible, , drop = FALSE])
pre_x <- as.matrix(pre_matrix_all[eligible, , drop = FALSE])
gia_x <- as.matrix(gia_matrix_all[eligible, , drop = FALSE])

pre_variable <- apply(pre_x, 2, function(x) length(unique(x)) > 1L)
if (!all(pre_variable)) pre_x <- pre_x[, pre_variable, drop = FALSE]

storage.mode(clinical_x) <- "double"
storage.mode(metabolite_x) <- "double"
storage.mode(pre_x) <- "double"
storage.mode(gia_x) <- "double"

# Predictor sets ----

model_order <- c(
  "Clinical + metabolites",
  "Clinical + metabolites + PRE",
  "Clinical + metabolites + GIA",
  "Clinical + metabolites + PRE + GIA"
)
if (
  nrow(model_data) != 117L ||
    sum(y == "TP") != 85L ||
    sum(y == "FP") != 32L
) {
  stop(
    "Cohort assertion failed after excluding newborns receiving TPN. Expected n=117 (85 TP, 32 FP), observed n=",
    nrow(model_data), " (", sum(y == "TP"), " TP, ", sum(y == "FP"), " FP)."
  )
}
if (!all(model_data$TPN_HYPERAL == 0)) stop("The final cohort contains TPN-exposed newborns.")
if (args$folds > min(table(y))) stop("Number of folds exceeds the smaller outcome class.")
if (args$top_metabolites > ncol(metabolite_x)) {
  stop("top_metabolites exceeds the metabolite-candidate count.")
}
if (ncol(metabolite_x) != 40L ||
    !all(c("FC", "C3_C2") %in% colnames(metabolite_x))) {
  stop("Expected a 40-metabolite candidate pool containing FC and C3/C2.")
}
minimum_predictors <- ncol(clinical_x) + args$top_metabolites
if (args$mtry > minimum_predictors) {
  stop("mtry exceeds the predictor count in the covariates-only model.")
}

sorted_ids <- sort(model_data[[id_col]])
analysis_id_map <- setNames(sprintf("MMA_TPN0_%03d", seq_along(sorted_ids)), sorted_ids)
analysis_ids <- unname(analysis_id_map[model_data[[id_col]]])
if (anyNA(analysis_ids) || anyDuplicated(analysis_ids)) {
  stop("Failed to create deidentified analysis IDs.")
}

modeling_dataset_internal <- as.data.table(cbind(
  data.frame(
    analysis_id = analysis_ids,
    outcome = as.character(y),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  as.data.frame(clinical_x, check.names = FALSE),
  as.data.frame(metabolite_x, check.names = FALSE),
  as.data.frame(pre_x, check.names = FALSE),
  as.data.frame(gia_x, check.names = FALSE)
))
if (anyDuplicated(modeling_dataset_internal$analysis_id)) {
  stop("Deidentified modeling dataset contains duplicate analysis IDs.")
}
fwrite(
  modeling_dataset_internal,
  file.path(table_dir, "modeling_dataset_deidentified_restricted_internal.csv")
)

model_specification <- rbindlist(lapply(model_order, function(model_name) {
  include_pre <- grepl("PRE", model_name, fixed = TRUE)
  include_gia <- grepl("GIA", model_name, fixed = TRUE)
  rbindlist(list(
    data.table(
      model = model_name,
      feature = colnames(clinical_x),
      feature_group = "Clinical",
      role = "always included",
      note = "Clinical predictor; no imputation or tuning"
    ),
    data.table(
      model = model_name,
      feature = colnames(metabolite_x),
      feature_group = "Metabolite candidate",
      role = "fold-wise candidate",
      note = paste0(
        "Eligible for training-fold selection; exactly ",
        args$top_metabolites,
        " candidates ranked by abs(univariate AUC - 0.5) are included"
      )
    ),
    if (include_pre) data.table(
      model = model_name,
      feature = colnames(pre_x),
      feature_group = "PRE",
      role = "always included",
      note = "One-hot indicators; White reference"
    ),
    if (include_gia) data.table(
      model = model_name,
      feature = colnames(gia_x),
      feature_group = "Continuous GIA",
      role = "always included",
      note = paste0(
        "Four of five continuous ADMIXTURE proportions; ",
        args$gia_reference,
        " is the implicit reference"
      )
    )
  ), use.names = TRUE, fill = TRUE)
}))
fwrite(model_specification, file.path(table_dir, "model_specification.csv"))

cohort_summary <- data.table(
  n = nrow(model_data),
  tp = sum(y == "TP"),
  fp = sum(y == "FP"),
  tpn_0 = sum(model_data$TPN_HYPERAL == 0),
  tpn_1 = sum(model_data$TPN_HYPERAL == 1),
  metabolite_candidates = ncol(metabolite_x),
  metabolites_selected_per_fold = args$top_metabolites,
  fc_candidate = "FC" %in% colnames(metabolite_x),
  c3_c2_candidate = "C3_C2" %in% colnames(metabolite_x),
  clinical_metabolite_predictors_per_fold = ncol(clinical_x) + args$top_metabolites,
  clinical_metabolite_pre_predictors_per_fold =
    ncol(clinical_x) + args$top_metabolites + ncol(pre_x),
  clinical_metabolite_gia_predictors_per_fold =
    ncol(clinical_x) + args$top_metabolites + ncol(gia_x),
  full_predictors_per_fold =
    ncol(clinical_x) + args$top_metabolites + ncol(pre_x) + ncol(gia_x),
  pre_reference = "White",
  omitted_gia_component = args$gia_reference
)
fwrite(cohort_summary, file.path(table_dir, "cohort_summary.csv"))

# Cross-validation folds ----

fold_list <- lapply(seq_len(args$repeats), function(repeat_id) {
  make_stratified_folds(y, args$folds, args$seed + repeat_id * 10000L)
})
if (any(vapply(fold_list, function(fold) {
  any(vapply(seq_len(args$folds), function(fold_id) {
    length(unique(y[fold == fold_id])) != 2L
  }, logical(1)))
}, logical(1)))) {
  stop("At least one held-out fold lacks one of the outcome classes.")
}
fold_assignments <- rbindlist(lapply(seq_len(args$repeats), function(repeat_id) {
  data.table(
    repeat_id = repeat_id,
    analysis_id = analysis_ids,
    outcome = as.character(y),
    fold = fold_list[[repeat_id]]
  )
}))
fwrite(fold_assignments, file.path(table_dir, "fold_assignments.csv"))

full_model_name <- model_order[[4]]
pre_feature_columns <- colnames(pre_x)
gia_feature_columns <- colnames(gia_x)

clinical_labels <- c(
  BIRTH_WT = "Birth weight",
  AGE_AT_COLCTN = "Age at collection",
  male = "Sex"
)
individual_clinical_groups <- lapply(colnames(clinical_x), function(column) column)
names(individual_clinical_groups) <- colnames(clinical_x)
individual_metabolite_groups <- lapply(
  colnames(metabolite_x),
  function(column) column
)
names(individual_metabolite_groups) <- colnames(metabolite_x)
gia_importance_groups <- if (args$gia_importance == "grouped") {
  list(GIA_grouped = gia_feature_columns)
} else {
  setNames(
    lapply(gia_feature_columns, function(column) column),
    paste0("GIA_", gia_feature_columns)
  )
}
importance_groups <- c(
  individual_clinical_groups,
  individual_metabolite_groups,
  list(PRE_grouped = pre_feature_columns),
  gia_importance_groups
)
importance_group_metadata <- rbindlist(lapply(names(importance_groups), function(group_id) {
  columns <- importance_groups[[group_id]]
  if (group_id == "PRE_grouped") {
    display_label <- "PRE (grouped)"
    feature_group <- "PRE"
  } else if (group_id == "GIA_grouped") {
    display_label <- "GIA (grouped)"
    feature_group <- "GIA"
  } else if (grepl("^GIA_", group_id)) {
    display_label <- sub("^GIA_", "", group_id)
    feature_group <- "GIA"
  } else if (group_id %in% names(clinical_labels)) {
    display_label <- unname(clinical_labels[group_id])
    feature_group <- "Clinical"
  } else {
    display_label <- if (group_id == "C3_C2") "C3/C2" else group_id
    feature_group <- "Metabolite"
  }
  data.table(
    group_id = group_id,
    display_label = display_label,
    feature_group = feature_group,
    columns = paste(columns, collapse = ";"),
    n_columns = length(columns)
  )
}))
fwrite(importance_group_metadata, file.path(table_dir, "importance_group_definition.csv"))

# Fold-wise selection and fitting ----

rank_metabolites_in_training <- function(train_index) {
  ranking <- rbindlist(lapply(colnames(metabolite_x), function(feature) {
    values <- metabolite_x[train_index, feature]
    auc_value <- if (length(unique(values)) < 2L) {
      NA_real_
    } else {
      tryCatch(
        as.numeric(pROC::auc(pROC::roc(
          response = y[train_index],
          predictor = values,
          levels = c("FP", "TP"),
          direction = "<",
          quiet = TRUE
        ))),
        error = function(error) NA_real_
      )
    }
    data.table(
      metabolite = feature,
      univariate_auc = auc_value,
      absolute_auc_distance = ifelse(
        is.finite(auc_value),
        abs(auc_value - 0.5),
        -Inf
      )
    )
  }))
  setorder(ranking, -absolute_auc_distance, metabolite)
  ranking[, rank := seq_len(.N)]
  ranking[, selected := rank <= args$top_metabolites]
  if (
    sum(ranking$selected) != args$top_metabolites ||
      any(!is.finite(ranking[selected == TRUE]$absolute_auc_distance))
  ) {
    stop("Could not select exactly ", args$top_metabolites, " finite metabolite candidates.")
  }
  ranking
}

fit_one_repeat <- function(repeat_id) {
  fold <- fold_list[[repeat_id]]
  predictions <- setNames(
    lapply(model_order, function(model_name) rep(NA_real_, length(y))),
    model_order
  )
  permuted_predictions <- setNames(
    lapply(names(importance_groups), function(group_id) rep(NA_real_, length(y))),
    names(importance_groups)
  )
  selection_rows <- vector("list", args$folds)
  fit_audit_rows <- vector("list", args$folds * length(model_order))
  fit_audit_index <- 0L

  for (fold_id in seq_len(args$folds)) {
    test_index <- which(fold == fold_id)
    train_index <- which(fold != fold_id)
    if (length(intersect(train_index, test_index)) != 0L) {
      stop("Training/test overlap in repeat ", repeat_id, ", fold ", fold_id, ".")
    }
    ranking <- rank_metabolites_in_training(train_index)
    selected_metabolites <- ranking[selected == TRUE]$metabolite
    ranking[, `:=`(
      repeat_id = repeat_id,
      fold = fold_id,
      n_train = length(train_index),
      train_tp = sum(y[train_index] == "TP"),
      train_fp = sum(y[train_index] == "FP")
    )]
    setcolorder(
      ranking,
      c(
        "repeat_id", "fold", "n_train", "train_tp", "train_fp",
        "metabolite", "univariate_auc", "absolute_auc_distance",
        "rank", "selected"
      )
    )
    selection_rows[[fold_id]] <- ranking

    core_x_fold <- cbind(
      clinical_x,
      metabolite_x[, selected_metabolites, drop = FALSE]
    )
    model_matrices_fold <- list(
      "Clinical + metabolites" = core_x_fold,
      "Clinical + metabolites + PRE" = cbind(core_x_fold, pre_x),
      "Clinical + metabolites + GIA" = cbind(core_x_fold, gia_x),
      "Clinical + metabolites + PRE + GIA" = cbind(core_x_fold, pre_x, gia_x)
    )
    model_matrices_fold <- model_matrices_fold[model_order]
    if (any(vapply(model_matrices_fold, nrow, integer(1)) != length(y))) {
      stop("Fold-specific model matrix row-count mismatch.")
    }

    for (model_index in seq_along(model_order)) {
      model_name <- model_order[[model_index]]
      x <- model_matrices_fold[[model_name]]
      fit_seed <- args$seed + repeat_id * 100000L + fold_id * 1000L + model_index * 10L
      set.seed(fit_seed)
      random_forest <- randomForest(
        x = x[train_index, , drop = FALSE],
        y = y[train_index],
        ntree = args$ntree,
        mtry = args$mtry,
        nodesize = 1L,
        importance = FALSE,
        keep.forest = TRUE
      )
      if (
        random_forest$ntree != args$ntree ||
          random_forest$mtry != args$mtry
      ) {
        stop("Fitted forest did not retain the requested ntree/mtry values.")
      }

      predictions[[model_name]][test_index] <- predict(
        random_forest,
        x[test_index, , drop = FALSE],
        type = "prob"
      )[, "TP"]

      fit_audit_index <- fit_audit_index + 1L
      fit_audit_rows[[fit_audit_index]] <- data.table(
        repeat_id = repeat_id,
        fold = fold_id,
        model = model_name,
        n_train = length(train_index),
        n_test = length(test_index),
        predictors = ncol(x),
        selected_metabolites = paste(selected_metabolites, collapse = ";"),
        fit_seed = fit_seed,
        requested_ntree = args$ntree,
        fitted_ntree = random_forest$ntree,
        requested_mtry = args$mtry,
        fitted_mtry = random_forest$mtry
      )

      if (identical(model_name, full_model_name)) {
        full_test <- x[test_index, , drop = FALSE]
        baseline_test_prediction <- predictions[[model_name]][test_index]
        for (group_index in seq_along(importance_groups)) {
          group_id <- names(importance_groups)[[group_index]]
          group_columns <- intersect(
            importance_groups[[group_id]],
            colnames(full_test)
          )
          # A metabolite not selected in this training fold contributes zero
          # to pipeline-level importance for the fold. Selection frequency is
          # reported separately.
          permuted_predictions[[group_id]][test_index] <- baseline_test_prediction
          if (length(group_columns) == 0L) next
          set.seed(
            args$seed + repeat_id * 100000L + fold_id * 1000L +
              500L + group_index
          )
          permutation <- sample(seq_along(test_index))
          permuted_test <- full_test
          permuted_test[, group_columns] <- full_test[permutation, group_columns, drop = FALSE]
          permuted_predictions[[group_id]][test_index] <- predict(
            random_forest,
            permuted_test,
            type = "prob"
          )[, "TP"]
        }
      }
    }
  }

  if (any(vapply(predictions, anyNA, logical(1)))) {
    stop("Missing out-of-fold model prediction in repeat ", repeat_id, ".")
  }
  if (any(vapply(permuted_predictions, anyNA, logical(1)))) {
    stop("Missing permuted out-of-fold prediction in repeat ", repeat_id, ".")
  }

  metrics <- rbindlist(lapply(model_order, function(model_name) {
    performance <- classification_metrics(y, predictions[[model_name]])
    data.table(
      repeat_id = repeat_id,
      model = model_name,
      auc = performance$auc,
      specificity_at_95_sensitivity = performance$specificity95,
      brier_score = performance$brier,
      operating_threshold = performance$operating_point$threshold,
      achieved_sensitivity = performance$operating_point$sensitivity
    )
  }))

  prediction_table <- rbindlist(lapply(model_order, function(model_name) {
    data.table(
      repeat_id = repeat_id,
      analysis_id = analysis_ids,
      outcome = as.character(y),
      model = model_name,
      probability = predictions[[model_name]]
    )
  }))

  outcome_numeric <- as.numeric(y == "TP")
  baseline_brier <- mean((predictions[[full_model_name]] - outcome_numeric)^2)
  importance <- rbindlist(lapply(names(importance_groups), function(group_id) {
    permuted_brier <- mean((permuted_predictions[[group_id]] - outcome_numeric)^2)
    data.table(
      repeat_id = repeat_id,
      model = full_model_name,
      group_id = group_id,
      baseline_brier = baseline_brier,
      permuted_brier = permuted_brier,
      importance_delta_brier = permuted_brier - baseline_brier
    )
  }))

  list(
    predictions = prediction_table,
    metrics = metrics,
    importance = importance,
    metabolite_selection = rbindlist(selection_rows),
    fit_audit = rbindlist(fit_audit_rows)
  )
}

cat(
  "Running ", args$repeats, " repeated ", args$folds,
  "-fold cross-validation runs with ", args$ntree,
  " trees per model; cores=", args$cores, "\n",
  sep = ""
)

if (.Platform$OS.type == "unix" && args$cores > 1L) {
  cv_results <- parallel::mclapply(
    seq_len(args$repeats),
    fit_one_repeat,
    mc.cores = args$cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  cv_results <- lapply(seq_len(args$repeats), function(repeat_id) {
    if (repeat_id %% 10L == 0L) cat("Finished repeat", repeat_id, "\n")
    fit_one_repeat(repeat_id)
  })
}

if (any(vapply(cv_results, inherits, logical(1), what = "try-error"))) {
  stop("One or more cross-validation workers failed.")
}

# Summaries across repeats ----

oof_predictions <- rbindlist(lapply(cv_results, `[[`, "predictions"))
repeat_metrics <- rbindlist(lapply(cv_results, `[[`, "metrics"))
importance_by_repeat <- rbindlist(lapply(cv_results, `[[`, "importance"))
metabolite_selection <- rbindlist(
  lapply(cv_results, `[[`, "metabolite_selection")
)
fit_audit <- rbindlist(lapply(cv_results, `[[`, "fit_audit"))
importance_by_repeat <- merge(
  importance_by_repeat,
  importance_group_metadata,
  by = "group_id",
  all.x = TRUE,
  sort = FALSE
)

fwrite(oof_predictions, file.path(table_dir, "oof_predictions_by_repeat.csv"))
fwrite(repeat_metrics, file.path(table_dir, "repeat_metrics.csv"))
fwrite(importance_by_repeat, file.path(table_dir, "permutation_importance_by_repeat.csv"))
fwrite(metabolite_selection, file.path(table_dir, "metabolite_selection_by_fold.csv"))
fwrite(fit_audit, file.path(table_dir, "forest_fit_audit.csv"))

expected_oof_rows <- nrow(model_data) * args$repeats * length(model_order)
expected_selection_rows <-
  args$repeats * args$folds * ncol(metabolite_x)
expected_fit_rows <- args$repeats * args$folds * length(model_order)
if (nrow(oof_predictions) != expected_oof_rows) {
  stop("Unexpected OOF-prediction row count.")
}
if (nrow(metabolite_selection) != expected_selection_rows) {
  stop("Unexpected metabolite-selection row count.")
}
if (nrow(fit_audit) != expected_fit_rows) {
  stop("Unexpected forest-fit audit row count.")
}
if (nrow(importance_by_repeat) != args$repeats * length(importance_groups)) {
  stop("Unexpected held-out permutation-importance row count.")
}
oof_count_check <- oof_predictions[, .N, by = .(analysis_id, model)]
if (any(oof_count_check$N != args$repeats)) {
  stop("A subject/model does not have exactly one OOF prediction per repeat.")
}
selection_count_check <- metabolite_selection[, .(
  candidates = .N,
  selected = sum(selected),
  unique_selected = uniqueN(metabolite[selected == TRUE])
), by = .(repeat_id, fold)]
if (any(selection_count_check$candidates != ncol(metabolite_x)) ||
    any(selection_count_check$selected != args$top_metabolites) ||
    any(selection_count_check$unique_selected != args$top_metabolites)) {
  stop("A repeat/fold did not select exactly ten unique metabolites.")
}
if (any(repeat_metrics$achieved_sensitivity < 0.95)) {
  stop("A repeat-level operating point failed to attain sensitivity >=0.95.")
}
if (any(!is.finite(oof_predictions$probability))) {
  stop("Non-finite out-of-fold probability.")
}
if (any(oof_predictions$probability < 0 | oof_predictions$probability > 1)) {
  stop("Out-of-fold probability outside [0, 1].")
}
if (any(fit_audit$fitted_ntree != args$ntree) ||
    any(fit_audit$fitted_mtry != args$mtry)) {
  stop("Forest-fit audit found an ntree/mtry mismatch.")
}
expected_predictor_counts <- c(
  "Clinical + metabolites" = ncol(clinical_x) + args$top_metabolites,
  "Clinical + metabolites + PRE" =
    ncol(clinical_x) + args$top_metabolites + ncol(pre_x),
  "Clinical + metabolites + GIA" =
    ncol(clinical_x) + args$top_metabolites + ncol(gia_x),
  "Clinical + metabolites + PRE + GIA" =
    ncol(clinical_x) + args$top_metabolites + ncol(pre_x) + ncol(gia_x)
)
if (any(fit_audit$predictors != expected_predictor_counts[fit_audit$model])) {
  stop("A forest used an unexpected number of predictors.")
}
panel_consistency <- fit_audit[, .(
  selected_panels = uniqueN(selected_metabolites),
  model_count = uniqueN(model)
), by = .(repeat_id, fold)]
if (any(panel_consistency$selected_panels != 1L) ||
    any(panel_consistency$model_count != length(model_order))) {
  stop("The selected metabolite panel was not shared by all four models.")
}

metabolite_selection_summary <- metabolite_selection[, .(
  folds_evaluated = .N,
  folds_selected = sum(selected),
  selection_fraction = mean(selected),
  mean_univariate_auc = mean(univariate_auc),
  mean_absolute_auc_distance = mean(absolute_auc_distance),
  median_rank = median(rank),
  minimum_rank = min(rank),
  maximum_rank = max(rank)
), by = metabolite]
setorder(
  metabolite_selection_summary,
  -selection_fraction,
  -mean_absolute_auc_distance,
  metabolite
)
fwrite(
  metabolite_selection_summary,
  file.path(table_dir, "metabolite_selection_summary.csv")
)

repeat_stability <- melt(
  repeat_metrics,
  id.vars = c("repeat_id", "model"),
  variable.name = "metric",
  value.name = "value"
)[, .(
  repeats = .N,
  mean = mean(value),
  sd = sd(value),
  median = median(value),
  q25 = quantile(value, 0.25),
  q75 = quantile(value, 0.75),
  min = min(value),
  max = max(value),
  interpretation = "Algorithmic stability across dependent CV repeats; not an inferential confidence interval"
), by = .(model, metric)]
fwrite(repeat_stability, file.path(table_dir, "repeat_stability_summary.csv"))

importance_summary <- importance_by_repeat[, .(
  repeats = .N,
  mean_delta_brier = mean(importance_delta_brier),
  sd_delta_brier = sd(importance_delta_brier),
  median_delta_brier = median(importance_delta_brier),
  q25_delta_brier = quantile(importance_delta_brier, 0.25),
  q75_delta_brier = quantile(importance_delta_brier, 0.75),
  proportion_positive = mean(importance_delta_brier > 0)
), by = .(group_id, display_label, feature_group, columns, n_columns)]
importance_summary <- merge(
  importance_summary,
  metabolite_selection_summary[, .(
    group_id = metabolite,
    folds_evaluated,
    folds_selected,
    selection_fraction
  )],
  by = "group_id",
  all.x = TRUE,
  sort = FALSE
)
importance_summary[is.na(selection_fraction), `:=`(
  folds_evaluated = args$repeats * args$folds,
  folds_selected = args$repeats * args$folds,
  selection_fraction = 1
)]
setorder(importance_summary, -mean_delta_brier)
fwrite(importance_summary, file.path(table_dir, "permutation_importance_summary.csv"))

mean_oof_predictions <- oof_predictions[, .(
  outcome = outcome[1],
  probability = mean(probability)
), by = .(analysis_id, model)]
fwrite(mean_oof_predictions, file.path(table_dir, "mean_oof_predictions.csv"))

point_results <- lapply(model_order, function(model_name) {
  model_predictions <- mean_oof_predictions[model == model_name]
  performance <- classification_metrics(model_predictions$outcome, model_predictions$probability)
  list(
    performance = data.table(
      model = model_name,
      metric = c("AUC", "Specificity at >=95% sensitivity", "Brier score"),
      estimate = c(performance$auc, performance$specificity95, performance$brier)
    ),
    operating_point = data.table(
      model = model_name,
      threshold = performance$operating_point$threshold,
      achieved_sensitivity = performance$operating_point$sensitivity,
      specificity = performance$operating_point$specificity,
      tp = performance$operating_point$tp,
      fn = performance$operating_point$fn,
      tn = performance$operating_point$tn,
      fp = performance$operating_point$fp
    )
  )
})
point_performance <- rbindlist(lapply(point_results, `[[`, "performance"))
point_operating_points <- rbindlist(lapply(point_results, `[[`, "operating_point"))
if (any(point_operating_points$achieved_sensitivity < 0.95)) {
  stop("A primary operating point failed to attain sensitivity >=0.95.")
}
fwrite(point_operating_points, file.path(table_dir, "primary_operating_points.csv"))

# Subject-level bootstrap inference ----

set.seed(args$seed + 900000L)
subject_table <- unique(mean_oof_predictions[, .(analysis_id, outcome)])
true_positive_index <- which(subject_table$outcome == "TP")
false_positive_index <- which(subject_table$outcome == "FP")
bootstrap_results <- lapply(seq_len(args$bootstrap), function(bootstrap_id) {
  sampled_rows <- c(
    sample(true_positive_index, length(true_positive_index), replace = TRUE),
    sample(false_positive_index, length(false_positive_index), replace = TRUE)
  )
  sampled_ids <- subject_table$analysis_id[sampled_rows]
  sampled_outcome <- subject_table$outcome[sampled_rows]

  metrics <- rbindlist(lapply(model_order, function(model_name) {
    probability_lookup <- setNames(
      mean_oof_predictions[model == model_name]$probability,
      mean_oof_predictions[model == model_name]$analysis_id
    )
    sampled_probability <- unname(probability_lookup[sampled_ids])
    performance <- classification_metrics(sampled_outcome, sampled_probability)
    data.table(
      bootstrap_id = bootstrap_id,
      model = model_name,
      auc = performance$auc,
      specificity_at_95_sensitivity = performance$specificity95,
      brier_score = performance$brier,
      operating_threshold = performance$operating_point$threshold,
      achieved_sensitivity = performance$operating_point$sensitivity
    )
  }))
  multiplicity <- tabulate(
    match(sampled_ids, subject_table$analysis_id),
    nbins = nrow(subject_table)
  )
  weights <- data.table(
    bootstrap_id = bootstrap_id,
    analysis_id = subject_table$analysis_id,
    outcome = subject_table$outcome,
    multiplicity = multiplicity
  )
  list(metrics = metrics, weights = weights)
})
bootstrap_metrics <- rbindlist(lapply(bootstrap_results, `[[`, "metrics"))
bootstrap_weights <- rbindlist(lapply(bootstrap_results, `[[`, "weights"))
fwrite(bootstrap_metrics, file.path(table_dir, "subject_bootstrap_metrics.csv"))
fwrite(bootstrap_weights, file.path(table_dir, "subject_bootstrap_weights.csv"))
if (any(bootstrap_metrics$achieved_sensitivity < 0.95)) {
  stop("A bootstrap operating point failed to attain sensitivity >=0.95.")
}
bootstrap_weight_check <- bootstrap_weights[, .(
  sampled_tp = sum(multiplicity[outcome == "TP"]),
  sampled_fp = sum(multiplicity[outcome == "FP"])
), by = bootstrap_id]
if (any(bootstrap_weight_check$sampled_tp != sum(y == "TP")) ||
    any(bootstrap_weight_check$sampled_fp != sum(y == "FP"))) {
  stop("Outcome-stratified bootstrap did not preserve class sample sizes.")
}

bootstrap_long <- melt(
  bootstrap_metrics,
  id.vars = c("bootstrap_id", "model"),
  measure.vars = c(
    "auc", "specificity_at_95_sensitivity", "brier_score"
  ),
  variable.name = "metric_code",
  value.name = "bootstrap_estimate"
)
metric_lookup <- c(
  auc = "AUC",
  specificity_at_95_sensitivity = "Specificity at >=95% sensitivity",
  brier_score = "Brier score"
)
bootstrap_long[, metric := unname(metric_lookup[metric_code])]
bootstrap_intervals <- bootstrap_long[, .(
  ci_low = quantile_ci(bootstrap_estimate)[1],
  ci_high = quantile_ci(bootstrap_estimate)[2]
), by = .(model, metric)]
primary_performance <- merge(
  point_performance,
  bootstrap_intervals,
  by = c("model", "metric"),
  all.x = TRUE,
  sort = FALSE
)
primary_performance[, uncertainty := paste0(
  args$bootstrap,
  "-replicate outcome-stratified subject bootstrap on mean OOF predictions"
)]
primary_performance[, model := factor(model, levels = model_order)]
setorder(primary_performance, model, metric)
primary_performance[, model := as.character(model)]
fwrite(primary_performance, file.path(table_dir, "primary_performance.csv"))

contrast_definitions <- data.table(
  contrast = c(
    "GIA increment without PRE",
    "GIA increment conditional on PRE",
    "PRE increment without GIA",
    "PRE increment conditional on GIA"
  ),
  augmented_model = c(model_order[3], model_order[4], model_order[2], model_order[4]),
  reference_model = c(model_order[1], model_order[2], model_order[1], model_order[3])
)

point_lookup <- dcast(point_performance, metric ~ model, value.var = "estimate")
bootstrap_wide <- dcast(
  bootstrap_long,
  bootstrap_id + metric ~ model,
  value.var = "bootstrap_estimate"
)
paired_effects <- rbindlist(lapply(seq_len(nrow(contrast_definitions)), function(index) {
  definition <- contrast_definitions[index]
  rbindlist(lapply(unique(point_performance$metric), function(metric_name) {
    point_row <- point_lookup[metric == metric_name]
    bootstrap_row <- bootstrap_wide[metric == metric_name]
    point_delta <- point_row[[definition$augmented_model]] - point_row[[definition$reference_model]]
    bootstrap_delta <- bootstrap_row[[definition$augmented_model]] -
      bootstrap_row[[definition$reference_model]]
    data.table(
      contrast = definition$contrast,
      metric = metric_name,
      augmented_model = definition$augmented_model,
      reference_model = definition$reference_model,
      estimate = point_delta,
      ci_low = quantile_ci(bootstrap_delta)[1],
      ci_high = quantile_ci(bootstrap_delta)[2],
      favorable_direction = ifelse(metric_name == "Brier score", "Negative", "Positive"),
      inference_note = "Paired outcome-stratified subject bootstrap; repeated CV runs were averaged before inference"
    )
  }))
}))
fwrite(paired_effects, file.path(table_dir, "paired_effects.csv"))

# Figure source tables ----

performance_plot_source <- melt(
  repeat_metrics,
  id.vars = c("repeat_id", "model"),
  measure.vars = c("auc", "specificity_at_95_sensitivity"),
  variable.name = "metric_code",
  value.name = "estimate"
)
performance_plot_source[, metric := fifelse(
  metric_code == "auc",
  "AUC",
  "Specificity at sensitivity >=0.95"
)]
performance_plot_source[, model := factor(model, levels = model_order)]
fwrite(performance_plot_source, file.path(table_dir, "figure_panel_a_source.csv"))

forced_groups <- c(
  "PRE_grouped",
  names(importance_groups)[grepl("^GIA_", names(importance_groups))]
)
other_top <- importance_summary[
  !group_id %in% forced_groups &
    (feature_group != "Metabolite" | folds_selected > 0)
][
  seq_len(min(max(args$top_importance - length(forced_groups), 0L), .N))
]$group_id
selected_importance_groups <- unique(c(forced_groups, other_top))
importance_plot_source <- importance_by_repeat[group_id %in% selected_importance_groups]
importance_order <- importance_summary[group_id %in% selected_importance_groups][
  order(mean_delta_brier)
]$display_label
importance_plot_source[, display_label := factor(display_label, levels = importance_order)]
fwrite(importance_plot_source, file.path(table_dir, "figure_panel_b_source.csv"))

model_colors <- c(
  "Clinical + metabolites" = "#7F7F7F",
  "Clinical + metabolites + PRE" = "#0072B2",
  "Clinical + metabolites + GIA" = "#D55E00",
  "Clinical + metabolites + PRE + GIA" = "#009E73"
)
model_plot_labels <- c(
  "Clinical + metabolites" = "Clinical + metabolites",
  "Clinical + metabolites + PRE" = "+ PRE",
  "Clinical + metabolites + GIA" = "+ GIA",
  "Clinical + metabolites + PRE + GIA" = "+ PRE + GIA"
)
importance_colors <- c(
  Clinical = "#6C757D",
  Metabolite = "#56B4E9",
  PRE = "#0072B2",
  GIA = "#D55E00"
)

# Diagnostic plot ----

theme_publication <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      plot.margin = margin(6, 8, 6, 6)
    )
}

panel_a <- ggplot(
  performance_plot_source,
  aes(x = model, y = estimate, fill = model, color = model)
) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.20, linewidth = 0.45) +
  geom_boxplot(
    width = 0.20,
    outlier.shape = NA,
    alpha = 0.72,
    linewidth = 0.45
  ) +
  geom_jitter(width = 0.075, height = 0, alpha = 0.18, size = 0.65, stroke = 0) +
  facet_wrap(~ metric, ncol = 2) +
  coord_flip() +
  scale_fill_manual(values = model_colors, labels = model_plot_labels) +
  scale_color_manual(values = model_colors, labels = model_plot_labels) +
  scale_x_discrete(labels = model_plot_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "A. Performance across repeated cross-validation runs",
    x = NULL,
    y = "Performance estimate"
  ) +
  theme_publication(9.5) +
  theme(legend.position = "none")

importance_means <- importance_plot_source[, .(
  mean_importance = mean(importance_delta_brier)
), by = .(display_label, feature_group)]

panel_b <- ggplot(
  importance_plot_source,
  aes(x = importance_delta_brier, y = display_label, fill = feature_group)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.45) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.62, linewidth = 0.42) +
  geom_point(
    data = importance_means,
    aes(x = mean_importance, y = display_label),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "black",
    size = 1.8,
    stroke = 0.45
  ) +
  scale_fill_manual(values = importance_colors) +
  labs(
    title = "B. Held-out permutation importance in the full model",
    x = "Increase in Brier score after permutation",
    y = NULL,
    fill = NULL
  ) +
  theme_publication(9.5) +
  theme(legend.position = "bottom")

figure_caption_text <- paste0(
  "Panel A distributions describe algorithmic stability; repeated CV runs are dependent and are not inferential confidence intervals. ",
  "Within each training fold, 10 of 40 metabolite candidates were selected by absolute univariate AUC distance from 0.5 and shared across all four models. ",
  "Panel B permutes each included clinical/metabolite predictor separately and PRE or GIA jointly within held-out folds; an unselected metabolite has zero pipeline-level contribution in that fold."
)

figure <- (panel_a | panel_b) +
  plot_layout(widths = c(1.22, 1)) +
  plot_annotation(
    subtitle = paste0(
      "MMA screen-positive newborns after excluding those receiving TPN: n=", nrow(model_data),
      " (", sum(y == "TP"), " TP, ", sum(y == "FP"), " FP); ",
      args$repeats, " repeated stratified ", args$folds, "-fold CV runs"
    ),
    caption = paste(strwrap(figure_caption_text, width = 180), collapse = "\n"),
    theme = theme(
      plot.subtitle = element_text(size = 10.5, hjust = 0),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )

pdf_path <- file.path(figure_dir, "mma_tpn0_four_model_performance_importance.pdf")
png_path <- file.path(figure_dir, "mma_tpn0_four_model_performance_importance.png")
ggsave(pdf_path, figure, width = 13.5, height = 7.4, device = cairo_pdf)
ggsave(png_path, figure, width = 13.5, height = 7.4, dpi = 300, bg = "white")

cat("\nCohort summary:\n")
print(cohort_summary)
cat("\nPrimary performance:\n")
print(primary_performance)
cat("\nGIA-related paired effects:\n")
print(paired_effects[grepl("GIA increment", contrast)])
cat("\nMetabolite selection frequency:\n")
print(metabolite_selection_summary[seq_len(min(15L, .N))])
cat("\nTop permutation-importance groups:\n")
print(importance_summary[seq_len(min(15L, .N))])
cat("\nFigures:\n", rel_path(pdf_path), "\n", rel_path(png_path), "\n", sep = "")
cat("End time:", format(Sys.time(), tz = "America/New_York"), "\n")
