#!/usr/bin/env Rscript

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
helper_dir <- file.path(project_root, "R")
source(file.path(helper_dir, "ancestry_helpers.R"))
source(file.path(helper_dir, "pre_helpers.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

input_dir <- normalizePath(
  Sys.getenv("HGG_DATA_DIR", file.path(project_root, "data", "raw")),
  mustWork = FALSE
)
results_dir <- normalizePath(
  Sys.getenv("HGG_RESULTS_DIR", file.path(project_root, "results")),
  mustWork = FALSE
)
analysis_dir <- file.path(results_dir, "cohort_characteristics")
table_dir <- file.path(analysis_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Read and match the cohort ----

input_files <- c(
  phenotype_workbook = file.path(
    input_dir, "Scharfelab-NBS1474samples-250207.xlsx"
  ),
  joint_fam = file.path(input_dir, "1000G_378.fam"),
  reference_selection = file.path(input_dir, "sample_pure.txt")
)

phenotype <- read_excel(
  input_files[["phenotype_workbook"]],
  sheet = "1474 DBS with full NBS data",
  .name_repair = "unique_quiet"
)
fam <- fread(input_files[["joint_fam"]], header = FALSE)
reference_selection <- fread(input_files[["reference_selection"]], header = FALSE)

stopifnot(
  !anyDuplicated(phenotype[["NBS-sample-ID"]]),
  !anyDuplicated(reference_selection$V1)
)

is_reference <- fam$V2 %in% reference_selection$V1
stopifnot(
  sum(is_reference) == nrow(reference_selection),
  sum(!is_reference) == 378L
)

phenotype_ids <- as.character(phenotype[["NBS-sample-ID"]])
study_ids <- as.character(fam$V2[!is_reference])

study_match <- match(study_ids, phenotype_ids)
unmatched <- is.na(study_match)
study_match[unmatched] <- match(
  canonical_sample_id(study_ids[unmatched]),
  canonical_sample_id(phenotype_ids)
)
stopifnot(
  !anyNA(study_match),
  !anyDuplicated(study_match),
  length(study_match) == 378L
)
cohort <- phenotype[study_match, , drop = FALSE]

# Derive table variables ----

race_columns <- paste0("RACE_ETH_", 1:4)
stopifnot(all(race_columns %in% names(cohort)))

pre_levels <- c(
  "Hispanic", "White", "East Asian", "Black", "Other/Unknown",
  "Middle Eastern", "South Asian", "Native American"
)

assigned_pre <- apply(
  cohort[, race_columns, drop = FALSE],
  1,
  assign_pre,
  hierarchy = c(
    "Hispanic", "Black", "East Asian", "South Asian", "Middle Eastern",
    "Native American", "White"
  ),
  east_asian_label = "East Asian",
  south_asian_label = "South Asian"
)
reported_pre_count <- apply(
  cohort[, race_columns, drop = FALSE],
  1,
  function(values) {
    values <- trimws(as.character(values))
    sum(!is.na(values) & nzchar(values))
  }
)

outcome <- ifelse(
  cohort[["Group (patient, controls, falsepos)"]] == "patients",
  "True positive",
  "False positive"
)
gender <- fcase(
  cohort[["GENDER"]] == "M", "Male",
  cohort[["GENDER"]] == "F", "Female",
  default = "Missing"
)

tpn_raw <- suppressWarnings(as.numeric(cohort[["TPN_HYPERAL"]]))
tpn <- fcase(
  tpn_raw == 0, "No",
  tpn_raw == 1, "Yes",
  default = "Missing/invalid"
)

birth_weight <- suppressWarnings(as.numeric(cohort[["BIRTH_WT"]]))
collection_age <- suppressWarnings(as.numeric(cohort[["AGE_AT_COLCTN"]]))
screening_year <- suppressWarnings(as.numeric(cohort[["YEAR"]]))
gestational_days <- suppressWarnings(as.numeric(cohort[["GA_DAYS"]]))

stopifnot(
  length(outcome) == 378L,
  sum(outcome == "True positive") == 235L,
  sum(outcome == "False positive") == 143L,
  sum(reported_pre_count >= 2L) == 52L,
  identical(
    as.integer(table(factor(assigned_pre, levels = pre_levels))),
    c(205L, 83L, 34L, 23L, 14L, 9L, 9L, 1L)
  ),
  sum(gender == "Male") == 214L,
  sum(gender == "Female") == 161L,
  sum(gender == "Missing") == 3L,
  sum(tpn == "No") == 303L,
  sum(tpn == "Yes") == 55L,
  sum(tpn == "Missing/invalid") == 20L,
  sum(!is.na(birth_weight)) == 375L,
  sum(!is.na(collection_age)) == 375L
)

format_count <- function(n, denominator = 378L) {
  sprintf("%d (%.1f%%)", n, 100 * n / denominator)
}

format_integer <- function(x) formatC(round(x), format = "f", digits = 0, big.mark = ",")

format_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  q <- unname(quantile(x, probs = c(0.25, 0.5, 0.75), type = 7))
  sprintf(
    "%s (%s-%s)",
    format_integer(q[[2L]]),
    format_integer(q[[1L]]),
    format_integer(q[[3L]])
  )
}

format_median_range <- function(x) {
  x <- x[is.finite(x)]
  sprintf(
    "%.0f (%.0f-%.0f)",
    median(x),
    min(x),
    max(x)
  )
}

section_row <- function(label) {
  data.frame(
    row_type = "section",
    characteristic = label,
    overall = "",
    available_n = NA_integer_,
    stringsAsFactors = FALSE
  )
}

data_row <- function(label, value, available_n = 378L) {
  data.frame(
    row_type = "data",
    characteristic = label,
    overall = value,
    available_n = as.integer(available_n),
    stringsAsFactors = FALSE
  )
}

pre_counts <- table(factor(assigned_pre, levels = pre_levels))

table_rows <- rbind(
  section_row("Outcome"),
  data_row("True positive", format_count(sum(outcome == "True positive"))),
  data_row("False positive", format_count(sum(outcome == "False positive"))),
  section_row("Parent-reported ethnicity (PRE)"),
  do.call(
    rbind,
    lapply(
      seq_along(pre_levels),
      function(index) data_row(pre_levels[[index]], format_count(pre_counts[[index]]))
    )
  ),
  section_row("PRE reporting"),
  data_row(
    "Single/no multiple report",
    format_count(sum(reported_pre_count < 2L))
  ),
  data_row(
    "Multiple reported categories",
    format_count(sum(reported_pre_count >= 2L))
  ),
  section_row("Sex"),
  data_row("Male", format_count(sum(gender == "Male"))),
  data_row("Female", format_count(sum(gender == "Female"))),
  data_row("Missing", format_count(sum(gender == "Missing"))),
  section_row("Clinical characteristics"),
  data_row(
    "Birth weight, g, median (IQR)",
    format_median_iqr(birth_weight),
    sum(!is.na(birth_weight))
  ),
  data_row(
    "Age at blood collection, h, median (IQR)",
    format_median_iqr(collection_age),
    sum(!is.na(collection_age))
  ),
  data_row("Total parenteral nutrition: no", format_count(sum(tpn == "No"))),
  data_row("Total parenteral nutrition: yes", format_count(sum(tpn == "Yes"))),
  data_row(
    "Total parenteral nutrition: missing/invalid",
    format_count(sum(tpn == "Missing/invalid"))
  ),
  data_row(
    "Screening year, median (range)",
    format_median_range(screening_year),
    sum(!is.na(screening_year))
  )
)

write.csv(
  table_rows,
  file.path(table_dir, "cohort_characteristics_table.csv"),
  row.names = FALSE
)

# Characteristics by referral outcome ----

gestational_week <- floor(gestational_days / 7)
gestational_week[
  is.na(gestational_days) | gestational_days < 140 | gestational_days > 315
] <- NA_real_
gestational_category <- fcase(
  is.na(gestational_week), "Unknown/invalid",
  gestational_week > 42, ">42",
  gestational_week == 42, "42",
  gestational_week == 41, "41",
  gestational_week >= 39, "39-40",
  gestational_week >= 37, "37-38",
  gestational_week >= 28, "28-36",
  default = "<28"
)

birth_weight_category <- fcase(
  is.na(birth_weight), "Unknown",
  birth_weight > 5000, ">5000",
  birth_weight >= 4001, "4001-5000",
  birth_weight >= 3501, "3501-4000",
  birth_weight >= 3001, "3001-3500",
  birth_weight >= 2500, "2500-3000",
  birth_weight >= 1000, "1000-2499",
  default = "<1000"
)

collection_age_category <- fcase(
  is.na(collection_age), "Unknown",
  collection_age < 12, "<12",
  collection_age < 24, "12-23",
  collection_age <= 48, "24-48",
  collection_age <= 168, "49-168",
  default = ">168"
)

pre_reporting <- ifelse(
  reported_pre_count >= 2L,
  "Multiple reported categories",
  "Single/no multiple report"
)
gender_stratified <- ifelse(gender == "Missing", "Unknown", gender)
tpn_stratified <- ifelse(tpn == "Missing/invalid", "Unknown", tpn)

fp_denominator <- sum(outcome == "False positive")
tp_denominator <- sum(outcome == "True positive")

stratified_section <- function(label) {
  data.frame(
    row_type = "section",
    label = label,
    false_positive_n = NA_integer_,
    false_positive_denominator = fp_denominator,
    false_positive_display = "",
    true_positive_n = NA_integer_,
    true_positive_denominator = tp_denominator,
    true_positive_display = "",
    stringsAsFactors = FALSE
  )
}

stratified_rows <- function(values, levels) {
  do.call(
    rbind,
    lapply(
      levels,
      function(level) {
        fp_n <- sum(values == level & outcome == "False positive")
        tp_n <- sum(values == level & outcome == "True positive")
        data.frame(
          row_type = "data",
          label = level,
          false_positive_n = fp_n,
          false_positive_denominator = fp_denominator,
          false_positive_display = format_count(fp_n, fp_denominator),
          true_positive_n = tp_n,
          true_positive_denominator = tp_denominator,
          true_positive_display = format_count(tp_n, tp_denominator),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

stratified_table_rows <- rbind(
  stratified_section("Gestational age (wk)"),
  stratified_rows(
    gestational_category,
    c(">42", "42", "41", "39-40", "37-38", "28-36", "<28", "Unknown/invalid")
  ),
  stratified_section("Birth weight (g)"),
  stratified_rows(
    birth_weight_category,
    c(">5000", "4001-5000", "3501-4000", "3001-3500", "2500-3000", "1000-2499", "<1000", "Unknown")
  ),
  stratified_section("Sex"),
  stratified_rows(gender_stratified, c("Male", "Female", "Unknown")),
  stratified_section("Parent-reported ethnicity (PRE)"),
  stratified_rows(assigned_pre, pre_levels),
  stratified_section("PRE reporting"),
  stratified_rows(
    pre_reporting,
    c("Single/no multiple report", "Multiple reported categories")
  ),
  stratified_section("Age at collection (h)"),
  stratified_rows(
    collection_age_category,
    c("<12", "12-23", "24-48", "49-168", ">168", "Unknown")
  ),
  stratified_section("Total parenteral nutrition (TPN)"),
  stratified_rows(tpn_stratified, c("No", "Yes", "Unknown"))
)

for (section_start in which(stratified_table_rows$row_type == "section")) {
  next_section <- which(
    stratified_table_rows$row_type == "section" &
      seq_len(nrow(stratified_table_rows)) > section_start
  )
  section_end <- if (length(next_section)) min(next_section) - 1L else nrow(stratified_table_rows)
  section_rows <- seq.int(section_start + 1L, section_end)
  stopifnot(
    sum(stratified_table_rows$false_positive_n[section_rows]) == fp_denominator,
    sum(stratified_table_rows$true_positive_n[section_rows]) == tp_denominator
  )
}

write.csv(
  stratified_table_rows,
  file.path(table_dir, "cohort_characteristics_stratified_table.csv"),
  row.names = FALSE
)

cat("Wrote aggregate cohort table to", table_dir, "\n")
