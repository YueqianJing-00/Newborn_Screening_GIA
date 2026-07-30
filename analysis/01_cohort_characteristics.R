#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(readxl)
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Could not determine script path from --file.")
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

script_path <- get_script_path()
source(file.path(dirname(script_path), "..", "R", "project_paths.R"))
paths <- get_release_paths(script_path)
project_root <- paths$root
input_dir <- paths$data
analysis_dir <- file.path(paths$results, "cohort_characteristics")
table_dir <- file.path(analysis_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  phenotype_workbook = file.path(
    input_dir, "Scharfelab-NBS1474samples-250207.xlsx"
  ),
  joint_fam = file.path(input_dir, "1000G_378.fam"),
  reference_selection = file.path(input_dir, "sample_pure.txt")
)
stopifnot(all(file.exists(input_files)))

input_checksums <- data.frame(
  input = names(input_files),
  file = basename(input_files),
  bytes = as.numeric(file.info(input_files)$size),
  sha256 = vapply(
    input_files,
    function(path) digest(file = path, algo = "sha256", serialize = FALSE),
    character(1)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  input_checksums,
  file.path(table_dir, "input_checksums_sha256.csv"),
  row.names = FALSE
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
stopifnot(sum(is_reference) == 2158L, sum(!is_reference) == 378L)

canonical_sample_key <- function(x) sub("^NBSfalsepos_", "", as.character(x))
phenotype_ids <- as.character(phenotype[["NBS-sample-ID"]])
study_ids <- as.character(fam$V2[!is_reference])

study_match <- match(study_ids, phenotype_ids)
unmatched <- is.na(study_match)
study_match[unmatched] <- match(
  canonical_sample_key(study_ids[unmatched]),
  canonical_sample_key(phenotype_ids)
)
stopifnot(
  !anyNA(study_match),
  !anyDuplicated(study_match),
  length(study_match) == 378L
)
cohort <- phenotype[study_match, , drop = FALSE]

race_columns <- paste0("RACE_ETH_", 1:4)
stopifnot(all(race_columns %in% names(cohort)))

eas_labels <- c("Japanese", "Chinese", "Laos", "Korean", "Vietnamese", "Filipino")
sas_labels <- "Asian East Indian"

clean_pre_values <- function(values) {
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  values[values %in% eas_labels] <- "East Asian"
  values[values %in% sas_labels] <- "South Asian"
  unique(values)
}

assign_pre_v8 <- function(values) {
  values <- clean_pre_values(values)
  if ("Hispanic" %in% values) return("Hispanic")
  if ("Black" %in% values) return("Black")
  if ("East Asian" %in% values) return("East Asian")
  if ("South Asian" %in% values) return("South Asian")
  if ("Middle Eastern" %in% values) return("Middle Eastern")
  if ("Native American" %in% values) return("Native American")
  if (identical(values, "White")) return("White")
  "Other/Unknown"
}

pre_levels <- c(
  "Hispanic", "White", "East Asian", "Black", "Other/Unknown",
  "Middle Eastern", "South Asian", "Native American"
)

assigned_pre <- apply(cohort[, race_columns, drop = FALSE], 1, assign_pre_v8)
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
gender <- ifelse(
  cohort[["GENDER"]] == "M",
  "Male",
  ifelse(cohort[["GENDER"]] == "F", "Female", "Missing")
)
gender[is.na(gender)] <- "Missing"

tpn_raw <- suppressWarnings(as.numeric(cohort[["TPN_HYPERAL"]]))
tpn <- ifelse(tpn_raw == 0, "No", ifelse(tpn_raw == 1, "Yes", "Missing/invalid"))
tpn[is.na(tpn)] <- "Missing/invalid"

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
  section_row("Parent-reported ethnicity (PRE), v8 hierarchy"),
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

# -----------------------------------------------------------------------------
# Reference-style subgroup table: false-positive versus true-positive referrals
# -----------------------------------------------------------------------------

gestational_week <- floor(gestational_days / 7)
gestational_week[
  is.na(gestational_days) | gestational_days < 140 | gestational_days > 315
] <- NA_real_
gestational_category <- ifelse(
  is.na(gestational_week),
  "Unknown/invalid",
  ifelse(
    gestational_week > 42,
    ">42",
    ifelse(
      gestational_week == 42,
      "42",
      ifelse(
        gestational_week == 41,
        "41",
        ifelse(
          gestational_week >= 39,
          "39-40",
          ifelse(
            gestational_week >= 37,
            "37-38",
            ifelse(gestational_week >= 28, "28-36", "<28")
          )
        )
      )
    )
  )
)

birth_weight_category <- ifelse(
  is.na(birth_weight),
  "Unknown",
  ifelse(
    birth_weight > 5000,
    ">5000",
    ifelse(
      birth_weight >= 4001,
      "4001-5000",
      ifelse(
        birth_weight >= 3501,
        "3501-4000",
        ifelse(
          birth_weight >= 3001,
          "3001-3500",
          ifelse(
            birth_weight >= 2500,
            "2500-3000",
            ifelse(birth_weight >= 1000, "1000-2499", "<1000")
          )
        )
      )
    )
  )
)

collection_age_category <- ifelse(
  is.na(collection_age),
  "Unknown",
  ifelse(
    collection_age < 12,
    "<12",
    ifelse(
      collection_age < 24,
      "12-23",
      ifelse(
        collection_age <= 48,
        "24-48",
        ifelse(collection_age <= 168, "49-168", ">168")
      )
    )
  )
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

format_subgroup_count <- function(n, denominator) {
  sprintf("%d (%.1f%%)", n, 100 * n / denominator)
}

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
          false_positive_display = format_subgroup_count(fp_n, fp_denominator),
          true_positive_n = tp_n,
          true_positive_denominator = tp_denominator,
          true_positive_display = format_subgroup_count(tp_n, tp_denominator),
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

run_manifest <- data.frame(
  item = c(
    "run_date",
    "cohort_n",
    "join_method",
    "pre_hierarchy",
    "continuous_summary",
    "percentage_denominator",
    "tpn_missing_invalid_definition"
  ),
  value = c(
    format(Sys.Date(), "%Y-%m-%d"),
    "378",
    "Validated sample-ID join with false-positive prefix normalization",
    "Hispanic > Black > East Asian > South Asian > Middle Eastern > Native American > White; otherwise Other/Unknown",
    "Median (interquartile range); screening year uses median (range)",
    "N = 378 for categorical characteristics",
    "Blank/NA or code 998"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  run_manifest,
  file.path(table_dir, "cohort_characteristics_run_manifest.csv"),
  row.names = FALSE
)

report_lines <- c(
  "# Cohort characteristics table report",
  "",
  paste0("Run date: ", format(Sys.Date(), "%Y-%m-%d")),
  "",
  "## Source and cohort",
  "",
  paste0(
    "The table was generated from the canonical phenotype workbook and the exact ",
    "378 study IDs in the joint FAM file. Samples were joined by validated IDs; ",
    "the historical false-positive prefix difference was normalized explicitly."
  ),
  "",
  "## Summary decisions",
  "",
  "- Categorical values are n (%) using N = 378.",
  "- Birth weight and age at blood collection are median (IQR) because age at collection is right-skewed.",
  "- Available-case denominators are 375 for birth weight and 375 for age at blood collection.",
  "- PRE uses the prespecified hierarchy and reproduces the manuscript counts.",
  "- TPN missing/invalid combines blank or NA values with code 998; no raw values were changed.",
  "- The stratified table reports percentages within the false-positive (n=143) and true-positive (n=235) columns.",
  "- Gestational age was derived from GA_DAYS as completed weeks. Blank values and four records outside 140-315 days were classified as unknown/invalid.",
  "- Age-at-collection bins are nonoverlapping: 12-23 hours and 24-48 hours.",
  "- The exported source table contains aggregate values only and no sample identifiers.",
  "",
  "## Validation",
  "",
  "- Total cohort: 378; true positive: 235; false positive: 143.",
  "- Multiple PRE selections: 52.",
  "- PRE counts: Hispanic 205; White 83; East Asian 34; Black 23; Other/Unknown 14; Middle Eastern 9; South Asian 9; Native American 1.",
  "- Gender: male 214; female 161; missing 3.",
  "- TPN: no 303; yes 55; missing/invalid 20."
)
writeLines(report_lines, file.path(analysis_dir, "analysis_report.md"))
capture.output(sessionInfo(), file = file.path(analysis_dir, "sessionInfo.txt"))

cat("Wrote aggregate cohort table to", table_dir, "\n")
