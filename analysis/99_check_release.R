#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Could not determine script path from --file.")
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

script_path <- get_script_path()
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

required_files <- c(
  "DESCRIPTION", ".gitignore", "README.md", "Makefile", "run_all.R",
  "R/project_paths.R", "config/input_manifest.csv",
  "config/software_versions.csv", "data/raw/README.md",
  "resources/figure1_layout.tex.in",
  "results/README.md", "docs/ANALYSIS_MAP.md", "docs/PRIVACY.md",
  "docs/RELEASE_CHECKLIST.md",
  sprintf("analysis/%02d_%s.R", 0:8, c(
    "descriptive_analysis", "cohort_characteristics",
    "figure1_pre_gia_concordance", "figure2_entropy_multiple_pre",
    "mma_random_forest", "validate_mma_random_forest",
    "figure3_model_performance", "figure4_prediction_shifts",
    "figure4_plot"
  ))
)

missing_required <- required_files[!file.exists(file.path(root, required_files))]
if (length(missing_required) > 0L) {
  stop("Missing required release files: ", paste(missing_required, collapse = ", "))
}

r_files <- list.files(
  root,
  pattern = "\\.[Rr]$",
  recursive = TRUE,
  full.names = TRUE
)
parse_errors <- character()
for (file in r_files) {
  tryCatch(
    parse(file = file, keep.source = FALSE),
    error = function(e) {
      parse_errors <<- c(
        parse_errors,
        paste0(sub(paste0("^", root, "/"), "", file), ": ", conditionMessage(e))
      )
    }
  )
}
if (length(parse_errors) > 0L) {
  stop("R syntax failures:\n", paste(parse_errors, collapse = "\n"))
}

all_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
all_files <- all_files[file.info(all_files)$isdir == FALSE]
relative_files <- sub(paste0("^", root, "/"), "", all_files)

public_scan <- !grepl("^(data/raw|results|\\.git)(/|$)", relative_files)
unsafe_names <- relative_files[public_scan & (
  grepl("\\.(xlsx?|fam|psam|Q|rds|RData|rdata)$", relative_files, ignore.case = TRUE) |
    grepl(
      "restricted_internal|fold_assignments|oof_predictions|subject_bootstrap_weights|modeling_dataset|subject_level",
      basename(relative_files),
      ignore.case = TRUE
    )
)]
if (length(unsafe_names) > 0L) {
  stop("Potentially restricted files outside ignored data/results directories: ",
       paste(unsafe_names, collapse = ", "))
}

text_extensions <- "\\.(R|md|csv|yml|yaml|txt|DESCRIPTION|gitignore)$"
text_files <- all_files[public_scan & (
  grepl(text_extensions, basename(all_files), ignore.case = TRUE) |
    basename(all_files) %in% c("DESCRIPTION", ".gitignore", "Makefile")
)]
text_files <- setdiff(normalizePath(text_files), normalizePath(script_path))
forbidden_patterns <- c(
  local_absolute_path = "/Users/",
  retired_tpn_term = "TPN-unexposed",
  legacy_project_path = "analysis/GWAS-metabolite"
)
text_violations <- character()
for (file in text_files) {
  lines <- readLines(file, warn = FALSE)
  for (pattern_name in names(forbidden_patterns)) {
    hits <- grep(forbidden_patterns[[pattern_name]], lines, fixed = TRUE)
    if (length(hits) > 0L) {
      text_violations <- c(
        text_violations,
        paste0(
          sub(paste0("^", root, "/"), "", file), ":",
          paste(hits, collapse = ","), " (", pattern_name, ")"
        )
      )
    }
  }
}
if (length(text_violations) > 0L) {
  stop("Nonportable or retired text detected:\n", paste(text_violations, collapse = "\n"))
}

ignore_lines <- readLines(file.path(root, ".gitignore"), warn = FALSE)
required_ignore_rules <- c(
  "data/raw/**", "results/**", "**/*restricted_internal*",
  "**/fold_assignments.csv", "**/oof_predictions*.csv",
  "**/subject_level*.csv"
)
missing_ignore_rules <- setdiff(required_ignore_rules, ignore_lines)
if (length(missing_ignore_rules) > 0L) {
  stop("Missing required .gitignore safeguards: ",
       paste(missing_ignore_rules, collapse = ", "))
}

description <- read.dcf(file.path(root, "DESCRIPTION"))
if (!identical(unname(description[1, "Version"]), "0.1.0")) {
  stop("Unexpected release version in DESCRIPTION.")
}

cat("PASS: required release structure is complete\n")
cat("PASS:", length(r_files), "R files parse successfully\n")
cat("PASS: no obvious restricted files occur outside ignored directories\n")
cat("PASS: no local absolute paths or retired TPN terminology were found\n")
cat("PASS: required Git privacy safeguards are present\n")
