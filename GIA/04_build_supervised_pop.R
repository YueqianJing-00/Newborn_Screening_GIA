#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 04_build_supervised_pop.R JOINT.fam REFERENCE.labels.tsv OUTPUT.pop",
    call. = FALSE
  )
}

fam_path <- args[[1L]]
labels_path <- args[[2L]]
output_path <- args[[3L]]

missing_paths <- c(fam_path, labels_path)[
  !file.exists(c(fam_path, labels_path))
]
if (length(missing_paths) > 0L) {
  stop(
    paste("Missing input file(s):", paste(missing_paths, collapse = ", ")),
    call. = FALSE
  )
}

fam <- read.table(
  fam_path,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
labels <- read.table(
  labels_path,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  col.names = c("IID", "SuperPop")
)

if (ncol(fam) < 2L) stop("FAM must contain FID and IID columns.", call. = FALSE)
if (ncol(labels) != 2L) stop("Reference labels must have two columns.", call. = FALSE)
if (anyDuplicated(fam[[2L]])) stop("Joint FAM contains duplicate IIDs.", call. = FALSE)
if (anyDuplicated(labels$IID)) stop("Reference labels contain duplicate IIDs.", call. = FALSE)

expected_groups <- c("AFR", "AMR", "EAS", "EUR", "SAS")
if (!setequal(unique(labels$SuperPop), expected_groups)) {
  stop("Reference labels must cover AFR, AMR, EAS, EUR, and SAS.", call. = FALSE)
}
if (!all(labels$IID %in% fam[[2L]])) {
  stop("Some labeled reference IIDs are absent from the joint FAM.", call. = FALSE)
}

label_index <- match(as.character(fam[[2L]]), labels$IID)
population <- rep("-", nrow(fam))
population[!is.na(label_index)] <- labels$SuperPop[label_index[!is.na(label_index)]]

if (!any(population == "-")) stop("Joint FAM contains no study samples.", call. = FALSE)
if (sum(population != "-") != nrow(labels)) {
  stop("Reference label count does not match the joint FAM.", call. = FALSE)
}

output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write.table(
  population,
  output_path,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message(
  sprintf(
    "Wrote %d reference labels and %d study placeholders.",
    sum(population != "-"),
    sum(population == "-")
  )
)
