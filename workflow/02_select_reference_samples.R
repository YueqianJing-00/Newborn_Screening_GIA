#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L || length(args) > 5L) {
  stop(
    paste(
      "Usage: 02_select_reference_samples.R",
      "REFERENCE.Q REFERENCE.fam REFERENCE.psam OUTPUT_PREFIX [THRESHOLD]"
    ),
    call. = FALSE
  )
}

q_path <- args[[1L]]
fam_path <- args[[2L]]
psam_path <- args[[3L]]
output_prefix <- args[[4L]]
threshold <- if (length(args) == 5L) as.numeric(args[[5L]]) else 0.80

if (!is.finite(threshold) || threshold < 0 || threshold >= 1) {
  stop("THRESHOLD must be a number in [0, 1).", call. = FALSE)
}

input_paths <- c(q_path, fam_path, psam_path)
missing_paths <- input_paths[!file.exists(input_paths)]
if (length(missing_paths) > 0L) {
  stop(
    paste("Missing input file(s):", paste(missing_paths, collapse = ", ")),
    call. = FALSE
  )
}

q <- as.matrix(read.table(q_path, header = FALSE, check.names = FALSE))
storage.mode(q) <- "double"
fam <- read.table(
  fam_path,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
psam <- read.table(
  psam_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)
names(psam) <- sub("^#", "", names(psam))

if (ncol(q) != 5L) stop("Expected five ancestry components in Q.", call. = FALSE)
if (ncol(fam) < 2L) stop("FAM must contain FID and IID columns.", call. = FALSE)
if (nrow(q) != nrow(fam)) stop("Q and FAM row counts differ.", call. = FALSE)
if (!all(c("IID", "SuperPop") %in% names(psam))) {
  stop("PSAM must contain IID and SuperPop columns.", call. = FALSE)
}
if (anyDuplicated(fam[[2L]])) stop("FAM contains duplicate IIDs.", call. = FALSE)
if (anyDuplicated(psam$IID)) stop("PSAM contains duplicate IIDs.", call. = FALSE)
if (any(!is.finite(q)) || any(q < 0) || any(q > 1)) {
  stop("Q must contain finite ancestry proportions between zero and one.", call. = FALSE)
}
if (max(abs(rowSums(q) - 1)) > 1e-3) {
  stop("Q rows do not sum to one within tolerance.", call. = FALSE)
}

metadata_index <- match(as.character(fam[[2L]]), as.character(psam$IID))
if (anyNA(metadata_index)) stop("Some FAM IIDs are absent from PSAM.", call. = FALSE)

selected <- apply(q, 1L, max) > threshold
if (!any(selected)) stop("No reference samples passed the threshold.", call. = FALSE)

selected_rows <- which(selected)
selected_labels <- data.frame(
  IID = as.character(fam[[2L]][selected_rows]),
  SuperPop = as.character(psam$SuperPop[metadata_index[selected_rows]]),
  stringsAsFactors = FALSE
)

expected_groups <- c("AFR", "AMR", "EAS", "EUR", "SAS")
if (!setequal(unique(selected_labels$SuperPop), expected_groups)) {
  stop("Selected references do not cover AFR, AMR, EAS, EUR, and SAS.", call. = FALSE)
}

historical_group_order <- c("AMR", "AFR", "EUR", "SAS", "EAS")
selected_order <- order(
  match(selected_labels$SuperPop, historical_group_order),
  selected_rows
)
selected_rows <- selected_rows[selected_order]
selected_labels <- selected_labels[selected_order, , drop = FALSE]
selected_keep <- fam[selected_rows, 1:2, drop = FALSE]

output_dir <- dirname(output_prefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.table(
  selected_keep,
  paste0(output_prefix, ".keep"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
write.table(
  selected_labels,
  paste0(output_prefix, ".labels.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message(
  sprintf(
    "Selected %d of %d reference samples with max(Q) > %.2f.",
    sum(selected),
    length(selected),
    threshold
  )
)
