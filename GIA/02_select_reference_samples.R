#!/usr/bin/env Rscript

# Edit these paths and parameters before running.
reference_output_dir <- "/secure/work/global_reference"
k <- 5L
q_path <- file.path(reference_output_dir, sprintf("reference_impact.%d.Q", k))
fam_path <- file.path(reference_output_dir, "reference_impact.fam")
psam_path <- "/path/to/1000g_phase3/all_phase3.psam"
keep_path <- file.path(reference_output_dir, "reference_selected.keep")
labels_path <- file.path(reference_output_dir, "reference_selected.labels.tsv")
threshold <- 0.80

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

write.table(
  selected_keep,
  keep_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
write.table(
  selected_labels,
  labels_path,
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
