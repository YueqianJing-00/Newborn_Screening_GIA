# Helpers for matching study samples and labeling ADMIXTURE components.

canonical_sample_id <- function(x) {
  sub("^NBSfalsepos_", "", as.character(x))
}

validate_component_mapping <- function(summary_table, ancestry_columns) {
  largest_component <- ancestry_columns[
    max.col(as.matrix(summary_table[, ancestry_columns, drop = FALSE]))
  ]
  if (!identical(as.character(summary_table$SuperPop), largest_component)) {
    stop(
      "Ancestry components could not be matched to the reference populations.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

read_global_ancestry <- function(q_path, fam_path, psam_path, study_n = 378L) {
  fam <- data.table::fread(fam_path, header = FALSE)
  q <- data.table::fread(q_path, header = FALSE)
  psam <- data.table::fread(psam_path)

  if (nrow(fam) != nrow(q)) stop("FAM and Q row counts differ.")
  if (ncol(q) != 5L) stop("Expected five ADMIXTURE components.")

  reference_n <- nrow(fam) - study_n
  if (reference_n < 1L) stop("Joint dataset contains no reference samples.")

  psam_id <- grep("IID$", names(psam), value = TRUE)
  if (length(psam_id) != 1L) stop("Could not identify a unique IID column in the PSAM file.")
  data.table::setnames(psam, psam_id, "IID")

  q_columns <- paste0("V", seq_len(5L))
  reference <- data.table::data.table(
    IID = fam$V2[seq_len(reference_n)],
    q[seq_len(reference_n)]
  )
  reference <- merge(
    reference,
    psam[, .(IID, SuperPop)],
    by = "IID",
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(reference$SuperPop)) stop("Some reference IDs are missing from the PSAM file.")

  reference_means <- reference[
    , lapply(.SD, mean),
    by = SuperPop,
    .SDcols = q_columns
  ]
  ancestry_components <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  if (!setequal(reference_means$SuperPop, ancestry_components)) {
    stop("Unexpected reference superpopulations.")
  }

  mapping <- reference_means[, {
    values <- unlist(.SD)
    index <- which.max(values)
    .(
      q_column = names(values)[index],
      mean_reference_membership = values[index]
    )
  }, by = SuperPop, .SDcols = q_columns]
  if (data.table::uniqueN(mapping$q_column) != 5L) {
    stop("ADMIXTURE component mapping is not one-to-one.")
  }
  if (any(mapping$mean_reference_membership < 0.95)) {
    stop("ADMIXTURE components are not well separated in the reference samples.")
  }

  study_rows <- (reference_n + 1L):nrow(fam)
  study_ids <- as.character(fam$V2[study_rows])
  unprefixed <- grepl("^p04w", study_ids)
  if (sum(unprefixed) != 48L) stop("Unexpected number of unprefixed study IDs.")
  study_ids[unprefixed] <- paste0("NBSfalsepos_", study_ids[unprefixed])
  if (anyDuplicated(study_ids)) stop("Duplicate study IDs after normalization.")

  study <- data.table::as.data.table(q[study_rows])
  q_to_label <- setNames(mapping$SuperPop, mapping$q_column)
  data.table::setnames(study, q_columns, unname(q_to_label[q_columns]))
  study[, raw_id := study_ids]
  data.table::setcolorder(study, c("raw_id", ancestry_components))

  maximum_sum_error <- max(
    abs(rowSums(study[, ..ancestry_components]) - 1)
  )
  if (maximum_sum_error > 5e-4) {
    stop("Study ancestry proportions do not sum to one within tolerance.")
  }

  list(
    study = study,
    mapping = mapping[order(SuperPop)],
    reference_means = reference_means[order(SuperPop)],
    reference_n = reference_n,
    maximum_sum_error = maximum_sum_error
  )
}
