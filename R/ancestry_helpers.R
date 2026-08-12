# Helpers for matching study samples and labeling ADMIXTURE components.

canonical_sample_id <- function(x) {
  # Remove the cohort prefix used only for false-positive sample files.
  sub("^NBSfalsepos_", "", as.character(x))
}

validate_component_mapping <- function(summary_table, ancestry_columns) {
  # The largest mean component must match each known reference population.
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
  # Read ADMIXTURE proportions, their sample order, and reference metadata.
  fam <- data.table::fread(fam_path, header = FALSE)
  q <- data.table::fread(q_path, header = FALSE)
  psam <- data.table::fread(psam_path)
  ancestry_components <- c("AFR", "AMR", "EAS", "EUR", "SAS")

  reference_n <- nrow(fam) - study_n
  valid_dimensions <-
    nrow(fam) == nrow(q) &&
    ncol(q) == length(ancestry_components) &&
    reference_n > 0L
  if (!valid_dimensions) {
    stop("FAM and five-component Q dimensions do not match.")
  }

  data.table::setnames(psam, sub("^#", "", names(psam)))

  q_columns <- paste0("V", seq_along(ancestry_components))

  # Join reference rows to known superpopulations without changing FAM order.
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

  # Average each raw Q column within the five reference superpopulations.
  reference_means <- reference[
    , lapply(.SD, mean),
    by = SuperPop,
    .SDcols = q_columns
  ]

  # Name each Q column for the population where its mean membership is largest.
  mapping <- reference_means[, {
    values <- unlist(.SD)
    index <- which.max(values)
    .(
      q_column = names(values)[index],
      mean_reference_membership = values[index]
    )
  }, by = SuperPop, .SDcols = q_columns]

  mapping_is_clear <-
    setequal(reference_means$SuperPop, ancestry_components) &&
    data.table::uniqueN(mapping$q_column) == length(ancestry_components) &&
    all(mapping$mean_reference_membership >= 0.95)
  if (!mapping_is_clear) {
    stop("Reference populations do not define five distinct ADMIXTURE components.")
  }

  # Apply the inferred component names to the study rows.
  study_rows <- (reference_n + 1L):nrow(fam)
  study_ids <- as.character(fam$V2[study_rows])
  unprefixed <- grepl("^p04w", study_ids)
  study_ids[unprefixed] <- paste0("NBSfalsepos_", study_ids[unprefixed])
  if (anyDuplicated(study_ids)) stop("Duplicate study IDs after normalization.")

  study <- data.table::as.data.table(q[study_rows])
  q_to_label <- setNames(mapping$SuperPop, mapping$q_column)
  data.table::setnames(study, q_columns, unname(q_to_label[q_columns]))
  study[, raw_id := study_ids]
  data.table::setcolorder(study, c("raw_id", ancestry_components))

  maximum_sum_error <- max(abs(rowSums(study[, ..ancestry_components]) - 1))
  if (maximum_sum_error > 5e-4) stop(
    "Study ancestry proportions do not sum to one within tolerance."
  )

  list(
    study = study,
    mapping = mapping[order(SuperPop)],
    reference_means = reference_means[order(SuperPop)],
    reference_n = reference_n,
    maximum_sum_error = maximum_sum_error
  )
}
