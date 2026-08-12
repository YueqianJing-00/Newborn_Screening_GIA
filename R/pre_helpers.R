# Parent-reported ethnicity (PRE) harmonization used across analyses.

east_asian_pre_labels <- c(
  "Japanese", "Chinese", "Laos", "Korean", "Vietnamese", "Filipino"
)
south_asian_pre_labels <- "Asian East Indian"

normalize_pre_values <- function(
    values,
    east_asian_label = "EAS",
    south_asian_label = "SAS") {
  # Drop empty entries and collapse detailed Asian labels into analysis groups.
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  values[values %in% east_asian_pre_labels] <- east_asian_label
  values[values %in% south_asian_pre_labels] <- south_asian_label
  unique(values)
}

assign_pre <- function(
    values,
    hierarchy = c(
      "Hispanic", "Black", "EAS", "SAS", "Middle Eastern",
      "Native American", "White"
    ),
    east_asian_label = "EAS",
    south_asian_label = "SAS") {
  # Apply the prespecified hierarchy when more than one category is reported.
  values <- normalize_pre_values(
    values,
    east_asian_label = east_asian_label,
    south_asian_label = south_asian_label
  )

  for (category in setdiff(hierarchy, "White")) {
    if (category %in% values) return(category)
  }

  # White is assigned only when it is the sole retained response.
  if (identical(values, "White")) return("White")
  "Other/Unknown"
}

harmonize_pre_value <- function(
    value,
    east_asian_label = "EAS",
    south_asian_label = "SAS") {
  # Convert one raw response to the categories retained in summaries and plots.
  value <- normalize_pre_values(
    value,
    east_asian_label = east_asian_label,
    south_asian_label = south_asian_label
  )
  if (!length(value)) return(NA_character_)

  retained <- c(
    "Hispanic", "White", "Black", east_asian_label, south_asian_label,
    "Middle Eastern", "Native American"
  )
  if (value[[1L]] %in% retained) value[[1L]] else "Other/Unknown"
}
