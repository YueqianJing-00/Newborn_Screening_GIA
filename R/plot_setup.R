# Shared constants and export helpers for single-panel plotting scripts.
# This support file does not create a plot; every executable plotting script
# exports exactly one subplot.

script_file <- function() {
  normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
    mustWork = TRUE
  )
}

project_root_from_script <- function() {
  normalizePath(file.path(dirname(script_file()), "..", ".."), mustWork = TRUE)
}

results_root <- function(project_root) {
  normalizePath(
    Sys.getenv("HGG_RESULTS_DIR", file.path(project_root, "results")),
    mustWork = FALSE
  )
}

save_panel <- function(plot, directory, stem, width, height) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file.path(directory, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(directory, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 300,
    device = ragg::agg_png,
    bg = "white"
  )
}

ancestry_levels <- c("AMR", "AFR", "EUR", "SAS", "EAS")
ancestry_palette <- c(
  AMR = "#E69F00",
  AFR = "#D55E00",
  EUR = "#0072B2",
  SAS = "#009E73",
  EAS = "#CC79A7"
)

theme_panel <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

group_boundaries <- function(data, group_column, index_column) {
  data |>
    dplyr::group_by(.data[[group_column]]) |>
    dplyr::summarise(
      start = min(.data[[index_column]]),
      end = max(.data[[index_column]]),
      center = (start + end) / 2,
      .groups = "drop"
    )
}
