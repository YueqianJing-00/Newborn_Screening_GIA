# Figure export helpers.

save_figure_pair <- function(plot, directory, stem, width, height) {
  # Use one stem for the editable vector and publication-resolution raster files.
  paths <- c(
    pdf = file.path(directory, paste0(stem, ".pdf")),
    png = file.path(directory, paste0(stem, ".png"))
  )

  # Cairo preserves vector text and shapes in the PDF export.
  ggplot2::ggsave(
    paths[["pdf"]],
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  # ragg renders an antialiased 300-dpi PNG for submission systems.
  ggplot2::ggsave(
    paths[["png"]],
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    device = ragg::agg_png,
    bg = "white"
  )
  invisible(paths)
}
