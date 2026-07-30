# Figure export helpers.

save_figure_pair <- function(plot, directory, stem, width, height) {
  paths <- c(
    pdf = file.path(directory, paste0(stem, ".pdf")),
    png = file.path(directory, paste0(stem, ".png"))
  )
  ggplot2::ggsave(
    paths[["pdf"]],
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
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
  paths
}
