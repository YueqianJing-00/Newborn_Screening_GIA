get_release_root <- function(script_path) {
  root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
  marker <- file.path(root, "DESCRIPTION")
  if (!file.exists(marker)) {
    stop("Could not locate the code-release root; missing DESCRIPTION at: ", marker)
  }
  root
}

get_release_paths <- function(script_path) {
  root <- get_release_root(script_path)

  data_dir <- Sys.getenv("HGG_DATA_DIR", unset = file.path(root, "data", "raw"))
  results_dir <- Sys.getenv("HGG_RESULTS_DIR", unset = file.path(root, "results"))

  list(
    root = root,
    data = normalizePath(data_dir, mustWork = FALSE),
    results = normalizePath(results_dir, mustWork = FALSE)
  )
}

relative_to_release <- function(path, root) {
  normalized <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
  if (startsWith(normalized, prefix)) {
    substring(normalized, nchar(prefix) + 1L)
  } else {
    file.path("external", basename(normalized))
  }
}
