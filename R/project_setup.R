# Paths and small setup helpers shared by the command-line scripts.

current_script <- function() {
  # Rscript supplies the active file as --file=<path>.
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )[[1L]]
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

find_project_root <- function(start = getwd()) {
  # Walk upward until the repository's analysis and helper folders are found.
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    is_project <- file.exists(file.path(path, "README.md")) &&
      dir.exists(file.path(path, "analysis")) &&
      dir.exists(file.path(path, "R"))
    if (is_project) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop("Could not find the project root.", call. = FALSE)
}

project_paths <- function(script = current_script()) {
  root <- find_project_root(dirname(script))

  # Environment variables can redirect controlled data and generated results.
  data_dir <- Sys.getenv("HGG_DATA_DIR", unset = file.path(root, "data", "raw"))
  results_dir <- Sys.getenv("HGG_RESULTS_DIR", unset = file.path(root, "results"))

  list(
    root = root,
    data = normalizePath(data_dir, mustWork = FALSE),
    results = normalizePath(results_dir, mustWork = FALSE)
  )
}

relative_to_project <- function(path, root) {
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)

  # Keep project paths reproducible and mask external directory structure.
  if (startsWith(path, prefix)) {
    substring(path, nchar(prefix) + 1L)
  } else {
    file.path("external", basename(path))
  }
}

require_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(packages)
}

require_files <- function(paths, label = "input") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Missing required ", label, if (length(missing) == 1L) ": " else "s: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(paths)
}

make_directories <- function(...) {
  # Create every requested output directory, including missing parents.
  paths <- unlist(list(...), use.names = FALSE)
  vapply(
    paths,
    dir.create,
    logical(1),
    recursive = TRUE,
    showWarnings = FALSE
  )
  invisible(paths)
}

is_run_name <- function(x) {
  length(x) == 1L && !is.na(x) && grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", x)
}
