# Paths and small setup helpers shared by the command-line scripts.

current_script <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Could not determine the current script from --file.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop("Could not find the project root (DESCRIPTION is missing).", call. = FALSE)
}

project_paths <- function(script = current_script()) {
  root <- find_project_root(dirname(script))
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

sha256_files <- function(paths) {
  vapply(
    paths,
    digest::digest,
    character(1),
    file = TRUE,
    algo = "sha256"
  )
}

is_run_name <- function(x) {
  length(x) == 1L && !is.na(x) && grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", x)
}
