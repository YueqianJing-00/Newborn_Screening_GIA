#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Could not determine run_all.R path from --file.")
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

root <- dirname(get_script_path())
old_working_directory <- setwd(root)
on.exit(setwd(old_working_directory), add = TRUE)
args <- commandArgs(trailingOnly = TRUE)
mode_arg <- grep("^--mode=", args, value = TRUE)

if (length(mode_arg) != 1L || length(args) != 1L) {
  stop("Usage: Rscript run_all.R --mode=check|full", call. = FALSE)
}

mode <- sub("^--mode=", "", mode_arg)
if (!mode %in% c("check", "full")) {
  stop("--mode must be check or full.", call. = FALSE)
}

run_step <- function(label, relative_script, script_args = character()) {
  script <- file.path(root, relative_script)
  if (!file.exists(script)) stop("Missing pipeline script: ", relative_script)

  message("\n[", label, "] ", relative_script)
  command_args <- c("--vanilla", relative_script, script_args)
  status <- system2(file.path(R.home("bin"), "Rscript"), command_args)
  if (!identical(status, 0L)) {
    stop("Pipeline step failed (exit ", status, "): ", relative_script, call. = FALSE)
  }
  invisible(TRUE)
}

if (mode == "check") {
  run_step("release checks", "analysis/99_check_release.R")
  quit(save = "no", status = 0L)
}

steps <- list(
  c("descriptive analysis", "analysis/00_descriptive_analysis.R"),
  c("cohort characteristics", "analysis/01_cohort_characteristics.R"),
  c("Figure 1", "analysis/02_figure1_pre_gia_concordance.R"),
  c("Figure 2", "analysis/03_figure2_entropy_multiple_pre.R"),
  c("MMA random forest", "analysis/04_mma_random_forest.R"),
  c("model validation", "analysis/05_validate_mma_random_forest.R"),
  c("Figure 3", "analysis/06_figure3_model_performance.R"),
  c("Figure 4 analysis", "analysis/07_figure4_prediction_shifts.R"),
  c("Figure 4 plot", "analysis/08_figure4_plot.R")
)

for (step in steps) run_step(step[[1]], step[[2]])

message("\nPipeline completed. Outputs are under: ",
        Sys.getenv("HGG_RESULTS_DIR", unset = file.path(root, "results")))
