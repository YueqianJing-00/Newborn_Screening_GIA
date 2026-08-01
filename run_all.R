#!/usr/bin/env Rscript

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "R", "project_setup.R"))
root <- find_project_root(dirname(script_path))
old_working_directory <- setwd(root)
on.exit(setwd(old_working_directory), add = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) stop("Usage: Rscript run_all.R", call. = FALSE)

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
