# Analysis workflow

`run_all.R` executes the numbered scripts in order. Scripts 00, 01, 04, and 07
read the controlled-access inputs; the figure scripts read tables created by an
earlier step.

| Step | Script | Purpose | Main output |
|---:|---|---|---|
| 00 | `analysis/00_descriptive_analysis.R` | Prepare the 378-subject PRE-GIA data and descriptive statistics | Descriptive source tables |
| 01 | `analysis/01_cohort_characteristics.R` | Summarize cohort characteristics | Cohort tables |
| 02 | `analysis/02_figure1_pre_gia_concordance.R` | Plot PRE-GIA concordance | Figure 1 |
| 03 | `analysis/03_figure2_entropy_multiple_pre.R` | Analyze ancestry entropy and multiple PRE selections | Figure 2 |
| 04 | `analysis/04_mma_random_forest.R` | Fit and evaluate the four MMA models | Model estimates and source tables |
| 05 | `analysis/05_validate_mma_random_forest.R` | Recalculate model metrics and check output consistency | Validation report |
| 06 | `analysis/06_figure3_model_performance.R` | Plot performance and permutation importance | Figure 3 |
| 07 | `analysis/07_figure4_prediction_shifts.R` | Calculate individual prediction shifts after adding GIA | Aggregate shift summaries |
| 08 | `analysis/08_figure4_plot.R` | Plot the final prediction-shift panels | Figure 4 |

Shared functions live in `R/`:

- `project_setup.R` handles paths, package checks, and file checksums.
- `ancestry_helpers.R` matches study IDs and labels ADMIXTURE components.
- `pre_helpers.R` applies the PRE harmonization and assignment hierarchy.
- `statistical_helpers.R` contains agreement and prediction metrics.
- `plot_helpers.R` exports paired PDF and PNG figures.

Local-ancestry and gestational-age sensitivity analyses are not part of the
reported workflow and are therefore not included here.
