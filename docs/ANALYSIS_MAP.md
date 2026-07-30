# Analysis map

The numbered scripts are the final code path for the current manuscript. Run them in order with `run_all.R`.

| Order | Manuscript component | Release script | Primary dependency | Original locked source |
|---:|---|---|---|---|
| 00 | Descriptive PRE–GIA statistics and reusable source tables | `analysis/00_descriptive_analysis.R` | Controlled-access inputs | `analysis/publication_v8/descriptive/generate_descriptive_figures.R` |
| 01 | Cohort-characteristics table | `analysis/01_cohort_characteristics.R` | Controlled-access inputs | `analysis/publication_v8/cohort_characteristics_20260717/generate_cohort_characteristics.R` |
| 02 | Figure 1: PRE–GIA concordance | `analysis/02_figure1_pre_gia_concordance.R` | Script 00 tables | `analysis/publication_v8/concordance_visualization_20260720/figure1_reference_layout_20260722/generate_figure1_reference_layout.R` |
| 03 | Figure 2: entropy and multiple PRE profiles | `analysis/03_figure2_entropy_multiple_pre.R` | Script 00 tables | `analysis/publication_v8/entropy_multiple_pre_v3_20260722/generate_entropy_multiple_pre_v3.R` |
| 04 | Main MMA random-forest analysis | `analysis/04_mma_random_forest.R` | Controlled-access inputs | `analysis/publication_v8/mma_tpn0_top10_fc_ratio_20260724/run_mma_tpn0_top10_fc_ratio.R` |
| 05 | Independent model-output validation | `analysis/05_validate_mma_random_forest.R` | Script 04 outputs | `analysis/publication_v8/mma_tpn0_top10_fc_ratio_20260724/validate_outputs.R` |
| 06 | Figure 3: performance and importance | `analysis/06_figure3_model_performance.R` | Script 04 outputs | `analysis/publication_v8/rf_figure3_v6_20260724/generate_figure3_rf_comparison_v6.R` |
| 07 | Figure 4 prediction-shift analysis | `analysis/07_figure4_prediction_shifts.R` | Inputs and script 04 outputs | `analysis/publication_v8/rf_figure4_case_patterns_20260722/generate_figure4_gia_case_patterns.R` |
| 08 | Final two-panel Figure 4 | `analysis/08_figure4_plot.R` | Script 07 outputs | `analysis/publication_v8/rf_figure4_case_patterns_20260722/redraw_figure4_two_panel_nature_style.R` |

## Deliberately excluded from the main release

- Earlier Figure 3 versions (`rf_figure3_20260722` through `rf_figure3_v5_20260722`) are superseded.
- The earlier fixed-panel four-model analysis is superseded by training-fold metabolite selection.
- Specificity-tuning experiments are exploratory and are not the reported model.
- Local-ancestry code is not included because local-ancestry results were removed from the current manuscript.
- Gestational-age complete-case and imputation analyses from 2026-07-28 are retained privately as sensitivity work but are not part of the current manuscript-v11 analysis path.
- Manuscript-building and reference-audit Python scripts are editorial utilities, not scientific analysis code.

If the manuscript changes, update this map before releasing the repository rather than adding every historical script.

