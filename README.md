# Newborn screening and GIA analysis

R code for the manuscript analyses comparing parent-reported ethnicity (PRE)
with genetically inferred ancestry (GIA) in newborn-screening referrals.


## Requirements

The analysis was run with R 4.5.1 and the following packages:

```r
install.packages(c(
  "data.table", "dplyr", "ggplot2", "pROC",
  "ragg", "randomForest", "readxl", "scales", "tidyr"
))
```

## Upstream ancestry estimation

The manuscript scripts consume precomputed ADMIXTURE ancestry proportions. The
code under [`GIA/`](GIA/) prepares the 1000 Genomes
reference, selects the reference panel, harmonizes reference and study
genotypes, and estimates ancestry in study samples with supervised ADMIXTURE.
Those workflows keep genotype data and participant-level outputs outside Git
and keep editable paths in configuration blocks at the top of each script.


## Data analysis

| Step | Script | Purpose |
| --- | --- | --- |
| 1 | `analysis/data_analysis/00_descriptive_analysis.R` | Prepare the 378-newborn ancestry/PRE dataset and descriptive source tables |
| 2 | `analysis/data_analysis/01_cohort_characteristics.R` | Generate cohort characteristics and exclusion summaries |
| 3 | `analysis/data_analysis/02_mma_random_forest.R` | Fit the four MMA random-forest models and save performance and importance tables |
| 4 | `analysis/data_analysis/03_prediction_shift_analysis.R` | Analyze individual prediction changes after adding GIA |

