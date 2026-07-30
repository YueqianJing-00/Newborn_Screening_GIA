# Newborn screening and genetic ancestry analysis

This repository contains the R code used to compare parent-reported ethnicity
(PRE) with genetically inferred ancestry (GIA) in newborn-screening referrals.
It also contains the MMA random-forest analysis reported in the manuscript.
Participant data and individual-level results are not included.

The analyses cover:

1. PRE-GIA concordance in 378 sequenced screen-positive newborns.
2. Cohort characteristics and Figures 1 and 2.
3. Random-forest models in 117 MMA screen-positive newborns after excluding
   newborns receiving total parenteral nutrition (TPN).
4. Model validation and Figures 3 and 4.

## Project layout

```text
analysis/   Scripts run in manuscript order
R/          Shared data, statistics, plotting, and path helpers
config/     Input-file and software-version records
data/raw/   Local input data (ignored by Git)
docs/       Workflow and data-release notes
resources/  Figure layout template
results/    Generated files (ignored by Git)
run_all.R   Pipeline entry point
```

The role of each numbered script is listed in
[`docs/ANALYSIS_MAP.md`](docs/ANALYSIS_MAP.md).

## Requirements

The manuscript analysis used R 4.5.1. Package versions are recorded in
[`config/software_versions.csv`](config/software_versions.csv), and the required
packages are listed in [`DESCRIPTION`](DESCRIPTION).

Place the six input files listed in
[`config/input_manifest.csv`](config/input_manifest.csv) in `data/raw/`. To keep
the data elsewhere, set:

```sh
export HGG_DATA_DIR="/secure/path/to/input/files"
```

Commands below assume the repository root is the working directory.

```sh
# Check file structure, R syntax, and privacy rules
Rscript run_all.R --mode=check

# Run all analyses in order
Rscript run_all.R --mode=full
```

Results are written to `results/` unless `HGG_RESULTS_DIR` is set. The full model
is the slowest step: it uses 100 repeated stratified 10-fold cross-validation
runs, 1,000 trees per forest, and 2,000 subject-level bootstrap samples.

Within each training fold, the model ranks 40 metabolite candidates by the
absolute distance of their univariate AUC from 0.5 and selects the top 10. FC and
C3/C2 are included in this candidate set but are not forced into the model.

## Data release

Raw data, fold assignments, subject-level predictions, and individual ancestry
profiles must remain outside Git. See [`docs/PRIVACY.md`](docs/PRIVACY.md) before
adding any generated result to the repository.
