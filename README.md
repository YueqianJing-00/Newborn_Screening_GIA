# HGG Advances newborn-screening analysis code

This directory is a GitHub-ready, privacy-conscious release of the analysis code used by the current manuscript. It separates the final analysis from superseded development versions and keeps controlled-access participant data and subject-level outputs outside version control.

The code reproduces:

1. PRE–GIA descriptive analyses in 378 sequenced screen-positive newborns.
2. Cohort characteristics and manuscript Figures 1 and 2.
3. The leakage-controlled random-forest analysis in 117 MMA screen-positive newborns after excluding newborns receiving TPN.
4. Manuscript Figures 3 and 4 and independent checks of the model outputs.

## Repository structure

```text
analysis/       Numbered analysis and figure scripts
R/              Shared path helpers
config/         Input and software manifests
data/raw/       Local controlled-access inputs; ignored by Git
docs/           Analysis map, privacy rules, and release checklist
results/        Generated outputs; ignored by Git by default
run_all.R       Ordered pipeline runner
```

## Requirements

The locked analyses used R 4.5.1. Direct package versions are recorded in `config/software_versions.csv` and minimum versions are listed in `DESCRIPTION`.

Place the six expected input files in `data/raw/`, or point to an existing controlled-access directory:

```sh
export HGG_DATA_DIR="/secure/path/to/input/files"
```

The required filenames and their roles are listed in `config/input_manifest.csv`. Raw inputs are deliberately not included.

## Running the code

Run static, privacy, and syntax checks:

```sh
Rscript run_all.R --mode=check
```

Run the complete pipeline in manuscript order:

```sh
Rscript run_all.R --mode=full
```

By default, outputs are written below `results/`. To use another private output location, set `HGG_RESULTS_DIR` before running.

The full model uses 100 repeated stratified 10-fold cross-validation runs, 1,000 trees per forest, training-fold selection of 10 from 40 metabolite candidates, and 2,000 outcome-stratified subject bootstrap samples. It is computationally heavier than the descriptive scripts.

## Data and privacy

Do not commit the phenotype workbook, FAM/Q/PSAM files, fold assignments, subject-level predictions, bootstrap weights, or any file containing `restricted_internal` in its name. Generated results remain ignored until each aggregate file is reviewed for public release. See `docs/PRIVACY.md`.

## Provenance

The scientific logic was copied from the locked manuscript scripts. Changes in this release are organizational: portable input/output paths, clearer filenames, removal of legacy version labels, and Git/privacy safeguards. The original-to-release mapping is in `docs/ANALYSIS_MAP.md`.

No Git repository or GitHub remote has been created yet. Complete `docs/RELEASE_CHECKLIST.md` before the first public push.

