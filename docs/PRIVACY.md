# Privacy and disclosure rules

This project combines newborn-screening phenotypes, reported population descriptors, ancestry estimates, and model predictions. Deidentification alone does not make individual-level derivatives appropriate for unrestricted release.

## Never commit

- Raw phenotype workbooks or files containing sample identifiers.
- Joint FAM/Q/PSAM inputs or other files that expose study ordering or ancestry profiles.
- Cross-validation fold assignments.
- Per-repeat or subject-averaged out-of-fold predictions.
- Subject bootstrap weights.
- Modeling datasets, even when direct identifiers have been removed.
- Individual profile or case-pattern files marked `restricted_internal`.

## Potentially releasable after review

- Final figure PDFs/PNGs/TIFFs already approved for the manuscript.
- Aggregate cohort counts.
- Aggregate performance, paired-effect, metabolite-selection-frequency, and permutation-importance summaries.
- Aggregate figure source tables that cannot be traced to individuals or very small cells.

## Before each public push

Run:

```sh
Rscript run_all.R --mode=check
git status --short
git diff --cached --name-only
```

Review every staged file. A passing automated check is a safeguard, not governance approval.

