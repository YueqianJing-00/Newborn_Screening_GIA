# Ancestry estimation workflows

These scripts generate the global-ancestry files used by the manuscript analyses. They replace exploratory notebooks and cluster jobs that contained private paths or participant identifiers. Keep controlled genotype data and all generated files outside the repository.

## Software

The historical analyses used PLINK 1.90b6.21, PLINK 2.00a3.7LM, and ADMIXTURE 1.3.0. Install these programs and R before running the workflow. The scripts assume that the configured input files are ready and do not perform software or file preflight checks.

Copy the example configuration and edit the paths:

```bash
cp GIA/config.example.env GIA/config.env
set -a
source GIA/config.env
set +a
```

## Global ancestry

### 1. Estimate ancestry in 1000 Genomes

The reference workflow retains autosomal, strict biallelic A/C/G/T SNPs with MAF at least 1%, applies LD pruning with a 1,000-variant window, 100-variant step, and r2 threshold of 0.2, then runs unsupervised ADMIXTURE at K=5.

```bash
WORK_DIR="$GLOBAL_REFERENCE_WORK_DIR" \
  bash GIA/01_prepare_1000g_reference.sh
```

The historical run started with 2,504 Phase 3 samples and retained 768,584 LD-pruned variants.

### 2. Select the reference panel

The manuscript reference panel includes samples with a maximum unsupervised K=5 ancestry proportion greater than 0.80. The threshold selected 2,158 individuals.

```bash
Rscript GIA/02_select_reference_samples.R \
  "$GLOBAL_REFERENCE_WORK_DIR/reference_ld_pruned.5.Q" \
  "$GLOBAL_REFERENCE_WORK_DIR/reference_ld_pruned.fam" \
  "$REFERENCE_PSAM" \
  "$GLOBAL_REFERENCE_WORK_DIR/reference_selected" \
  0.80
```

The script writes a PLINK keep file and a two-column IID/superpopulation label file. It requires all five superpopulations and orders references as AMR, AFR, EUR, SAS, and EAS, matching the historical joint analysis.

### 3. Harmonize reference and study genotypes

The historical joint analysis used 5,378 candidate ancestry-informative SNPs. Allele-compatible data were available for 5,375 markers in 2,158 reference samples and 378 study samples.

```bash
REFERENCE_KEEP="$GLOBAL_REFERENCE_WORK_DIR/reference_selected.keep" \
WORK_DIR="$GLOBAL_JOINT_WORK_DIR" \
  bash GIA/03_prepare_joint_dataset.sh

Rscript GIA/04_build_supervised_pop.R \
  "$GLOBAL_JOINT_WORK_DIR/joint.fam" \
  "$GLOBAL_REFERENCE_WORK_DIR/reference_selected.labels.tsv" \
  "$GLOBAL_JOINT_WORK_DIR/joint.pop"
```

The harmonization script assigns `chromosome:position:reference:alternate` variant IDs to both datasets, retains shared markers, and merges the selected references before the study samples.

### 4. Estimate ancestry in study samples

```bash
JOINT_BFILE="$GLOBAL_JOINT_WORK_DIR/joint" \
  bash GIA/05_run_supervised_admixture.sh
```

ADMIXTURE reads the reference labels from `joint.pop` and estimates the five ancestry proportions for study rows marked with `-`. The historical run did not retain ADMIXTURE's terminal output. The saved `.pop`, FAM, and Q files support this supervised K=5 reconstruction.

Map the generated files to the manuscript inputs as follows:

| Generated file | Manuscript input name |
| --- | --- |
| `reference_ld_pruned.5.Q` | `gwas_ld_pruned.5.Q` |
| `reference_selected.labels.tsv` | `sample_pure.txt` |
| `joint.5.Q` | `1000G_378.5.Q` |
| `joint.fam` | `1000G_378.fam` |
| Phase 3 PSAM metadata | `all_phase3.psam` |

## Privacy and review

The scripts accept paths to controlled inputs but never include study IDs. Before a push, inspect `git status`, review staged files, and confirm that no VCF/PLINK files, sample maps, Q matrices, logs, or results are tracked.
