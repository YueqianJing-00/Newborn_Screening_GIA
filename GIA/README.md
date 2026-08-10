# Ancestry estimation workflows

These scripts generate the global-ancestry files used by the manuscript analyses. They replace exploratory notebooks and cluster jobs that contained private paths or participant identifiers. Keep controlled genotype data and all generated files outside the repository.

## Software

The historical analyses used PLINK 1.90b6.21, PLINK 2.00a3.7LM, and ADMIXTURE 1.3.0. Install these programs and R before running the workflow, and make `plink`, `plink2`, `admixture`, and `Rscript` available on `PATH`. The scripts receive all paths and parameters as arguments; they do not perform software or file preflight checks or create output directories.

Copy the example configuration and edit the paths:

```bash
cp GIA/config.example.env GIA/config.env
set -a
source GIA/config.env
set +a
mkdir -p "$REFERENCE_OUTPUT_DIR" "$JOINT_OUTPUT_DIR"
```

## Global ancestry

### 1. Build the reference-selection matrix

The reference-selection workflow retains autosomal, strict biallelic A/C/G/T SNPs with MAF at least 1%, applies LD pruning with a 1,000-variant window, 100-variant step, and r2 threshold of 0.2, then runs unsupervised ADMIXTURE at K=5. This LD-pruned matrix is used only to select homogeneous 1000 Genomes reference samples; it is not the marker set used for final study GIA estimation.

```bash
bash GIA/01_prepare_1000g_reference.sh \
  "$REFERENCE_BFILE" \
  "$AUTOSOMAL_PREFIX" \
  "$PRUNE_PREFIX" \
  "$LD_PRUNED_PREFIX" \
  "$REFERENCE_MAF" \
  "$LD_WINDOW" \
  "$LD_STEP" \
  "$LD_R2" \
  "$K"
```

The historical run started with 2,504 Phase 3 samples and retained 768,584 LD-pruned variants.

### 2. Select homogeneous reference samples

The manuscript reference panel includes samples with a maximum unsupervised K=5 ancestry proportion greater than 0.80. The threshold selected 2,158 individuals.

```bash
Rscript GIA/02_select_reference_samples.R \
  "$REFERENCE_Q" \
  "$REFERENCE_FAM" \
  "$REFERENCE_PSAM" \
  "$REFERENCE_KEEP" \
  "$REFERENCE_LABELS" \
  "$REFERENCE_THRESHOLD"
```

The script writes a PLINK keep file and a two-column IID/superpopulation label file. It requires all five superpopulations and orders references as AMR, AFR, EUR, SAS, and EAS, matching the historical joint analysis.

### 3. Build the MSK-IMPACT joint dataset

Final study GIA estimation uses the 5,378-SNP MSK-IMPACT panel rather than the LD-pruned reference-selection marker set. Allele-compatible data were available for 5,375 of these markers in 2,158 reference samples and 378 study samples.
Set `STUDY_INPUT_OPTION` to `--vcf` for a VCF or `--bfile` for a PLINK binary-file prefix.

```bash
bash GIA/03_prepare_joint_dataset.sh \
  "$REFERENCE_BFILE" \
  "$MSK_IMPACT_SNP_LIST" \
  "$REFERENCE_KEEP" \
  "$STUDY_INPUT_OPTION" \
  "$STUDY_INPUT" \
  "$VARIANT_ID_TEMPLATE" \
  "$SELECTED_REFERENCE_IMPACT_PREFIX" \
  "$REFERENCE_CANONICAL_PREFIX" \
  "$STUDY_CANONICAL_PREFIX" \
  "$REFERENCE_IDS" \
  "$STUDY_IDS" \
  "$SHARED_IDS" \
  "$REFERENCE_SHARED_PREFIX" \
  "$STUDY_SHARED_PREFIX" \
  "$JOINT_PREFIX"

Rscript GIA/04_build_supervised_pop.R \
  "$JOINT_FAM" \
  "$REFERENCE_LABELS" \
  "$JOINT_POP"
```

The harmonization script assigns `chromosome:position:reference:alternate` variant IDs to both datasets, retains shared markers, and merges the selected references before the study samples.

### 4. Estimate study GIA with MSK-IMPACT

```bash
bash GIA/05_run_supervised_admixture.sh \
  "$JOINT_PREFIX" \
  "$K"
```

ADMIXTURE reads the reference labels from `joint.pop` and estimates five ancestry proportions from the shared MSK-IMPACT markers for study rows marked with `-`. The historical run did not retain ADMIXTURE's terminal output. The saved `.pop`, FAM, and Q files support this supervised K=5 reconstruction.

Map the generated files to the manuscript inputs as follows:

| Generated file | Role | Manuscript input name |
| --- | --- | --- |
| `reference_ld_pruned.5.Q` | Reference selection only | `gwas_ld_pruned.5.Q` |
| `reference_selected.labels.tsv` | Selected reference labels | `sample_pure.txt` |
| `joint.5.Q` | Final MSK-IMPACT GIA estimates | `1000G_378.5.Q` |
| `joint.fam` | Final reference and study sample order | `1000G_378.fam` |
| Phase 3 PSAM metadata | Reference superpopulation metadata | `all_phase3.psam` |

## Privacy and review

The scripts accept paths to controlled inputs but never include study IDs. Before a push, inspect `git status`, review staged files, and confirm that no VCF/PLINK files, sample maps, Q matrices, logs, or results are tracked.
