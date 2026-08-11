# Ancestry estimation workflows

These scripts generate the global-ancestry files used by the manuscript analyses. They replace exploratory notebooks and cluster jobs that contained private paths or participant identifiers. Keep controlled genotype data and all generated files outside the repository.

## Software

The workflow uses [PLINK 1.9 stable build 1.9.0-b.7.11 (19 August 2025)](https://www.cog-genomics.org/plink/1.9/) for every genotype operation and ADMIXTURE 1.3.0 for ancestry estimation. Install these programs and R before running the workflow, and make `plink`, `admixture`, and `Rscript` available on `PATH`. PLINK 2 is not required.

Each script contains an editable path and parameter block near the top and accepts no command-line arguments. The scripts do not perform software or file preflight checks or create output directories.

Edit the paths in each script, then create the configured output directories. For the example paths currently shown in the scripts:

```bash
mkdir -p /secure/work/global_reference /secure/work/global_joint
```

## Global ancestry

### 1. Estimate ancestry in 1000 Genomes with MSK-IMPACT

The reference workflow extracts the 5,378-SNP MSK-IMPACT panel from the 1000 Genomes dataset and runs unsupervised ADMIXTURE at K=5. The resulting ancestry matrix is used to select homogeneous reference samples.

```bash
bash GIA/01_prepare_1000g_reference.sh
```

The panel is used consistently for both reference selection and final study GIA estimation. Run all subsequent steps again whenever the reference selection changes so the keep file, joint dataset, and supervised ADMIXTURE output remain synchronized.

### 2. Select homogeneous reference samples

Reference samples are selected when their maximum unsupervised K=5 ancestry proportion exceeds the configured threshold. The script records the resulting count rather than assuming a fixed number of references.

```bash
Rscript GIA/02_select_reference_samples.R
```

The script writes a PLINK keep file and a two-column IID/superpopulation label file. It requires all five superpopulations and orders references as AMR, AFR, EUR, SAS, and EAS, matching the historical joint analysis.

### 3. Build the MSK-IMPACT joint dataset

The workflow reuses `reference_impact` from Step 1, then keeps and orders the references selected in Step 2. It does not extract the MSK-IMPACT SNPs from the original 1000 Genomes dataset again. PLINK then retains the allele-compatible panel markers shared by the selected references and study genotypes.
Set `study_input_option` in `03_prepare_joint_dataset.sh` to `--vcf` for a VCF or `--bfile` for a PLINK binary-file prefix.

```bash
bash GIA/03_prepare_joint_dataset.sh

Rscript GIA/04_build_supervised_pop.R
```

The harmonization script preserves allele order, assigns `chromosome:position:A2:A1` variant IDs to both datasets by updating their `.bim` files, retains shared markers, and merges the selected references before the study samples. For VCF input, PLINK stores REF as A2 and ALT as A1.

### 4. Estimate study GIA with MSK-IMPACT

```bash
bash GIA/05_run_supervised_admixture.sh
```

ADMIXTURE reads the reference labels from `joint.pop` and estimates five ancestry proportions from the shared MSK-IMPACT markers for study rows marked with `-`. The historical run did not retain ADMIXTURE's terminal output. The saved `.pop`, FAM, and Q files support this supervised K=5 reconstruction.

Map the generated files to the manuscript inputs as follows:

| Generated file | Role | Manuscript input name |
| --- | --- | --- |
| `reference_impact.5.Q` | MSK-IMPACT reference estimation and selection | `1000G_impact.5.Q` |
| `reference_selected.labels.tsv` | Selected reference labels | `sample_pure.txt` |
| `joint.5.Q` | Final MSK-IMPACT GIA estimates | `1000G_378.5.Q` |
| `joint.fam` | Final reference and study sample order | `1000G_378.fam` |
| Phase 3 PSAM metadata | Reference superpopulation metadata | `all_phase3.psam` |

## Privacy and review

The committed scripts contain placeholder paths and no study IDs. If you edit them with private paths, restore the placeholders before committing. Before a push, inspect `git status`, review staged files, and confirm that no VCF/PLINK files, sample maps, Q matrices, logs, or results are tracked.
