#!/usr/bin/env bash
set -euo pipefail

# Edit these paths and parameters before running.
reference_bfile="/path/to/1000g_phase3/all_phase3"
msk_impact_snp_list="/path/to/msk_impact_5378_snps.txt"
reference_output_dir="/secure/work/global_reference"
joint_output_dir="/secure/work/global_joint"
reference_keep="${reference_output_dir}/reference_selected.keep"
study_input_option=--vcf
study_input="/secure/path/to/study_genotypes.vcf.gz"
variant_id_template='@:#:$r:$a'
selected_reference_impact="${joint_output_dir}/selected_reference_impact"
reference_canonical="${joint_output_dir}/reference_canonical"
study_canonical="${joint_output_dir}/study_canonical"
reference_ids="${joint_output_dir}/reference_variant_ids.txt"
study_ids="${joint_output_dir}/study_variant_ids.txt"
shared_ids="${joint_output_dir}/shared_variant_ids.txt"
reference_shared="${joint_output_dir}/reference_shared"
study_shared="${joint_output_dir}/study_shared"
joint_prefix="${joint_output_dir}/joint"

# Extract the MSK-IMPACT SNPs and selected 1000 Genomes references.
plink \
  --bfile "$reference_bfile" \
  --allow-extra-chr \
  --extract "$msk_impact_snp_list" \
  --keep "$reference_keep" \
  --indiv-sort f "$reference_keep" \
  --keep-allele-order \
  --make-bed \
  --out "$selected_reference_impact"

# Give reference variants coordinate-and-allele IDs for dataset matching.
plink2 \
  --bfile "$selected_reference_impact" \
  --set-all-var-ids "$variant_id_template" \
  --make-bed \
  --out "$reference_canonical"

# Convert study genotypes to PLINK with the same variant-ID convention.
plink2 \
  "$study_input_option" "$study_input" \
  --set-all-var-ids "$variant_id_template" \
  --make-bed \
  --out "$study_canonical"

# Identify variant IDs shared by the reference and study datasets.
cut -f2 "$reference_canonical.bim" | LC_ALL=C sort -u >"$reference_ids"
cut -f2 "$study_canonical.bim" | LC_ALL=C sort -u >"$study_ids"
comm -12 "$reference_ids" "$study_ids" >"$shared_ids"

# Restrict both datasets to the same shared variants and allele order.
plink \
  --bfile "$reference_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_shared"

plink \
  --bfile "$study_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$study_shared"

# Merge references first and study samples second for supervised ADMIXTURE.
plink \
  --bfile "$reference_shared" \
  --bmerge "$study_shared" \
  --indiv-sort 0 \
  --keep-allele-order \
  --make-bed \
  --out "$joint_prefix"
