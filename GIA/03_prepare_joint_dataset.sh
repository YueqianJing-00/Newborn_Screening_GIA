#!/usr/bin/env bash
set -euo pipefail

# Edit these paths and parameters before running.
reference_output_dir="/secure/work/global_reference"
joint_output_dir="/secure/work/global_joint"
reference_impact_prefix="${reference_output_dir}/reference_impact"
reference_keep="${reference_output_dir}/reference_selected.keep"
study_input_option=--vcf
study_input="/secure/path/to/study_genotypes.vcf.gz"
selected_reference_impact="${joint_output_dir}/selected_reference_impact"
reference_canonical="${joint_output_dir}/reference_canonical"
study_canonical="${joint_output_dir}/study_canonical"
reference_ids="${joint_output_dir}/reference_variant_ids.txt"
study_ids="${joint_output_dir}/study_variant_ids.txt"
shared_ids="${joint_output_dir}/shared_variant_ids.txt"
reference_shared="${joint_output_dir}/reference_shared"
study_shared="${joint_output_dir}/study_shared"
joint_prefix="${joint_output_dir}/joint"

# Rename variants as chromosome:position:A2:A1; VCF imports place REF in A2.
canonicalize_bim_ids() {
  local bim_path="$1"
  awk 'BEGIN { OFS = "\t" } { $2 = $1 ":" $4 ":" $6 ":" $5; print }' \
    "$bim_path" >"${bim_path}.tmp"
  mv "${bim_path}.tmp" "$bim_path"
}

# Keep and order the selected references from the panel prepared in Step 1.
plink \
  --bfile "$reference_impact_prefix" \
  --keep "$reference_keep" \
  --indiv-sort f "$reference_keep" \
  --keep-allele-order \
  --make-bed \
  --out "$selected_reference_impact"

# Give reference variants coordinate-and-allele IDs for dataset matching.
plink \
  --bfile "$selected_reference_impact" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_canonical"
canonicalize_bim_ids "$reference_canonical.bim"

# Convert study genotypes to PLINK with the same variant-ID convention.
plink \
  "$study_input_option" "$study_input" \
  --keep-allele-order \
  --make-bed \
  --out "$study_canonical"
canonicalize_bim_ids "$study_canonical.bim"

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
