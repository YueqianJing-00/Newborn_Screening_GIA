#!/usr/bin/env bash
set -euo pipefail

# Edit these paths and parameters before running.
reference_bfile="/path/to/1000g_phase3/all_phase3"
msk_impact_snp_list="/path/to/msk_impact_5378_snps.txt"
reference_output_dir="/secure/work/global_reference"
reference_impact_prefix="${reference_output_dir}/reference_impact"
k=5

# Extract the 5,378 MSK-IMPACT SNPs into a PLINK reference dataset.
plink \
  --bfile "$reference_bfile" \
  --allow-extra-chr \
  --extract "$msk_impact_snp_list" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_impact_prefix"

# Estimate K ancestry components in the 1000 Genomes reference samples.
(
  cd "$(dirname "$reference_impact_prefix")"
  admixture "$(basename "$reference_impact_prefix").bed" "$k"
)
