#!/usr/bin/env bash
set -euo pipefail

reference_bfile=$1
msk_impact_snp_list=$2
reference_impact_prefix=$3
k=$4

plink \
  --bfile "$reference_bfile" \
  --allow-extra-chr \
  --extract "$msk_impact_snp_list" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_impact_prefix"

(
  cd "$(dirname "$reference_impact_prefix")"
  admixture "$(basename "$reference_impact_prefix").bed" "$k"
)
