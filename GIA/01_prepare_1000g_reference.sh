#!/usr/bin/env bash
set -euo pipefail

reference_bfile=$1
autosomal_prefix=$2
prune_prefix=$3
ld_pruned_prefix=$4
reference_maf=$5
ld_window=$6
ld_step=$7
ld_r2=$8
k=$9

plink \
  --bfile "$reference_bfile" \
  --allow-extra-chr \
  --autosome \
  --biallelic-only strict \
  --snps-only just-acgt \
  --maf "$reference_maf" \
  --make-bed \
  --out "$autosomal_prefix"

plink \
  --bfile "$autosomal_prefix" \
  --indep-pairwise "$ld_window" "$ld_step" "$ld_r2" \
  --out "$prune_prefix"

plink \
  --bfile "$autosomal_prefix" \
  --extract "$prune_prefix.prune.in" \
  --make-bed \
  --out "$ld_pruned_prefix"

(
  cd "$(dirname "$ld_pruned_prefix")"
  admixture "$(basename "$ld_pruned_prefix").bed" "$k"
)
