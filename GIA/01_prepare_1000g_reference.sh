#!/usr/bin/env bash
set -euo pipefail

PLINK=${PLINK:-plink}
ADMIXTURE=${ADMIXTURE:-admixture}
K=${K:-5}
REFERENCE_MAF=${REFERENCE_MAF:-0.01}
LD_WINDOW=${LD_WINDOW:-1000}
LD_STEP=${LD_STEP:-100}
LD_R2=${LD_R2:-0.2}

mkdir -p "$WORK_DIR"

autosomal_prefix="$WORK_DIR/reference_autosomal"
prune_prefix="$WORK_DIR/reference_prune"
ld_pruned_prefix="$WORK_DIR/reference_ld_pruned"

"$PLINK" \
  --bfile "$REFERENCE_BFILE" \
  --allow-extra-chr \
  --autosome \
  --biallelic-only strict \
  --snps-only just-acgt \
  --maf "$REFERENCE_MAF" \
  --make-bed \
  --out "$autosomal_prefix"

"$PLINK" \
  --bfile "$autosomal_prefix" \
  --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
  --out "$prune_prefix"

"$PLINK" \
  --bfile "$autosomal_prefix" \
  --extract "$prune_prefix.prune.in" \
  --make-bed \
  --out "$ld_pruned_prefix"

(
  cd "$WORK_DIR"
  "$ADMIXTURE" "$(basename "$ld_pruned_prefix").bed" "$K"
)
