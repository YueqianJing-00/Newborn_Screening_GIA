#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../lib/common.sh"

require_env REFERENCE_BFILE
require_env WORK_DIR

PLINK=${PLINK:-plink}
ADMIXTURE=${ADMIXTURE:-admixture}
K=${K:-5}
REFERENCE_MAF=${REFERENCE_MAF:-0.01}
LD_WINDOW=${LD_WINDOW:-1000}
LD_STEP=${LD_STEP:-100}
LD_R2=${LD_R2:-0.2}

[[ $K =~ ^[1-9][0-9]*$ ]] || die "K must be a positive integer"

require_command "$PLINK"
require_command "$ADMIXTURE"
require_bfile "$REFERENCE_BFILE"
prepare_output_directory "$WORK_DIR"

autosomal_prefix="$WORK_DIR/reference_autosomal"
prune_prefix="$WORK_DIR/reference_prune"
ld_pruned_prefix="$WORK_DIR/reference_ld_pruned"

run_command "$PLINK" \
  --bfile "$REFERENCE_BFILE" \
  --allow-extra-chr \
  --autosome \
  --biallelic-only strict \
  --snps-only just-acgt \
  --maf "$REFERENCE_MAF" \
  --make-bed \
  --out "$autosomal_prefix"

run_command "$PLINK" \
  --bfile "$autosomal_prefix" \
  --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
  --out "$prune_prefix"

run_command "$PLINK" \
  --bfile "$autosomal_prefix" \
  --extract "$prune_prefix.prune.in" \
  --make-bed \
  --out "$ld_pruned_prefix"

run_in_directory "$WORK_DIR" \
  "$ADMIXTURE" "$(basename "$ld_pruned_prefix").bed" "$K"
