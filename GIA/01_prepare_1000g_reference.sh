#!/usr/bin/env bash
set -euo pipefail

DRY_RUN=${DRY_RUN:-0}

fail() {
  printf 'ERROR: %s\n' "$1" >&2
  exit 1
}

run() {
  printf 'RUN'
  printf ' %s' "$@"
  printf '\n'
  [[ $DRY_RUN == 1 ]] || "$@"
}

for required_variable in REFERENCE_BFILE WORK_DIR; do
  [[ -n ${!required_variable:-} ]] || fail "$required_variable is required"
done

PLINK=${PLINK:-plink}
ADMIXTURE=${ADMIXTURE:-admixture}
K=${K:-5}
REFERENCE_MAF=${REFERENCE_MAF:-0.01}
LD_WINDOW=${LD_WINDOW:-1000}
LD_STEP=${LD_STEP:-100}
LD_R2=${LD_R2:-0.2}

[[ $K =~ ^[1-9][0-9]*$ ]] || fail "K must be a positive integer"

if [[ $DRY_RUN != 1 ]]; then
  for required_command in "$PLINK" "$ADMIXTURE"; do
    command -v "$required_command" >/dev/null 2>&1 || fail "required command not found: $required_command"
  done
  for extension in bed bim fam; do
    [[ -f ${REFERENCE_BFILE}.${extension} ]] || fail "missing PLINK file: ${REFERENCE_BFILE}.${extension}"
  done
  mkdir -p "$WORK_DIR"
fi

autosomal_prefix="$WORK_DIR/reference_autosomal"
prune_prefix="$WORK_DIR/reference_prune"
ld_pruned_prefix="$WORK_DIR/reference_ld_pruned"

run "$PLINK" \
  --bfile "$REFERENCE_BFILE" \
  --allow-extra-chr \
  --autosome \
  --biallelic-only strict \
  --snps-only just-acgt \
  --maf "$REFERENCE_MAF" \
  --make-bed \
  --out "$autosomal_prefix"

run "$PLINK" \
  --bfile "$autosomal_prefix" \
  --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
  --out "$prune_prefix"

run "$PLINK" \
  --bfile "$autosomal_prefix" \
  --extract "$prune_prefix.prune.in" \
  --make-bed \
  --out "$ld_pruned_prefix"

admixture_input="$(basename "$ld_pruned_prefix").bed"
printf 'RUN cd %s && %s %s %s\n' \
  "$WORK_DIR" "$ADMIXTURE" "$admixture_input" "$K"
if [[ $DRY_RUN != 1 ]]; then
  (
    cd "$WORK_DIR"
    "$ADMIXTURE" "$admixture_input" "$K"
  )
fi
