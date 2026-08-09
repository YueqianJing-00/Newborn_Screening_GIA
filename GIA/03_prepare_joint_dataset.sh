#!/usr/bin/env bash
set -euo pipefail

PLINK=${PLINK:-plink}
PLINK2=${PLINK2:-plink2}
VARIANT_ID_TEMPLATE='@:#:$r:$a'

mkdir -p "$WORK_DIR"

reference_aims="$WORK_DIR/reference_aims"
reference_canonical="$WORK_DIR/reference_canonical"
study_canonical="$WORK_DIR/study_canonical"
reference_shared="$WORK_DIR/reference_shared"
study_shared="$WORK_DIR/study_shared"
joint_prefix="$WORK_DIR/joint"
reference_ids="$WORK_DIR/reference_variant_ids.txt"
study_ids="$WORK_DIR/study_variant_ids.txt"
shared_ids="$WORK_DIR/shared_variant_ids.txt"

"$PLINK" \
  --bfile "$REFERENCE_BFILE" \
  --allow-extra-chr \
  --extract "$AIM_LIST" \
  --keep "$REFERENCE_KEEP" \
  --indiv-sort f "$REFERENCE_KEEP" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_aims"

"$PLINK2" \
  --bfile "$reference_aims" \
  --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
  --make-bed \
  --out "$reference_canonical"

if [[ -n ${STUDY_VCF:-} ]]; then
  "$PLINK2" \
    --vcf "$STUDY_VCF" \
    --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
    --make-bed \
    --out "$study_canonical"
else
  "$PLINK2" \
    --bfile "$STUDY_BFILE" \
    --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
    --make-bed \
    --out "$study_canonical"
fi

cut -f2 "$reference_canonical.bim" | LC_ALL=C sort -u >"$reference_ids"
cut -f2 "$study_canonical.bim" | LC_ALL=C sort -u >"$study_ids"
comm -12 "$reference_ids" "$study_ids" >"$shared_ids"

"$PLINK" \
  --bfile "$reference_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_shared"

"$PLINK" \
  --bfile "$study_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$study_shared"

"$PLINK" \
  --bfile "$reference_shared" \
  --bmerge "$study_shared" \
  --indiv-sort 0 \
  --keep-allele-order \
  --make-bed \
  --out "$joint_prefix"
