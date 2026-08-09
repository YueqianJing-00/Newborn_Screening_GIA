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

for required_variable in REFERENCE_BFILE AIM_LIST REFERENCE_KEEP WORK_DIR; do
  [[ -n ${!required_variable:-} ]] || fail "$required_variable is required"
done

if [[ -n ${STUDY_VCF:-} && -n ${STUDY_BFILE:-} ]]; then
  fail "set only one of STUDY_VCF or STUDY_BFILE"
fi
if [[ -z ${STUDY_VCF:-} && -z ${STUDY_BFILE:-} ]]; then
  fail "set STUDY_VCF or STUDY_BFILE"
fi

PLINK=${PLINK:-plink}
PLINK2=${PLINK2:-plink2}
VARIANT_ID_TEMPLATE='@:#:$r:$a'

if [[ $DRY_RUN != 1 ]]; then
  for required_command in "$PLINK" "$PLINK2"; do
    command -v "$required_command" >/dev/null 2>&1 || {
      fail "required command not found: $required_command"
    }
  done
  for extension in bed bim fam; do
    [[ -f ${REFERENCE_BFILE}.${extension} ]] || {
      fail "missing PLINK file: ${REFERENCE_BFILE}.${extension}"
    }
  done
  [[ -f $AIM_LIST ]] || fail "missing file: $AIM_LIST"
  [[ -f $REFERENCE_KEEP ]] || fail "missing file: $REFERENCE_KEEP"
  if [[ -n ${STUDY_VCF:-} ]]; then
    [[ -f $STUDY_VCF ]] || fail "missing file: $STUDY_VCF"
  else
    for extension in bed bim fam; do
      [[ -f ${STUDY_BFILE}.${extension} ]] || {
        fail "missing PLINK file: ${STUDY_BFILE}.${extension}"
      }
    done
  fi
  mkdir -p "$WORK_DIR"
fi

reference_aims="$WORK_DIR/reference_aims"
reference_canonical="$WORK_DIR/reference_canonical"
study_canonical="$WORK_DIR/study_canonical"
reference_shared="$WORK_DIR/reference_shared"
study_shared="$WORK_DIR/study_shared"
joint_prefix="$WORK_DIR/joint"
reference_ids="$WORK_DIR/reference_variant_ids.txt"
study_ids="$WORK_DIR/study_variant_ids.txt"
shared_ids="$WORK_DIR/shared_variant_ids.txt"

run "$PLINK" \
  --bfile "$REFERENCE_BFILE" \
  --allow-extra-chr \
  --extract "$AIM_LIST" \
  --keep "$REFERENCE_KEEP" \
  --indiv-sort f "$REFERENCE_KEEP" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_aims"

run "$PLINK2" \
  --bfile "$reference_aims" \
  --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
  --make-bed \
  --out "$reference_canonical"

if [[ -n ${STUDY_VCF:-} ]]; then
  run "$PLINK2" \
    --vcf "$STUDY_VCF" \
    --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
    --make-bed \
    --out "$study_canonical"
else
  run "$PLINK2" \
    --bfile "$STUDY_BFILE" \
    --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
    --make-bed \
    --out "$study_canonical"
fi

printf 'RUN cut -f2 %s.bim | sort -u > %s\n' "$reference_canonical" "$reference_ids"
printf 'RUN cut -f2 %s.bim | sort -u > %s\n' "$study_canonical" "$study_ids"
printf 'RUN comm -12 %s %s > %s\n' "$reference_ids" "$study_ids" "$shared_ids"
if [[ $DRY_RUN != 1 ]]; then
  cut -f2 "$reference_canonical.bim" | LC_ALL=C sort -u >"$reference_ids"
  cut -f2 "$study_canonical.bim" | LC_ALL=C sort -u >"$study_ids"
  comm -12 "$reference_ids" "$study_ids" >"$shared_ids"
  [[ -s $shared_ids ]] || fail "reference and study datasets have no shared variants"
fi

run "$PLINK" \
  --bfile "$reference_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_shared"

run "$PLINK" \
  --bfile "$study_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$study_shared"

run "$PLINK" \
  --bfile "$reference_shared" \
  --bmerge "$study_shared" \
  --indiv-sort 0 \
  --keep-allele-order \
  --make-bed \
  --out "$joint_prefix"
