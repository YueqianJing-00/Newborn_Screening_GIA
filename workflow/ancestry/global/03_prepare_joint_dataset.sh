#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../lib/common.sh"

require_env REFERENCE_BFILE
require_env AIM_LIST
require_env REFERENCE_KEEP
require_env WORK_DIR

if [[ -n ${STUDY_VCF:-} && -n ${STUDY_BFILE:-} ]]; then
  die "set only one of STUDY_VCF or STUDY_BFILE"
fi
if [[ -z ${STUDY_VCF:-} && -z ${STUDY_BFILE:-} ]]; then
  die "set STUDY_VCF or STUDY_BFILE"
fi

PLINK=${PLINK:-plink}
PLINK2=${PLINK2:-plink2}
VARIANT_ID_TEMPLATE='@:#:$r:$a'

require_command "$PLINK"
require_command "$PLINK2"
require_bfile "$REFERENCE_BFILE"
require_file "$AIM_LIST"
require_file "$REFERENCE_KEEP"
if [[ -n ${STUDY_VCF:-} ]]; then
  require_file "$STUDY_VCF"
else
  require_bfile "$STUDY_BFILE"
fi
prepare_output_directory "$WORK_DIR"

reference_aims="$WORK_DIR/reference_aims"
reference_canonical="$WORK_DIR/reference_canonical"
study_canonical="$WORK_DIR/study_canonical"
reference_shared="$WORK_DIR/reference_shared"
study_shared="$WORK_DIR/study_shared"
joint_prefix="$WORK_DIR/joint"
reference_ids="$WORK_DIR/reference_variant_ids.txt"
study_ids="$WORK_DIR/study_variant_ids.txt"
shared_ids="$WORK_DIR/shared_variant_ids.txt"

run_command "$PLINK" \
  --bfile "$REFERENCE_BFILE" \
  --allow-extra-chr \
  --extract "$AIM_LIST" \
  --keep "$REFERENCE_KEEP" \
  --indiv-sort f "$REFERENCE_KEEP" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_aims"

run_command "$PLINK2" \
  --bfile "$reference_aims" \
  --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
  --make-bed \
  --out "$reference_canonical"

if [[ -n ${STUDY_VCF:-} ]]; then
  run_command "$PLINK2" \
    --vcf "$STUDY_VCF" \
    --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
    --make-bed \
    --out "$study_canonical"
else
  run_command "$PLINK2" \
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
  [[ -s $shared_ids ]] || die "reference and study datasets have no shared variants"
fi

run_command "$PLINK" \
  --bfile "$reference_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_shared"

run_command "$PLINK" \
  --bfile "$study_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$study_shared"

run_command "$PLINK" \
  --bfile "$reference_shared" \
  --bmerge "$study_shared" \
  --indiv-sort 0 \
  --keep-allele-order \
  --make-bed \
  --out "$joint_prefix"
