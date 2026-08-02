#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../lib/common.sh"

require_env QUERY_BFILE
require_env WORK_DIR

PLINK=${PLINK:-plink}
TABIX=${TABIX:-tabix}
QUERY_MAX_MISSING=${QUERY_MAX_MISSING:-0.05}
QUERY_MAF=${QUERY_MAF:-0.01}

require_command "$PLINK"
require_command "$TABIX"
require_bfile "$QUERY_BFILE"
prepare_output_directory "$WORK_DIR"

qc_prefix="$WORK_DIR/query_qc"
common_prefix="$WORK_DIR/query_common"

run_command "$PLINK" \
  --bfile "$QUERY_BFILE" \
  --geno "$QUERY_MAX_MISSING" \
  --keep-allele-order \
  --make-bed \
  --out "$qc_prefix"

run_command "$PLINK" \
  --bfile "$qc_prefix" \
  --maf "$QUERY_MAF" \
  --keep-allele-order \
  --make-bed \
  --out "$common_prefix"

run_command "$PLINK" \
  --bfile "$common_prefix" \
  --snps-only \
  --keep-allele-order \
  --recode vcf bgz \
  --out "$common_prefix"

run_command "$TABIX" -f -p vcf "$common_prefix.vcf.gz"
