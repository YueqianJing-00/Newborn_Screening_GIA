#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../lib/common.sh"

require_env REFERENCE_VCF
require_env REFERENCE_LABELS
require_env WORK_DIR

BCFTOOLS=${BCFTOOLS:-bcftools}
TABIX=${TABIX:-tabix}
REFERENCE_MIN_AF=${REFERENCE_MIN_AF:-0.01}
BCFTOOLS_THREADS=${BCFTOOLS_THREADS:-20}

require_command "$BCFTOOLS"
require_command "$TABIX"
require_file "$REFERENCE_VCF"
require_file "$REFERENCE_LABELS"
prepare_output_directory "$WORK_DIR"

sample_ids="$WORK_DIR/reference_sample_ids.txt"
sample_map="$WORK_DIR/reference_sample_map.txt"
common_vcf="$WORK_DIR/reference_common.vcf.gz"
biallelic_vcf="$WORK_DIR/reference_common_biallelic_snps.vcf.gz"
selected_vcf="$WORK_DIR/rfmix_reference.vcf.gz"

printf 'RUN validate-reference-labels %s\n' "$REFERENCE_LABELS"
printf 'RUN write-reference-ids %s > %s\n' "$REFERENCE_LABELS" "$sample_ids"
printf 'RUN write-reference-map %s > %s\n' "$REFERENCE_LABELS" "$sample_map"
if [[ $DRY_RUN != 1 ]]; then
  awk '
    BEGIN {
      valid["AFR"] = valid["AMR"] = valid["EAS"] = 1
      valid["EUR"] = valid["SAS"] = 1
    }
    NF != 2 { print "reference labels must have two columns" > "/dev/stderr"; exit 1 }
    !($2 in valid) { print "unexpected superpopulation: " $2 > "/dev/stderr"; exit 1 }
    seen[$1]++ { print "duplicate reference IID: " $1 > "/dev/stderr"; exit 1 }
    { groups[$2]++ }
    END {
      for (group in valid) {
        if (!groups[group]) {
          print "missing superpopulation: " group > "/dev/stderr"
          exit 1
        }
      }
    }
  ' "$REFERENCE_LABELS"
  awk '{print $1}' "$REFERENCE_LABELS" >"$sample_ids"
  awk 'BEGIN {OFS="\t"} {print $1, $2}' "$REFERENCE_LABELS" >"$sample_map"
fi

run_command "$BCFTOOLS" view \
  -q "$REFERENCE_MIN_AF" \
  -Oz \
  -o "$common_vcf" \
  --threads "$BCFTOOLS_THREADS" \
  "$REFERENCE_VCF"

run_command "$BCFTOOLS" view \
  -m2 -M2 -v snps \
  -Oz \
  -o "$biallelic_vcf" \
  --threads "$BCFTOOLS_THREADS" \
  "$common_vcf"

run_command "$BCFTOOLS" view \
  -S "$sample_ids" \
  -Oz \
  -o "$selected_vcf" \
  --threads "$BCFTOOLS_THREADS" \
  "$biallelic_vcf"

run_command "$TABIX" -f -p vcf "$selected_vcf"
