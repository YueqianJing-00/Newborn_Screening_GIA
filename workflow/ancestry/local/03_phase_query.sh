#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../lib/common.sh"

require_env QUERY_VCF
require_env CHROMOSOME
require_env GENETIC_MAP_PATTERN
require_env WORK_DIR

SHAPEIT4=${SHAPEIT4:-shapeit4}
SHAPEIT4_THREADS=${SHAPEIT4_THREADS:-8}

[[ $CHROMOSOME =~ ^([1-9]|1[0-9]|2[0-2]|X)$ ]] || {
  die "CHROMOSOME must be 1-22 or X"
}
[[ $GENETIC_MAP_PATTERN == *'%s'* ]] || {
  die "GENETIC_MAP_PATTERN must contain %s for the chromosome"
}

genetic_map=${GENETIC_MAP_PATTERN//%s/$CHROMOSOME}
output_vcf="$WORK_DIR/phased_chr${CHROMOSOME}.vcf.gz"

require_command "$SHAPEIT4"
require_file "$QUERY_VCF"
require_file "$genetic_map"
prepare_output_directory "$WORK_DIR"

run_command "$SHAPEIT4" \
  --input "$QUERY_VCF" \
  --map "$genetic_map" \
  --region "$CHROMOSOME" \
  --output "$output_vcf" \
  --thread "$SHAPEIT4_THREADS" \
  --sequencing
