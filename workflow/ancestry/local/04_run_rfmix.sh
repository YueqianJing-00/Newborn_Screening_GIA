#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../lib/common.sh"

require_env PHASED_QUERY_PATTERN
require_env REFERENCE_VCF
require_env REFERENCE_SAMPLE_MAP
require_env GENETIC_MAP
require_env CHROMOSOME
require_env WORK_DIR

RFMIX=${RFMIX:-rfmix}
RFMIX_THREADS=${RFMIX_THREADS:-8}

[[ $CHROMOSOME =~ ^([1-9]|1[0-9]|2[0-2]|X)$ ]] || {
  die "CHROMOSOME must be 1-22 or X"
}
[[ $PHASED_QUERY_PATTERN == *'%s'* ]] || {
  die "PHASED_QUERY_PATTERN must contain %s for the chromosome"
}

phased_query=${PHASED_QUERY_PATTERN//%s/$CHROMOSOME}
output_prefix="$WORK_DIR/chr${CHROMOSOME}_local_ancestry"

require_command "$RFMIX"
require_file "$phased_query"
require_file "$REFERENCE_VCF"
require_file "$REFERENCE_SAMPLE_MAP"
require_file "$GENETIC_MAP"
prepare_output_directory "$WORK_DIR"

run_command "$RFMIX" \
  -f "$phased_query" \
  -r "$REFERENCE_VCF" \
  -m "$REFERENCE_SAMPLE_MAP" \
  -g "$GENETIC_MAP" \
  -o "$output_prefix" \
  --chromosome="$CHROMOSOME" \
  --n-threads="$RFMIX_THREADS"
