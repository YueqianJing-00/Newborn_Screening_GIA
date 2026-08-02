#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../lib/common.sh"

require_env JOINT_BFILE

ADMIXTURE=${ADMIXTURE:-admixture}
K=${K:-5}

[[ $K =~ ^[1-9][0-9]*$ ]] || die "K must be a positive integer"

require_command "$ADMIXTURE"
require_bfile "$JOINT_BFILE"
require_file "$JOINT_BFILE.pop"

if [[ $DRY_RUN != 1 ]]; then
  fam_rows=$(wc -l <"$JOINT_BFILE.fam")
  pop_rows=$(wc -l <"$JOINT_BFILE.pop")
  [[ $fam_rows -eq $pop_rows ]] || {
    die "joint FAM and POP row counts differ ($fam_rows vs $pop_rows)"
  }
fi

run_in_directory "$(dirname "$JOINT_BFILE")" \
  "$ADMIXTURE" --supervised "$(basename "$JOINT_BFILE").bed" "$K"
