#!/usr/bin/env bash
set -euo pipefail

DRY_RUN=${DRY_RUN:-0}

fail() {
  printf 'ERROR: %s\n' "$1" >&2
  exit 1
}

[[ -n ${JOINT_BFILE:-} ]] || fail "JOINT_BFILE is required"

ADMIXTURE=${ADMIXTURE:-admixture}
K=${K:-5}

[[ $K =~ ^[1-9][0-9]*$ ]] || fail "K must be a positive integer"

if [[ $DRY_RUN != 1 ]]; then
  command -v "$ADMIXTURE" >/dev/null 2>&1 || {
    fail "required command not found: $ADMIXTURE"
  }
  for extension in bed bim fam; do
    [[ -f ${JOINT_BFILE}.${extension} ]] || {
      fail "missing PLINK file: ${JOINT_BFILE}.${extension}"
    }
  done
  [[ -f $JOINT_BFILE.pop ]] || fail "missing file: $JOINT_BFILE.pop"

  fam_rows=$(wc -l <"$JOINT_BFILE.fam")
  pop_rows=$(wc -l <"$JOINT_BFILE.pop")
  [[ $fam_rows -eq $pop_rows ]] || {
    fail "joint FAM and POP row counts differ ($fam_rows vs $pop_rows)"
  }
fi

joint_directory=$(dirname "$JOINT_BFILE")
joint_input="$(basename "$JOINT_BFILE").bed"
printf 'RUN cd %s && %s --supervised %s %s\n' \
  "$joint_directory" "$ADMIXTURE" "$joint_input" "$K"
if [[ $DRY_RUN != 1 ]]; then
  (
    cd "$joint_directory"
    "$ADMIXTURE" --supervised "$joint_input" "$K"
  )
fi
