#!/usr/bin/env bash
set -euo pipefail

ADMIXTURE=${ADMIXTURE:-admixture}
K=${K:-5}

joint_directory=$(dirname "$JOINT_BFILE")
joint_input="$(basename "$JOINT_BFILE").bed"
(
  cd "$joint_directory"
  "$ADMIXTURE" --supervised "$joint_input" "$K"
)
