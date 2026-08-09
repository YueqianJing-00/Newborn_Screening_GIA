#!/usr/bin/env bash
set -euo pipefail

joint_bfile=$1
k=$2

joint_directory=$(dirname "$joint_bfile")
joint_input="$(basename "$joint_bfile").bed"
(
  cd "$joint_directory"
  admixture --supervised "$joint_input" "$k"
)
