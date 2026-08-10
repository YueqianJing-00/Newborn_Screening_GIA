#!/usr/bin/env bash
set -euo pipefail

# Edit this path and parameter before running.
joint_output_dir="/secure/work/global_joint"
joint_bfile="${joint_output_dir}/joint"
k=5

joint_directory=$(dirname "$joint_bfile")
joint_input="$(basename "$joint_bfile").bed"
(
  cd "$joint_directory"
  admixture --supervised "$joint_input" "$k"
)
