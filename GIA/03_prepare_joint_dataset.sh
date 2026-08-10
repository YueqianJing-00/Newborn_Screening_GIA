#!/usr/bin/env bash
set -euo pipefail

reference_bfile=${1}
msk_impact_snp_list=${2}
reference_keep=${3}
study_input_option=${4}
study_input=${5}
variant_id_template=${6}
selected_reference_impact=${7}
reference_canonical=${8}
study_canonical=${9}
reference_ids=${10}
study_ids=${11}
shared_ids=${12}
reference_shared=${13}
study_shared=${14}
joint_prefix=${15}

plink \
  --bfile "$reference_bfile" \
  --allow-extra-chr \
  --extract "$msk_impact_snp_list" \
  --keep "$reference_keep" \
  --indiv-sort f "$reference_keep" \
  --keep-allele-order \
  --make-bed \
  --out "$selected_reference_impact"

plink2 \
  --bfile "$selected_reference_impact" \
  --set-all-var-ids "$variant_id_template" \
  --make-bed \
  --out "$reference_canonical"

plink2 \
  "$study_input_option" "$study_input" \
  --set-all-var-ids "$variant_id_template" \
  --make-bed \
  --out "$study_canonical"

cut -f2 "$reference_canonical.bim" | LC_ALL=C sort -u >"$reference_ids"
cut -f2 "$study_canonical.bim" | LC_ALL=C sort -u >"$study_ids"
comm -12 "$reference_ids" "$study_ids" >"$shared_ids"

plink \
  --bfile "$reference_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$reference_shared"

plink \
  --bfile "$study_canonical" \
  --extract "$shared_ids" \
  --keep-allele-order \
  --make-bed \
  --out "$study_shared"

plink \
  --bfile "$reference_shared" \
  --bmerge "$study_shared" \
  --indiv-sort 0 \
  --keep-allele-order \
  --make-bed \
  --out "$joint_prefix"
