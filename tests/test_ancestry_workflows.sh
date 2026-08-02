#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

fail() {
  printf 'FAIL: %s\n' "$1" >&2
  exit 1
}

assert_contains() {
  local haystack=$1
  local needle=$2
  [[ $haystack == *"$needle"* ]] || fail "expected output to contain: $needle"
}

reference_script="$repo_root/workflow/ancestry/global/01_prepare_1000g_reference.sh"
[[ -f $reference_script ]] || fail "missing reference workflow: $reference_script"

reference_output=$(
  DRY_RUN=1 \
  REFERENCE_BFILE=/secure/1000g/all_phase3 \
  WORK_DIR=/tmp/ancestry-reference \
  bash "$reference_script"
)

assert_contains "$reference_output" "--autosome"
assert_contains "$reference_output" "--allow-extra-chr"
assert_contains "$reference_output" "--biallelic-only strict"
assert_contains "$reference_output" "--snps-only just-acgt"
assert_contains "$reference_output" "--maf 0.01"
assert_contains "$reference_output" "--indep-pairwise 1000 100 0.2"
assert_contains "$reference_output" "reference_ld_pruned.bed 5"

printf 'PASS: reference global-ancestry dry run\n'

RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}
selector_script="$repo_root/workflow/ancestry/global/02_select_reference_samples.R"
[[ -f $selector_script ]] || fail "missing reference selector: $selector_script"

fixture_dir=$(mktemp -d)
trap 'rm -rf "$fixture_dir"' EXIT

cat >"$fixture_dir/reference.fam" <<'EOF'
F1 S1 0 0 1 -9
F2 S2 0 0 2 -9
F3 S3 0 0 1 -9
F4 S4 0 0 2 -9
F5 S5 0 0 1 -9
F6 S6 0 0 2 -9
EOF

cat >"$fixture_dir/reference.Q" <<'EOF'
0.90 0.03 0.02 0.03 0.02
0.05 0.80 0.05 0.05 0.05
0.04 0.81 0.05 0.05 0.05
0.02 0.02 0.92 0.02 0.02
0.01 0.01 0.01 0.95 0.02
0.03 0.03 0.03 0.03 0.88
EOF

cat >"$fixture_dir/reference.psam" <<'EOF'
#FID IID SuperPop
F1 S1 AFR
F2 S2 AMR
F3 S3 AMR
F4 S4 EAS
F5 S5 EUR
F6 S6 SAS
EOF

"$RSCRIPT_BIN" "$selector_script" \
  "$fixture_dir/reference.Q" \
  "$fixture_dir/reference.fam" \
  "$fixture_dir/reference.psam" \
  "$fixture_dir/selected" \
  0.80

expected_keep=$'F3\tS3\nF1\tS1\nF5\tS5\nF6\tS6\nF4\tS4'
actual_keep=$(cat "$fixture_dir/selected.keep")
[[ $actual_keep == "$expected_keep" ]] || fail "reference keep file did not preserve the historical superpopulation order"

expected_labels=$'S3\tAMR\nS1\tAFR\nS5\tEUR\nS6\tSAS\nS4\tEAS'
actual_labels=$(cat "$fixture_dir/selected.labels.tsv")
[[ $actual_labels == "$expected_labels" ]] || fail "reference labels were not written in the historical superpopulation order"

printf 'PASS: strict reference selection\n'

joint_script="$repo_root/workflow/ancestry/global/03_prepare_joint_dataset.sh"
[[ -f $joint_script ]] || fail "missing joint-dataset workflow: $joint_script"

joint_output=$(
  DRY_RUN=1 \
  REFERENCE_BFILE=/secure/1000g/all_phase3 \
  STUDY_VCF=/secure/study/cohort_aims.vcf.gz \
  AIM_LIST=/public/ancestry_informative_snps.txt \
  REFERENCE_KEEP=/secure/reference.keep \
  WORK_DIR=/tmp/ancestry-joint \
  bash "$joint_script"
)

assert_contains "$joint_output" "--extract /public/ancestry_informative_snps.txt"
assert_contains "$joint_output" "--allow-extra-chr"
assert_contains "$joint_output" "--indiv-sort f /secure/reference.keep"
assert_contains "$joint_output" '--set-all-var-ids @:#:$r:$a'
assert_contains "$joint_output" "comm -12"
assert_contains "$joint_output" "--bmerge /tmp/ancestry-joint/study_shared"
assert_contains "$joint_output" "--out /tmp/ancestry-joint/joint"

printf 'PASS: joint ancestry dataset dry run\n'

pop_script="$repo_root/workflow/ancestry/global/04_build_supervised_pop.R"
[[ -f $pop_script ]] || fail "missing supervised-pop builder: $pop_script"

cat >"$fixture_dir/joint.fam" <<'EOF'
F1 S1 0 0 1 -9
F3 S3 0 0 1 -9
F4 S4 0 0 2 -9
F5 S5 0 0 1 -9
F6 S6 0 0 2 -9
T1 T1 0 0 0 -9
T2 T2 0 0 0 -9
EOF

"$RSCRIPT_BIN" "$pop_script" \
  "$fixture_dir/joint.fam" \
  "$fixture_dir/selected.labels.tsv" \
  "$fixture_dir/joint.pop"

expected_pop=$'AFR\nAMR\nEAS\nEUR\nSAS\n-\n-'
actual_pop=$(cat "$fixture_dir/joint.pop")
[[ $actual_pop == "$expected_pop" ]] || fail "supervised .pop labels do not match joint FAM order"

printf 'PASS: supervised population labels\n'

supervised_script="$repo_root/workflow/ancestry/global/05_run_supervised_admixture.sh"
[[ -f $supervised_script ]] || fail "missing supervised ADMIXTURE workflow: $supervised_script"

supervised_output=$(
  DRY_RUN=1 \
  JOINT_BFILE=/tmp/ancestry-joint/joint \
  bash "$supervised_script"
)

assert_contains "$supervised_output" "--supervised"
assert_contains "$supervised_output" "joint.bed 5"

printf 'PASS: supervised ADMIXTURE dry run\n'

query_script="$repo_root/workflow/ancestry/local/01_prepare_query_vcf.sh"
[[ -f $query_script ]] || fail "missing local-ancestry query workflow: $query_script"

query_output=$(
  DRY_RUN=1 \
  QUERY_BFILE=/secure/study/cohort_wgs \
  WORK_DIR=/tmp/local-ancestry-query \
  bash "$query_script"
)

assert_contains "$query_output" "--geno 0.05"
assert_contains "$query_output" "--maf 0.01"
assert_contains "$query_output" "--snps-only"
assert_contains "$query_output" "--recode vcf bgz"
assert_contains "$query_output" "tabix -f -p vcf /tmp/local-ancestry-query/query_common.vcf.gz"

printf 'PASS: local-ancestry query preparation dry run\n'

rfmix_reference_script="$repo_root/workflow/ancestry/local/02_prepare_rfmix_reference.sh"
[[ -f $rfmix_reference_script ]] || fail "missing RFMix reference workflow: $rfmix_reference_script"

rfmix_reference_output=$(
  DRY_RUN=1 \
  REFERENCE_VCF=/secure/1000g/phase3_phased.vcf.gz \
  REFERENCE_LABELS=/secure/reference.labels.tsv \
  WORK_DIR=/tmp/rfmix-reference \
  bash "$rfmix_reference_script"
)

assert_contains "$rfmix_reference_output" "view -q 0.01"
assert_contains "$rfmix_reference_output" "view -m2 -M2 -v snps"
assert_contains "$rfmix_reference_output" "view -S /tmp/rfmix-reference/reference_sample_ids.txt"
assert_contains "$rfmix_reference_output" "tabix -f -p vcf"

printf 'PASS: RFMix reference preparation dry run\n'

phase_script="$repo_root/workflow/ancestry/local/03_phase_query.sh"
[[ -f $phase_script ]] || fail "missing SHAPEIT4 workflow: $phase_script"

phase_output=$(
  DRY_RUN=1 \
  QUERY_VCF=/secure/study/query_common.vcf.gz \
  CHROMOSOME=1 \
  GENETIC_MAP_PATTERN=/public/maps/chr%s.b37.gmap.gz \
  WORK_DIR=/tmp/local-ancestry-phased \
  bash "$phase_script"
)

assert_contains "$phase_output" "--map /public/maps/chr1.b37.gmap.gz"
assert_contains "$phase_output" "--region 1"
assert_contains "$phase_output" "--output /tmp/local-ancestry-phased/phased_chr1.vcf.gz"
assert_contains "$phase_output" "--thread 8"
assert_contains "$phase_output" "--sequencing"

printf 'PASS: SHAPEIT4 phasing dry run\n'

rfmix_script="$repo_root/workflow/ancestry/local/04_run_rfmix.sh"
[[ -f $rfmix_script ]] || fail "missing RFMix inference workflow: $rfmix_script"

rfmix_output=$(
  DRY_RUN=1 \
  PHASED_QUERY_PATTERN=/secure/study/phased_chr%s.vcf.gz \
  REFERENCE_VCF=/secure/reference/rfmix_reference.vcf.gz \
  REFERENCE_SAMPLE_MAP=/secure/reference/reference_sample_map.txt \
  GENETIC_MAP=/public/maps/genetic_map_b37.txt \
  CHROMOSOME=1 \
  WORK_DIR=/tmp/rfmix-output \
  bash "$rfmix_script"
)

assert_contains "$rfmix_output" "-f /secure/study/phased_chr1.vcf.gz"
assert_contains "$rfmix_output" "-r /secure/reference/rfmix_reference.vcf.gz"
assert_contains "$rfmix_output" "-m /secure/reference/reference_sample_map.txt"
assert_contains "$rfmix_output" "-g /public/maps/genetic_map_b37.txt"
assert_contains "$rfmix_output" "--chromosome=1"
assert_contains "$rfmix_output" "--n-threads=8"

printf 'PASS: RFMix inference dry run\n'

while IFS= read -r shell_file; do
  bash -n "$shell_file"
done < <(find "$repo_root/workflow" "$repo_root/tests" -type f -name '*.sh' | sort)

while IFS= read -r r_file; do
  "$RSCRIPT_BIN" -e 'parse(file = commandArgs(trailingOnly = TRUE)[[1L]])' "$r_file" >/dev/null
done < <(find "$repo_root/R" "$repo_root/analysis" "$repo_root/workflow" -type f -name '*.R' | sort)

printf 'PASS: shell and R syntax checks\n'
