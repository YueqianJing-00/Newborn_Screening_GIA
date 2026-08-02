#!/usr/bin/env bash

DRY_RUN=${DRY_RUN:-0}

die() {
  printf 'ERROR: %s\n' "$1" >&2
  exit 1
}

require_env() {
  local name=$1
  [[ -n ${!name:-} ]] || die "$name is required"
}

require_command() {
  local command_name=$1
  [[ $DRY_RUN == 1 ]] || command -v "$command_name" >/dev/null 2>&1 || {
    die "required command not found: $command_name"
  }
}

require_file() {
  local path=$1
  [[ $DRY_RUN == 1 || -f $path ]] || die "missing file: $path"
}

require_bfile() {
  local prefix=$1
  [[ $DRY_RUN == 1 ]] && return
  local extension
  for extension in bed bim fam; do
    [[ -f ${prefix}.${extension} ]] || die "missing PLINK file: ${prefix}.${extension}"
  done
}

print_command() {
  printf 'RUN'
  printf ' %s' "$@"
  printf '\n'
}

run_command() {
  print_command "$@"
  [[ $DRY_RUN == 1 ]] || "$@"
}

run_in_directory() {
  local directory=$1
  shift
  printf 'RUN cd %s &&' "$directory"
  printf ' %s' "$@"
  printf '\n'
  if [[ $DRY_RUN != 1 ]]; then
    (
      cd "$directory"
      "$@"
    )
  fi
}

prepare_output_directory() {
  local directory=$1
  [[ $DRY_RUN == 1 ]] || mkdir -p "$directory"
}
