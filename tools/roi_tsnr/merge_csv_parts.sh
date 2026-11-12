#!/usr/bin/env bash
# Concatenate the per-task CSV fragments the job array wrote.
#   bash merge_csv_parts.sh <out_dir>
# Each task writes its own file because appending to one CSV from parallel
# tasks on a network filesystem loses lines. Keeps the first header only.
set -euo pipefail

OUT_DIR="${1:?usage: merge_csv_parts.sh <out_dir>}"
PARTS="${OUT_DIR}/csv_parts"
[[ -d "$PARTS" ]] || { echo "No such directory: $PARTS" >&2; exit 1; }

shopt -s nullglob

merge_one() {
    local name="$1"; shift
    (( $# )) || return 0
    awk 'FNR==1 && NR!=1 {next} {print}' "$@" > "${OUT_DIR}/${name}.csv"
    echo "  $(( $(wc -l < "${OUT_DIR}/${name}.csv") - 1 )) rows from $# parts -> ${OUT_DIR}/${name}.csv"
}

mains=( "${PARTS}"/roi_tsnr_*.csv )

# "${arr[@]}" on an EMPTY array is an unbound-variable error under set -u in
# bash before 4.4, and a dry run leaves it empty.
merge_one roi_tsnr ${mains[@]+"${mains[@]}"}
