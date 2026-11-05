#!/usr/bin/env bash
# Concatenate the per-subject CSV fragments that the job array wrote.
#
#   bash merge_csv_parts.sh <out_dir>
#
# Each array task writes its own file under <out_dir>/csv_parts/ because
# appending to one CSV from parallel tasks on a network filesystem loses
# lines. This keeps the first header and drops the rest.
set -euo pipefail

OUT_DIR="${1:?usage: merge_csv_parts.sh <out_dir>}"
PARTS="${OUT_DIR}/csv_parts"
[[ -d "$PARTS" ]] || { echo "No such directory: $PARTS" >&2; exit 1; }

shopt -s nullglob

merge_one() {   # merge_one <target name> <file>...
    local name="$1"; shift
    (( $# )) || return 0
    awk 'FNR==1 && NR!=1 {next} {print}' "$@" > "${OUT_DIR}/${name}.csv"
    echo "  $(( $(wc -l < "${OUT_DIR}/${name}.csv") - 1 )) rows from $# parts -> ${OUT_DIR}/${name}.csv"
}

maps=( "${PARTS}"/maps_*.csv )
summaries=( "${PARTS}"/clusters_*_summary.csv )

# every clusters_* that is not a _summary
clusters=()
for f in "${PARTS}"/clusters_*.csv; do
    [[ "$f" == *_summary.csv ]] || clusters+=( "$f" )
done

# "${arr[@]}" on an EMPTY array is an unbound-variable error under set -u in
# bash before 4.4, and --what maps leaves one of these empty by design.
merge_one maps              ${maps[@]+"${maps[@]}"}
merge_one clusters          ${clusters[@]+"${clusters[@]}"}
merge_one clusters_summary  ${summaries[@]+"${summaries[@]}"}
