#!/usr/bin/env bash
# =============================================================================
#  Run the whole pipeline
# =============================================================================
#  Discovers subjects and filters them here, then queues the three cluster
#  stages as one SLURM dependency chain so they run unattended:
#
#      00_list_subjects.sh      (now, locally)
#         |
#      01_preprocess            -.
#      02_prepare_gridcat        |-- chained with --dependency=CHAIN_DEPENDENCY
#      03_gridcat                -'
#
#  Usage:
#    bash run_all.sh                  queue everything
#    bash run_all.sh --check          run the preflight checks and stop
#    bash run_all.sh --skip-list      keep the existing subses_list.txt
#
#  Settings come from config.sh. To run a single stage, call its script
#  directly - each one takes the same options on its own:
#
#    bash submit_preproc.sh --stages smooth
#    bash submit_preproc.sh --sub sub-01s13 --ses ses-01 --local
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SELF_DIR}/config.sh"

CHECK_ONLY=false
SKIP_LIST=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --check)     CHECK_ONLY=true; shift ;;
        --skip-list) SKIP_LIST=true;  shift ;;
        -h|--help)   sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config.sh not found at $CONFIG_FILE" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"
# shellcheck disable=SC1091
source "${SELF_DIR}/bids_common.sh"

# -----------------------------------------------------------------------------
#  Preflight
# -----------------------------------------------------------------------------
#  Only facts the code cannot know on its own: does this path exist, is that
#  directory writable, is the toolbox where the config says it is. Checking
#  that a setting is one of a list of allowed strings is left to the code
#  that consumes it, which is where the list actually lives.

fail=0
check_dir() {
    if [[ -d "$2" ]]; then
        printf '  ok    %-14s %s\n' "$1" "$2"
    else
        printf '  MISS  %-14s %s\n' "$1" "${2:-<not set>}"; fail=1
    fi
}
check_file() {
    if [[ -f "$2" ]]; then
        printf '  ok    %-14s %s\n' "$1" "$2"
    else
        printf '  MISS  %-14s %s\n' "$1" "${2:-<not found>}"; fail=1
    fi
}

echo "Checking configuration"
check_dir  BIDS_ROOT    "${BIDS_ROOT:-}"
check_dir  PIPELINE_DIR "${PIPELINE_DIR:-}"
check_dir  SPM_DIR      "${SPM_DIR:-}"
check_dir  GRIDCAT_DIR  "${GRIDCAT_DIR:-}"
check_file spm.m        "${SPM_DIR:-}/spm.m"
check_file specifyGLM.m "${GRIDCAT_DIR:-}/specifyGLM.m"

# OUTPUT_ROOT may not exist yet, but something above it must be writable.
out_probe="${OUTPUT_ROOT:-}"
while [[ -n "$out_probe" && ! -e "$out_probe" ]]; do out_probe="$(dirname "$out_probe")"; done
if [[ -w "$out_probe" ]]; then
    printf '  ok    %-14s %s (writable)\n' OUTPUT_ROOT "${OUTPUT_ROOT}"
else
    printf '  FAIL  %-14s %s is not writable\n' OUTPUT_ROOT "$out_probe"; fail=1
fi

# The one mistake that destroys data rather than just failing.
if [[ -n "${DERIV_ROOT:-}" && "${DERIV_ROOT}" == "${BIDS_ROOT}" ]]; then
    echo "  FAIL  DERIV_ROOT equals BIDS_ROOT - the pipeline would write into your raw data"
    fail=1
fi

if [[ "$fail" -ne 0 ]]; then
    echo
    echo "Fix the above in config.sh, then run again." >&2
    exit 1
fi
echo "  Mode: ${PREPROC_MODE:-<not set>}   Stages: ${PREPROC_STAGES:-all}   Variant: ${RUN_VARIANT:-<none>}"
echo "  Preprocessed BOLD will carry the prefix: $(func_prefix)"
echo

[[ "$CHECK_ONLY" == "true" ]] && { echo "Checks passed (--check given, nothing submitted)."; exit 0; }

# -----------------------------------------------------------------------------
#  Local stage: build the subject list
# -----------------------------------------------------------------------------
if [[ "$SKIP_LIST" == "true" ]]; then
    echo "Keeping the existing subject list (--skip-list)."
else
    bash "${SELF_DIR}/00_list_subjects.sh"
fi
echo

# -----------------------------------------------------------------------------
#  Cluster stages: queue as one dependency chain
# -----------------------------------------------------------------------------
CHAIN="${CHAIN_DEPENDENCY:-afterok}"

job_id_from() { grep -oE 'Submitted batch job [0-9]+' | grep -oE '[0-9]+' | tail -1; }

echo "Submitting preprocessing"
preproc_id=$(bash "${SELF_DIR}/submit_preproc.sh" | tee /dev/stderr | job_id_from)
[[ -n "$preproc_id" ]] || { echo "ERROR: preprocessing did not submit." >&2; exit 1; }

echo
echo "Submitting GridCAT preparation (waits ${CHAIN}:${preproc_id})"
prep_id=$(sbatch \
    --export=ALL,PIPELINE_SCRIPT_DIR="${SELF_DIR}" \
    --partition="${SLURM_PARTITION}" \
    --time="${PREP_TIME}" \
    --cpus-per-task="${PREP_CPUS}" \
    --mem="${PREP_MEM}" \
    --dependency="${CHAIN}:${preproc_id}" \
    "${SELF_DIR}/02_prepare_gridcat.sh" | tee /dev/stderr | job_id_from)
if [[ -z "$prep_id" ]]; then
    echo "ERROR: prep did not submit. Cancel the queued job with: scancel ${preproc_id}" >&2
    exit 1
fi

echo
echo "Submitting GridCAT analysis (waits ${CHAIN}:${prep_id})"
gridcat_id=$(GRIDCAT_DEPENDENCY="${CHAIN}:${prep_id}" \
    bash "${SELF_DIR}/03_gridcat.sh" | tee /dev/stderr | job_id_from)
if [[ -z "$gridcat_id" ]]; then
    echo "ERROR: GridCAT did not submit. Cancel with: scancel ${preproc_id} ${prep_id}" >&2
    exit 1
fi

cat <<SUMMARY

Queued as a dependency chain - SLURM starts each stage when the previous
one finishes, so there is nothing to come back and launch.

  preprocessing  job ${preproc_id}
  prep           job ${prep_id}   (waits ${CHAIN}:${preproc_id})
  GridCAT GLM    job ${gridcat_id}   (waits ${CHAIN}:${prep_id})

  watch:   squeue -u ${USER}
  cancel:  scancel ${preproc_id} ${prep_id} ${gridcat_id}
  results: ${OUTPUT_ROOT}/GLM_output${RUN_VARIANT:+_${RUN_VARIANT}}/
SUMMARY
