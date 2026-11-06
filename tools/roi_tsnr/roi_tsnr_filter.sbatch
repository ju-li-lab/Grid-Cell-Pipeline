#!/usr/bin/env bash
# =============================================================================
#  roi_tsnr_filter.sbatch - one SLURM array task per subject x SESSION
# =============================================================================
#  Do not call this directly; run_roi_tsnr_filter.sh builds the work list,
#  sizes the array and submits the merge. Call it by hand only to rerun one
#  task, with WORK_LIST set and SLURM_ARRAY_TASK_ID naming the line.
#
#  Sessions are the unit of parallelism, as asked: each one reads its own three
#  4D runs and its own masks, and shares nothing with any other. Each task
#  writes its own CSV under <out>/csv_parts/ - appending to one CSV from
#  parallel tasks on a network filesystem loses lines.
#
#  Memory: one 4D run is 96x96x30x347 single = ~380 MB, held one at a time,
#  plus a 20000-voxel working block. 8G is comfortable.
# =============================================================================
#SBATCH --job-name=roitsnr
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

BOOT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
CONFIG_FILE="${PIPELINE_CONFIG:-${BOOT_DIR}/config.sh}"
# shellcheck disable=SC1090
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE"

if [[ -f "${BOOT_DIR}/bids_common.sh" ]]; then
    # shellcheck disable=SC1090
    source "${BOOT_DIR}/bids_common.sh"
    load_modules "${MATLAB_MODULE:-matlab}"
elif command -v module >/dev/null 2>&1; then
    module purge
    module load "${MATLAB_MODULE:-matlab}"
fi

: "${DERIV_ROOT_IN:?not set - submit through run_roi_tsnr_filter.sh}"
: "${ROI_DIR:?not set - submit through run_roi_tsnr_filter.sh}"
: "${OUT_DIR:?not set - submit through run_roi_tsnr_filter.sh}"
: "${WORK_LIST:?not set - submit through run_roi_tsnr_filter.sh}"

# Defaults here, not inline: a ${VAR:-{'x'}} default holding braces is
# mis-parsed by the shell, and the MATLAB cell syntax is all braces.
#
# NAMES MUST NOT COLLIDE WITH config.sh. This script inherits the launcher's
# exports and THEN sources config.sh, so a shared name is silently overwritten
# by the config value. RTF_TASKS is prefixed for exactly that reason: config.sh
# defines TASKS=1,2,3, which arrives here as "1,2,3" and is not MATLAB syntax.
: "${KEEP_PCT:=[50 50]}"
: "${MIN_TSNR:=[]}"
: "${FUNC_PREFIX:=u}"
: "${RTF_TASKS:=[1 2 3]}"
: "${COMBINE:=mean}"
: "${HPF_SEC:=128}"
: "${TR_SEC:=2.06}"
: "${MIN_VOXELS:=5}"
: "${DRY_RUN:=false}"
: "${TRIAL:=$(date +%Y%m%d_%H%M%S)}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
LINE="$(sed -n "${TASK_ID}p" "$WORK_LIST")"
[[ -n "$LINE" ]] || { echo "No work on line ${TASK_ID} of ${WORK_LIST}" >&2; exit 1; }
SUB="${LINE%% *}"
SES="${LINE##* }"

SELF_DIR="${ROI_TSNR_DIR:-${BOOT_DIR}/tools/roi_tsnr}"
PARTS="${OUT_DIR}/csv_parts"
mkdir -p "$PARTS"

echo "task ${TASK_ID}: ${SUB} ${SES}"
echo "  functionals : ${DERIV_ROOT_IN} (${FUNC_PREFIX}sub-*)"
echo "  ROI masks   : ${ROI_DIR}"
echo "  keep        : ${KEEP_PCT}   min tSNR: ${MIN_TSNR}"
echo "  out         : ${OUT_DIR}"

matlab -nodisplay -nosplash -singleCompThread -batch "\
    addpath('${SELF_DIR}'); \
    try, \
        roi_tsnr_filter('${DERIV_ROOT_IN}', \
            'RoiDir',    '${ROI_DIR}', \
            'OutDir',    '${OUT_DIR}', \
            'Subjects',  {'${SUB}'}, \
            'Sessions',  {'${SES}'}, \
            'KeepPct',   ${KEEP_PCT}, \
            'MinTSNR',   ${MIN_TSNR}, \
            'Prefix',    '${FUNC_PREFIX}', \
            'Tasks',     ${RTF_TASKS}, \
            'Combine',   '${COMBINE}', \
            'HPFSec',    ${HPF_SEC}, \
            'TR',        ${TR_SEC}, \
            'MinVoxels', ${MIN_VOXELS}, \
            'DryRun',    ${DRY_RUN}, \
            'Trial',     '${TRIAL}', \
            'CSV',       '${PARTS}/roi_tsnr_${SUB}_${SES}.csv'); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
        exit(1); \
    end"

echo "task ${TASK_ID} done: ${SUB} ${SES}"
