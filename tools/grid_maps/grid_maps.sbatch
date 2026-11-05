#!/usr/bin/env bash
# =============================================================================
#  grid_maps.sbatch - voxel-wise grid maps, one SLURM array task per subject
# =============================================================================
#  Do not call this directly; run_grid_maps.sh builds the subject list, sizes
#  the array and submits the merge job. Call it by hand only to rerun one
#  task, with SUBJECT_LIST pointing at a file of subject ids.
#
#  One subject per task because that is the natural grain: nothing in
#  grid_maps or grid_clusters crosses a subject, and each task is one MATLAB
#  start (~30 s) plus a few seconds of arithmetic. Six subjects finish in the
#  time one takes.
#
#  Each task writes its OWN CSV under <out>/csv_parts/. Appending to one CSV
#  from parallel tasks on a network filesystem is a lost-line bug waiting to
#  happen, so the merge is a separate dependent job.
# =============================================================================
#SBATCH --job-name=gridmaps
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
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

: "${GLM_DIR:?GLM_DIR not set - submit through run_grid_maps.sh}"
: "${OUT_DIR:?OUT_DIR not set - submit through run_grid_maps.sh}"
: "${SUBJECT_LIST:?SUBJECT_LIST not set - submit through run_grid_maps.sh}"

# Defaults set here, not inline: a ${VAR:-{'x'}} default holding braces is
# mis-parsed by the shell, and the MATLAB cell syntax is all braces.
: "${WHAT:=both}"
: "${SESSIONS_ARG:={'ses-01'}}"
: "${MASK_DIR:=}"
: "${ANAT_ROOT:=}"
: "${RADIUS_MM:=5}"
: "${CLUST_RUN:=0}"
: "${ANGLE_TOL:=5}"
: "${MIN_SIZE:=10}"
: "${PERMS:=500}"
: "${TRIAL:=$(date +%Y%m%d_%H%M%S)}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
SUB="$(sed -n "${TASK_ID}p" "$SUBJECT_LIST")"
[[ -n "$SUB" ]] || { echo "No subject on line ${TASK_ID} of ${SUBJECT_LIST}" >&2; exit 1; }

SELF_DIR="${GRID_MAPS_DIR:-${BOOT_DIR}/tools/grid_maps}"
PARTS="${OUT_DIR}/csv_parts"
mkdir -p "$PARTS"

# MATLAB spawns threads per core by default; the array already uses the node.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

echo "task ${TASK_ID}: ${SUB}"
echo "  GLM  : ${GLM_DIR}"
echo "  out  : ${OUT_DIR}"
echo "  what : ${WHAT}"

run_matlab() {
    matlab -nodisplay -nosplash -singleCompThread -batch "\
        addpath('${SELF_DIR}'); \
        try, $1 catch ME, \
            fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
            exit(1); \
        end"
}

if [[ "$WHAT" == "maps" || "$WHAT" == "both" ]]; then
    run_matlab "grid_maps('${GLM_DIR}', \
        'Subjects', {'${SUB}'}, \
        'Sessions', ${SESSIONS_ARG}, \
        'MaskDir',  '${MASK_DIR}', \
        'AnatRoot', '${ANAT_ROOT}', \
        'OutDir',   '${OUT_DIR}/maps', \
        'RadiusMM', ${RADIUS_MM}, \
        'Trial',    '${TRIAL}', \
        'CSV',      '${PARTS}/maps_${SUB}.csv');"
fi

if [[ "$WHAT" == "clusters" || "$WHAT" == "both" ]]; then
    run_matlab "grid_clusters('${GLM_DIR}', \
        'Subjects', {'${SUB}'}, \
        'Sessions', ${SESSIONS_ARG}, \
        'MaskDir',  '${MASK_DIR}', \
        'OutDir',   '${OUT_DIR}/clusters', \
        'Run',      ${CLUST_RUN}, \
        'AngleTol', ${ANGLE_TOL}, \
        'MinSize',  ${MIN_SIZE}, \
        'Perms',    ${PERMS}, \
        'Trial',    '${TRIAL}', \
        'CSV',      '${PARTS}/clusters_${SUB}.csv');"
fi

echo "task ${TASK_ID} done: ${SUB}"
