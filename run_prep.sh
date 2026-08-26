#!/usr/bin/env bash
# =============================================================================
#  GridCAT Data Preparation — SLURM Job
# =============================================================================
#  Prepares preprocessed data for GridCAT analysis:
#    - Splits 4D functional files into 3D volumes
#    - Copies motion regressors and event tables
#    - Creates bilateral ROI masks
#
#  Reads the preprocessing derivatives (<DERIV_ROOT>/<DERIV_PREPROC>), not the
#  raw BIDS directory, and drops the run entity from every name it writes so
#  the event tables line up. Which run each file came from is recorded in
#  GLM_runauto/run_manifest.tsv.
#
#  All settings are read from pipeline_config.cfg — no editing needed here.
#
#  Usage:
#    sbatch run_prep.sh
#
#  Or simply use run_pipeline.sh, which handles this for you.
# =============================================================================

# --- SLURM settings ---
#SBATCH --job-name=gridcat_prep
#SBATCH --partition=compute
#SBATCH --output=logs/gridcat_prep_%j.out
#SBATCH --error=logs/gridcat_prep_%j.err

set -euo pipefail

# ---- Resolve script directory ----
# PIPELINE_SCRIPT_DIR is passed via --export when submitted from run_pipeline.sh.
# Falls back to auto-detect for manual runs (won't work inside SLURM spool).
SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

source "$CONFIG_FILE"
source "${SCRIPT_DIR}/bids_common.sh"

# ---- Validate required config variables ----
: "${BIDS_ROOT:?ERROR: BIDS_ROOT not set in pipeline_config.cfg}"
: "${OUTPUT_ROOT:?ERROR: OUTPUT_ROOT not set in pipeline_config.cfg}"
: "${SPM_DIR:?ERROR: SPM_DIR not set in pipeline_config.cfg}"
: "${MATLAB_MODULE:=matlab}"

PREPROC_ROOT="$(deriv_preproc)"
if [[ ! -d "$PREPROC_ROOT" ]]; then
    echo "ERROR: no preprocessing derivatives at: $PREPROC_ROOT"
    echo "       Run the preprocessing first (menu option 2, or submit_preproc.sh)."
    exit 1
fi

# ---- Create logs directory ----
mkdir -p "${SCRIPT_DIR}/logs"

echo "============================================"
echo "  GridCAT Data Preparation"
echo "  Derivatives: $PREPROC_ROOT"
echo "  Output Root: $OUTPUT_ROOT"
echo "  Config:      $CONFIG_FILE"
echo "============================================"
echo ""

# ---- Load MATLAB ----
module purge
module load "${MATLAB_MODULE}"

# ---- Determine functional file suffix based on preprocessing mode ----
# In topup+applytopup mode, SPM writes u*_bold_dc.nii instead of u*_bold.nii
FUNC_SUFFIX="_bold"
if [[ "${PREPROC_MODE:-}" == "topup" && "${TOPUP_APPLY_METHOD:-}" == "applytopup" ]]; then
    FUNC_SUFFIX="_bold_dc"
    echo "  Func suffix: _bold_dc (topup+applytopup mode)"
else
    echo "  Func suffix: _bold (standard mode)"
fi

# ---- Build exclude tasks list for MATLAB ----
# Convert comma-separated EXCLUDE_TASKS to MATLAB cell array
EXCLUDE="${EXCLUDE_TASKS:-}"

# ---- Determine subject list path ----
SUBSES_LIST="${SCRIPT_DIR}/subses_list.txt"
if [[ -f "$SUBSES_LIST" ]]; then
    echo "  Subject list: $SUBSES_LIST ($(wc -l < "$SUBSES_LIST") entries)"
    SUBJECT_LIST_ARG="'SubjectList', '${SUBSES_LIST}',"
else
    echo "  No subses_list.txt found — will scan BIDS directory"
    SUBJECT_LIST_ARG=""
fi

matlab -nodisplay -nosplash -r "\
    addpath('${SCRIPT_DIR}'); \
    try, \
        excludeStr = '${EXCLUDE}'; \
        if ~isempty(excludeStr), \
            parts = strsplit(excludeStr, ','); \
            excludeTasks = cellfun(@(x) ['task-' strtrim(x)], parts, 'UniformOutput', false); \
        else, \
            excludeTasks = {}; \
        end; \
        prepare_gridcat_directory( \
            '${PREPROC_ROOT}', \
            '${OUTPUT_ROOT}', \
            'ExcludeTasks', excludeTasks, \
            'FuncPrefix', '${FUNC_PREFIX:-u}', \
            'FuncSuffix', '${FUNC_SUFFIX}', \
            'ROIPrefix', '${ROI_PREFIX:-r}', \
            'CopyMode', '${COPY_MODE:-copy}', \
            'SPMPath', '${SPM_DIR}', \
            'MaxWorkers', ${PREP_CPUS:-5}, \
            ${SUBJECT_LIST_ARG} \
            'Verbose', true \
        ); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
        exit(1); \
    end; \
    exit(0);"

echo ""
echo "Data preparation complete."
echo "Output: ${OUTPUT_ROOT}/GLM_runauto/"
