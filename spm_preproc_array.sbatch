#!/usr/bin/env bash
# =============================================================================
#  SPM Preprocessing — SLURM Array Job
# =============================================================================
#  Submits one preprocessing job per subject-session.
#  All settings are read from pipeline_config.cfg — no editing needed here.
#
#  Usage:
#    N=$(wc -l < subses_list.txt)
#    sbatch --array=1-"$N" spm_preproc_array.sbatch
#
#  Or simply use run_pipeline.sh, which handles this for you.
# =============================================================================

# --- SLURM defaults (overridden at submission via run_pipeline.sh) ---
#SBATCH --job-name=spm_preproc
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# ---- Resolve script directory ----
# PIPELINE_SCRIPT_DIR is passed via --export when submitted from run_pipeline.sh.
# Falls back to auto-detect for manual runs (won't work inside SLURM spool).
SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    echo "       Make sure pipeline_config.cfg is in the same folder as this script."
    exit 1
fi

source "$CONFIG_FILE"

# ---- Validate required config variables ----
: "${BIDS_ROOT:?ERROR: BIDS_ROOT not set in pipeline_config.cfg}"
: "${SCRIPT_DIR:?ERROR: SCRIPT_DIR not set in pipeline_config.cfg}"
: "${MATLAB_MODULE:=matlab}"

# ---- Create logs directory ----
mkdir -p "${SCRIPT_DIR}/logs"

# ---- Determine subject-session list ----
LIST="${SCRIPT_DIR}/subses_list.txt"
if [[ ! -f "$LIST" ]]; then
    echo "ERROR: Subject list not found: $LIST"
    echo "       Run step0_make_subses_list.sh first."
    exit 1
fi

# ---- Read the Nth line (SLURM array task ID) ----
SUB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | awk '{print $1}')
SES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | awk '{print $2}')

if [[ -z "${SUB:-}" || -z "${SES:-}" ]]; then
    echo "No SUB/SES for array task ${SLURM_ARRAY_TASK_ID}"
    echo "  (check --array range vs number of lines in subses_list.txt)"
    exit 1
fi

echo "============================================"
echo "  SPM Preprocessing"
echo "  Subject: $SUB"
echo "  Session: $SES"
echo "  Config:  $CONFIG_FILE"
echo "============================================"
echo ""

# ---- Load modules ----
module purge
module load "${MATLAB_MODULE}"

# Load FSL if topup mode is selected
if [[ "${PREPROC_MODE:-}" == "topup" ]]; then
    FSL_MODULE="${FSL_MODULE:-fsl}"
    echo "Loading FSL module: ${FSL_MODULE}"
    module load "${FSL_MODULE}"
    echo "FSL directory: ${FSLDIR:-not set}"
fi

# ---- Run preprocessing ----
# Pass the config file path so MATLAB reads all settings from it
matlab -nodisplay -nosplash -r "\
    addpath('${SCRIPT_DIR}'); \
    try, \
        run_spm_preproc('${BIDS_ROOT}', '${SUB}', '${SES}', '${CONFIG_FILE}'); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
        exit(1); \
    end; \
    exit(0);"
