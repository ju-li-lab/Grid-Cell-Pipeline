#!/usr/bin/env bash
# =============================================================================
#  SPM Preprocessing — SLURM Array Job (also runnable directly)
# =============================================================================
#  Runs one preprocessing job per subject-session. All settings are read from
#  pipeline_config.cfg — no editing needed here.
#
#  Usage (SLURM array, one task per line of subses_list.txt):
#    N=$(wc -l < subses_list.txt)
#    sbatch --array=1-"$N" spm_preproc_array.sbatch
#    sbatch --array=1-"$N" spm_preproc_array.sbatch --stages realign,coreg,smooth
#
#  Usage (one subject, right here, no SLURM):
#    bash spm_preproc_array.sbatch --sub sub-01s13 --ses ses-01 --stages smooth
#
#  Or use submit_preproc.sh / run_pipeline.sh, which fill in the SLURM resource
#  flags from pipeline_config.cfg for you.
#
#  Output goes to <DERIV_ROOT>/<DERIV_PREPROC>/<sub>/<ses>/, never into
#  BIDS_ROOT. A stage list that starts the preprocessing from its beginning
#  clears that directory first; see DERIV_RESET in pipeline_config.cfg.
#
#  Options:
#    --stages LIST   Preprocessing stages to run (default: PREPROC_STAGES from
#                    the config, or "all"). See PREPROC_STAGES in the config for
#                    the stage list and what each one writes.
#    --sub SUB       Process this subject instead of reading the list.
#    --ses SES       Session that goes with --sub.
#    --list FILE     Subject-session list to index with SLURM_ARRAY_TASK_ID
#                    (default: <script dir>/subses_list.txt).
#    --config FILE   Config file (default: <script dir>/pipeline_config.cfg).
#    --task-id N     Line of the list to process (default: SLURM_ARRAY_TASK_ID).
#
#  Environment overrides (used when submitting via --export):
#    PIPELINE_SCRIPT_DIR, PIPELINE_CONFIG, PREPROC_STAGES, SUBSES_LIST
# =============================================================================

# --- SLURM defaults (overridden at submission via run_pipeline.sh) ---
#SBATCH --job-name=spm_preproc
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# ---- Locate the config file ----
# PIPELINE_SCRIPT_DIR is passed via --export when submitted from run_pipeline.sh.
# Falls back to auto-detect for manual runs (won't work inside SLURM spool).
BOOT_DIR="${PIPELINE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

CONFIG_FILE="${PIPELINE_CONFIG:-${BOOT_DIR}/pipeline_config.cfg}"
LIST="${SUBSES_LIST:-}"
STAGES="${PREPROC_STAGES:-}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-}"
SUB=""
SES=""

# ---- Parse command-line options ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stages)  STAGES="$2";      shift 2 ;;
        --sub)     SUB="$2";         shift 2 ;;
        --ses)     SES="$2";         shift 2 ;;
        --list)    LIST="$2";        shift 2 ;;
        --config)  CONFIG_FILE="$2"; shift 2 ;;
        --task-id) TASK_ID="$2";     shift 2 ;;
        -h|--help)
            awk 'NR>1 { if ($0 !~ /^#/) exit; print }' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "       Run with --help to see the accepted options."
            exit 1 ;;
    esac
done

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    echo "       Make sure pipeline_config.cfg is in the same folder as this script."
    exit 1
fi

# Sourcing sets SCRIPT_DIR (and PREPROC_STAGES) from the config file.
source "$CONFIG_FILE"

# Derivative path helpers, for the banner below. The MATLAB side resolves the
# same paths itself from the config, so this is only for the log.
if [[ -f "${BOOT_DIR}/bids_common.sh" ]]; then
    # shellcheck disable=SC1091
    source "${BOOT_DIR}/bids_common.sh"
fi

# ---- Validate required config variables ----
: "${BIDS_ROOT:?ERROR: BIDS_ROOT not set in pipeline_config.cfg}"
: "${SCRIPT_DIR:?ERROR: SCRIPT_DIR not set in pipeline_config.cfg}"
: "${MATLAB_MODULE:=matlab}"

# Subject list: --list / SUBSES_LIST wins, otherwise the one next to the scripts.
if [[ -z "$LIST" ]]; then
    LIST="${SCRIPT_DIR}/subses_list.txt"
fi

# Stage list: --stages wins, then the PREPROC_STAGES environment variable, then
# PREPROC_STAGES from the config (set by the source above), then "all".
if [[ -z "$STAGES" ]]; then
    STAGES="${PREPROC_STAGES:-all}"
fi
STAGES="${STAGES// /}"

# ---- Create logs directory ----
mkdir -p "${SCRIPT_DIR}/logs"

# ---- Determine which subject-session to process ----
if [[ -n "$SUB" || -n "$SES" ]]; then
    # Explicit subject/session — no list needed
    if [[ -z "$SUB" || -z "$SES" ]]; then
        echo "ERROR: --sub and --ses must be given together."
        exit 1
    fi
else
    if [[ ! -f "$LIST" ]]; then
        echo "ERROR: Subject list not found: $LIST"
        echo "       Run step0_make_subses_list.sh first, or pass --sub/--ses."
        exit 1
    fi

    if [[ -z "$TASK_ID" ]]; then
        echo "ERROR: No SLURM_ARRAY_TASK_ID and no --sub/--ses given."
        echo "       Submit as an array job, or pass --task-id N, or --sub/--ses."
        exit 1
    fi

    SUB=$(sed -n "${TASK_ID}p" "$LIST" | awk '{print $1}')
    SES=$(sed -n "${TASK_ID}p" "$LIST" | awk '{print $2}')

    if [[ -z "${SUB:-}" || -z "${SES:-}" ]]; then
        echo "No SUB/SES for array task ${TASK_ID}"
        echo "  (check --array range vs number of lines in ${LIST})"
        exit 1
    fi
fi

echo "============================================"
echo "  SPM Preprocessing"
echo "  Subject: $SUB"
echo "  Session: $SES"
echo "  Stages:  $STAGES"
echo "  Mode:    ${PREPROC_MODE:-<not set>}"
echo "  Config:  $CONFIG_FILE"
echo "  Reads:   ${BIDS_ROOT}/${SUB}/${SES}   (never written to)"
if declare -f deriv_session >/dev/null 2>&1; then
    echo "  Writes:  $(deriv_session "$(deriv_preproc)" "$SUB" "$SES" 2>/dev/null || echo '<derivatives root not resolvable>')"
fi
echo "============================================"
echo ""

# ---- Load modules ----
if command -v module >/dev/null 2>&1; then
    module purge
    module load "${MATLAB_MODULE}"

    # FSL is needed by topup estimation, by converting a topup field to a VDM,
    # and by applytopup during realignment. Only coreg and smooth never call it,
    # so FSL is loaded unless every selected stage is one of those two.
    if [[ "${PREPROC_MODE:-}" == "topup" ]] && [[ ! "$STAGES" =~ ^(coreg|smooth)(,(coreg|smooth))*$ ]]; then
        FSL_MODULE="${FSL_MODULE:-fsl}"
        echo "Loading FSL module: ${FSL_MODULE}"
        module load "${FSL_MODULE}"
        echo "FSL directory: ${FSLDIR:-not set}"
    else
        echo "Skipping FSL module (not needed for stages: ${STAGES})"
    fi
else
    echo "NOTE: no 'module' command available — using matlab/fsl from PATH."
    command -v matlab >/dev/null 2>&1 || {
        echo "ERROR: matlab not found on PATH and no module system available."
        exit 1
    }
fi

# ---- Run preprocessing ----
# Pass the config file path so MATLAB reads all settings from it, and the stage
# list so it knows which parts of the preprocessing to execute.
matlab -nodisplay -nosplash -r "\
    addpath('${SCRIPT_DIR}'); \
    try, \
        run_spm_preproc('${BIDS_ROOT}', '${SUB}', '${SES}', '${CONFIG_FILE}', '${STAGES}'); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
        exit(1); \
    end; \
    exit(0);"
