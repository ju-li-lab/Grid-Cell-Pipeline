#!/usr/bin/env bash
# =============================================================================
#  GridCAT GLM Analysis — dual-mode launcher + SLURM job
# =============================================================================
#  Preferred usage:
#    bash run_gridcat.sh                       — submit the job (snapshots config)
#    sbatch run_gridcat.sh                     — also works; SLURM runs the
#                                                 submitter, which then submits
#                                                 the real analysis job.
#
#  On submission this script takes a timestamped snapshot of pipeline_config.cfg
#  and subses_list.txt into runs/gridcat_<variant>_<timestamp>/ and tells SLURM
#  to read from that snapshot. After submission you can freely edit
#  pipeline_config.cfg (e.g. set a new RUN_VARIANT and resubmit) without
#  disturbing the job that is already queued or running.
#
#  All settings are read from pipeline_config.cfg — no editing needed here.
# =============================================================================

# --- SLURM settings (used when this script runs AS the job) ---
#SBATCH --job-name=gridcat_glm
#SBATCH --partition=compute
#SBATCH --output=logs/gridcat_glm_%j.out
#SBATCH --error=logs/gridcat_glm_%j.err

set -euo pipefail

ORIG_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
#  SUBMITTER MODE — not yet inside a SLURM job
# =============================================================================
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    CONFIG_FILE="${ORIG_SCRIPT_DIR}/pipeline_config.cfg"
    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "ERROR: Config file not found: $CONFIG_FILE" >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "$CONFIG_FILE"

    : "${OUTPUT_ROOT:?ERROR: OUTPUT_ROOT not set in pipeline_config.cfg}"
    : "${SLURM_PARTITION:=compute}"
    : "${GRIDCAT_TIME:=01:00:00}"
    : "${GRIDCAT_CPUS:=30}"
    : "${GRIDCAT_MEM:=80G}"

    # Sanitize variant for filesystem use (allow letters, digits, _ - .)
    VARIANT_RAW="${RUN_VARIANT:-}"
    SAFE_VARIANT="$(printf '%s' "$VARIANT_RAW" | tr -c 'A-Za-z0-9._-' '_' | sed -E 's/_+/_/g; s/^_//; s/_$//')"

    TS="$(date +%Y%m%d_%H%M%S)"
    if [[ -n "$SAFE_VARIANT" ]]; then
        SNAP_NAME="gridcat_${SAFE_VARIANT}_${TS}"
    else
        SNAP_NAME="gridcat_${TS}"
    fi
    SNAP_DIR="${ORIG_SCRIPT_DIR}/runs/${SNAP_NAME}"

    mkdir -p "${SNAP_DIR}" "${ORIG_SCRIPT_DIR}/logs"

    cp "$CONFIG_FILE" "${SNAP_DIR}/pipeline_config.cfg"
    if [[ -f "${ORIG_SCRIPT_DIR}/subses_list.txt" ]]; then
        cp "${ORIG_SCRIPT_DIR}/subses_list.txt" "${SNAP_DIR}/subses_list.txt"
    fi

    {
        echo "snapshot_dir: ${SNAP_DIR}"
        echo "variant_raw:  ${VARIANT_RAW}"
        echo "variant_safe: ${SAFE_VARIANT}"
        echo "submitted_by: ${USER:-unknown}"
        echo "submitted_at: ${TS}"
        echo "orig_scripts: ${ORIG_SCRIPT_DIR}"
        echo "output_dir:   ${OUTPUT_ROOT}/GLM_output${SAFE_VARIANT:+_${SAFE_VARIANT}}"
    } > "${SNAP_DIR}/submission_info.txt"

    echo "============================================"
    echo "  GridCAT GLM — submitting job"
    echo "  Variant:   ${VARIANT_RAW:-<none>}"
    echo "  Snapshot:  ${SNAP_DIR}"
    echo "  Output:    ${OUTPUT_ROOT}/GLM_output${SAFE_VARIANT:+_${SAFE_VARIANT}}"
    echo "============================================"

    GRES_FLAG=()
    if [[ "${USE_LOCAL_SCRATCH:-false}" == "true" && -n "${GRIDCAT_LOCAL_TMP:-}" ]]; then
        GRES_FLAG=(--gres="tmp:${GRIDCAT_LOCAL_TMP}")
    fi

    # Optional SLURM dependency (e.g. afterok:1234) so run_pipeline.sh can chain
    # this GLM job after the prep job in a single "Run ALL steps" submission.
    DEP_FLAG=()
    if [[ -n "${GRIDCAT_DEPENDENCY:-}" ]]; then
        DEP_FLAG=(--dependency="${GRIDCAT_DEPENDENCY}")
    fi

    # Email once when the GridCAT analysis finishes (or fails).
    MAIL_FLAG=()
    if [[ -n "${NOTIFY_EMAIL:-}" ]]; then
        MAIL_FLAG=(--mail-user="${NOTIFY_EMAIL}" --mail-type=END,FAIL)
    fi

    exec sbatch \
        --export=ALL,PIPELINE_SCRIPT_DIR="${SNAP_DIR}",GRIDCAT_ORIG_SCRIPT_DIR="${ORIG_SCRIPT_DIR}" \
        --partition="${SLURM_PARTITION}" \
        --time="${GRIDCAT_TIME}" \
        --cpus-per-task="${GRIDCAT_CPUS}" \
        --mem="${GRIDCAT_MEM}" \
        "${DEP_FLAG[@]}" \
        "${MAIL_FLAG[@]}" \
        "${GRES_FLAG[@]}" \
        "${ORIG_SCRIPT_DIR}/run_gridcat.sh"
fi

# =============================================================================
#  EXECUTION MODE — running inside a SLURM job
# =============================================================================

# Resolve the snapshot directory that the submitter created.
SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-$ORIG_SCRIPT_DIR}"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

if [[ -z "${PIPELINE_SCRIPT_DIR:-}" ]]; then
    echo "WARNING: PIPELINE_SCRIPT_DIR not set — using live config at $CONFIG_FILE."
    echo "         Editing that file while this job runs may corrupt the run."
    echo "         Prefer: bash run_gridcat.sh   (handles snapshotting automatically)."
fi

# Fall back to ORIG_SCRIPT_DIR for logging directory (shared across runs).
LOG_DIR="${GRIDCAT_ORIG_SCRIPT_DIR:-$ORIG_SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

source "$CONFIG_FILE"

# ---- Validate required config variables ----
: "${SPM_DIR:?ERROR: SPM_DIR not set in pipeline_config.cfg}"
: "${GRIDCAT_DIR:?ERROR: GRIDCAT_DIR not set in pipeline_config.cfg}"
: "${OUTPUT_ROOT:?ERROR: OUTPUT_ROOT not set in pipeline_config.cfg}"
: "${MATLAB_MODULE:=matlab}"

# ---- Determine output suffix from RUN_VARIANT ----
VARIANT_RAW="${RUN_VARIANT:-}"
SAFE_VARIANT="$(printf '%s' "$VARIANT_RAW" | tr -c 'A-Za-z0-9._-' '_' | sed -E 's/_+/_/g; s/^_//; s/_$//')"
OUTPUT_SUFFIX=""
if [[ -n "$SAFE_VARIANT" ]]; then
    OUTPUT_SUFFIX="_${SAFE_VARIANT}"
fi

# ---- Determine subject list path ----
SUBSES_LIST="${SCRIPT_DIR}/subses_list.txt"
if [[ -f "$SUBSES_LIST" ]]; then
    echo "  Subject list: $SUBSES_LIST ($(wc -l < "$SUBSES_LIST") entries)"
    SUBJECT_LIST_ARG="'${SUBSES_LIST}'"
else
    echo "  No subses_list.txt found — will process all EventData files"
    SUBJECT_LIST_ARG="''"
fi

echo "============================================"
echo "  GridCAT GLM Analysis"
echo "  Variant: ${VARIANT_RAW:-<none>}"
echo "  Input:   ${OUTPUT_ROOT}/GLM_runauto"
echo "  Output:  ${OUTPUT_ROOT}/GLM_output${OUTPUT_SUFFIX}"
echo "  Config:  $CONFIG_FILE"
echo "============================================"
echo ""

# ---- Load MATLAB ----
module purge
module load "${MATLAB_MODULE}"

# ---- Run GridCAT analysis ----
# Pass config file and subject list so run_gridcat_analysis filters correctly.
# The MATLAB side finds its helper functions in GRIDCAT_ORIG_SCRIPT_DIR (the
# untouched checkout). Only config + subject list live in the snapshot.
MATLAB_CODE_DIR="${GRIDCAT_ORIG_SCRIPT_DIR:-$ORIG_SCRIPT_DIR}"

matlab -nodisplay -nosplash -r "\
    addpath('${MATLAB_CODE_DIR}'); \
    try, \
        run_gridcat_analysis('${CONFIG_FILE}', ${SUBJECT_LIST_ARG}); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
        exit(1); \
    end; \
    exit(0);"

echo ""
echo "GridCAT analysis complete."
echo "Results: ${OUTPUT_ROOT}/GLM_output${OUTPUT_SUFFIX}/"
