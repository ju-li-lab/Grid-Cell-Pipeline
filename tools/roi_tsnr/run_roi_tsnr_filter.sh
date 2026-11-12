#!/usr/bin/env bash
# =============================================================================
#  Filter ROI voxels by tSNR - one job per subject x session
# =============================================================================
#  Rule: within each subject x session, drop every ROI voxel with
#  tSNR < Q1 - 1.5*IQR of the pooled (left + right) ROI tSNR distribution.
#  Runs on all sessions by default.
#
#    bash tools/roi_tsnr/run_roi_tsnr_filter.sh --dry-run --local   # look first
#    bash tools/roi_tsnr/run_roi_tsnr_filter.sh                     # write masks
#
#  Output defaults to <rois>_tsnr (GLM_runauto/ROI_tsnr), which is what
#  03_gridcat.sh reads when ROI_TSNR_FILTER=1 in config.sh.
#
#  Submits a SLURM array (one task per subject x session) plus a merge job.
#  --local runs them serially instead.
#
#  Options:
#    --dry-run        report only, write nothing.
#    --deriv DIR      derivatives tree with the functionals.
#                     Default OUTPUT_ROOT/derivatives/<DERIV_PREPROC>.
#    --rois DIR       flat dir of r<sub>_<ses>_hemi-*_label-*_mask.nii.
#                     Default OUTPUT_ROOT/GLM_runauto/ROI.
#    --out DIR        default <rois>_tsnr.
#    --prefix P       functional prefix, default 'u' (unwarped, UNSMOOTHED).
#                     Use 'ua' where slice timing ran first. Never 's*'.
#    --combine WHICH  mean (default) or min across runs.
#    --sessions LIST  comma-separated, default all sessions.
#    --sub SUB        limit to one subject (repeatable).
#    --tasks LIST     comma-separated task numbers, default 1,2,3.
#    --hpf SEC        high-pass, default 128, matching GLM.HPF_cutoff.
#    --min-voxels N   refuse to write a mask below this, default 5.
#    --local          run here, serially, no SLURM.
#    -h, --help       this text.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
CONFIG_FILE="${ROOT_DIR}/config.sh"

TRIAL="$(date +%Y%m%d_%H%M%S)"
DERIV_IN=""; ROIS=""; OUT=""
PREFIX="u"; COMBINE="mean"
SESSIONS=""; TASKS_IN="1,2,3"; HPF=128; MINVOX=5
DRY="false"; LOCAL="false"
SUBS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)      DRY="true";     shift ;;
        --deriv)        DERIV_IN="$2";  shift 2 ;;
        --rois)         ROIS="$2";      shift 2 ;;
        --out)          OUT="$2";       shift 2 ;;
        --prefix)       PREFIX="$2";    shift 2 ;;
        --combine)      COMBINE="$2";   shift 2 ;;
        --sessions)     SESSIONS="$2";  shift 2 ;;
        --sub)          SUBS+=("$2");   shift 2 ;;
        --tasks)        TASKS_IN="$2";  shift 2 ;;
        --hpf)          HPF="$2";       shift 2 ;;
        --min-voxels)   MINVOX="$2";    shift 2 ;;
        --local)        LOCAL="true";   shift ;;
        -h|--help)
            sed -n '2,36p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
    esac
done

case "$COMBINE" in
    mean|min) ;;
    *) echo "ERROR: --combine must be mean or min (got: $COMBINE)" >&2; exit 1 ;;
esac
case "$PREFIX" in
    s*) echo "ERROR: --prefix '$PREFIX' is a SMOOTHED series. tSNR on smoothed" >&2
        echo "       data measures distance to the tissue edge, not the voxel." >&2
        echo "       Use 'u' or 'ua'." >&2
        exit 1 ;;
esac

[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config.sh not found: $CONFIG_FILE" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"
# shellcheck disable=SC1091
source "${ROOT_DIR}/bids_common.sh"

: "${MATLAB_MODULE:=matlab}"
: "${OUTPUT_ROOT:?ERROR: OUTPUT_ROOT not set in config.sh}"

[[ -z "$DERIV_IN" ]] && DERIV_IN="${OUTPUT_ROOT}/derivatives/${DERIV_PREPROC:-spm-preproc}"
[[ -z "$ROIS"     ]] && ROIS="${OUTPUT_ROOT}/GLM_runauto/ROI"
[[ -z "$OUT"      ]] && OUT="${ROIS}_tsnr"
RTF_TASKS="[${TASKS_IN//,/ }]"

[[ -d "$DERIV_IN" ]] || { echo "ERROR: no derivatives tree at: $DERIV_IN" >&2; exit 1; }
[[ -d "$ROIS"     ]] || { echo "ERROR: no ROI directory at: $ROIS" >&2; exit 1; }

# ---- work list: one line per subject x session -------------------------
WORK=()
shopt -s nullglob
for f in "$ROIS"/r*_hemi-left_label-*_mask.nii; do
    b="$(basename "$f")"
    sub="$(printf '%s' "$b" | sed -E 's/^r(sub-[^_]+)_(ses-[^_]+)_hemi-left.*/\1/')"
    ses="$(printf '%s' "$b" | sed -E 's/^r(sub-[^_]+)_(ses-[^_]+)_hemi-left.*/\2/')"
    [[ "$sub" == sub-* && "$ses" == ses-* ]] || continue
    if [[ ${#SUBS[@]} -gt 0 ]]; then
        printf '%s\n' "${SUBS[@]}" | grep -qxF "$sub" || continue
    fi
    if [[ -n "$SESSIONS" ]]; then
        printf '%s' "$SESSIONS" | tr ',' '\n' | grep -qxF "$ses" || continue
    fi
    WORK+=("$sub $ses")
done
[[ ${#WORK[@]} -gt 0 ]] || { echo "ERROR: no masks under $ROIS matched the selection" >&2; exit 1; }

mkdir -p "$OUT" "${OUT}/csv_parts"
WORK_LIST="${OUT}/worklist_${TRIAL}.txt"
printf '%s\n' "${WORK[@]}" | sort > "$WORK_LIST"

echo "============================================"
echo "  ROI tSNR filter"
echo "  Trial:       ${TRIAL}"
echo "  Functionals: ${DERIV_IN}  (${PREFIX}sub-*, unsmoothed)"
echo "  ROI masks:   ${ROIS}"
echo "  Rule:        drop voxels with tSNR < Q1 - 1.5*IQR (per subject x session)"
echo "  Runs:        tasks ${TASKS_IN}, combined by ${COMBINE}"
echo "  Jobs:        ${#WORK[@]} (subject x session)"
[[ "$DRY" == "true" ]] && echo "  DRY RUN - nothing will be written"
echo "  Out:         ${OUT}"
echo "============================================"

export DERIV_ROOT_IN="$DERIV_IN" ROI_DIR="$ROIS" OUT_DIR="$OUT" WORK_LIST
export RTF_TASKS TRIAL
export FUNC_PREFIX="$PREFIX" COMBINE HPF_SEC="$HPF" MIN_VOXELS="$MINVOX" DRY_RUN="$DRY"
export PIPELINE_SCRIPT_DIR="$ROOT_DIR" ROI_TSNR_DIR="$SELF_DIR"

# ---- local, serial (also: inside an interactive srun) ------------------
if [[ "$LOCAL" == "true" ]]; then
    load_modules "${MATLAB_MODULE}"
    require_command matlab "running roi_tsnr_filter"
    for i in $(seq 1 ${#WORK[@]}); do
        SLURM_ARRAY_TASK_ID="$i" bash "${SELF_DIR}/roi_tsnr_filter.sbatch"
    done
    if [[ "$DRY" != "true" ]]; then
        echo; echo "merging:"
        bash "${SELF_DIR}/merge_csv_parts.sh" "$OUT"
    fi
    exit 0
fi

# ---- SLURM: array over subject x session, then a dependent merge -------
require_command sbatch "submitting the job array"
mkdir -p "${ROOT_DIR}/logs"

cd "$ROOT_DIR"
ARRAY_ID="$(sbatch --parsable \
    --array="1-${#WORK[@]}" \
    --partition="${SLURM_PARTITION:-compute}" \
    --export=ALL \
    "${SELF_DIR}/roi_tsnr_filter.sbatch")"
echo "array job: ${ARRAY_ID} (${#WORK[@]} tasks)"

if [[ "$DRY" != "true" ]]; then
    MERGE_ID="$(sbatch --parsable \
        --job-name=roitsnr_merge \
        --partition="${SLURM_PARTITION:-compute}" \
        --dependency="afterany:${ARRAY_ID}" \
        --time=00:10:00 --mem=2G \
        --output="logs/%x_%j.out" --error="logs/%x_%j.err" \
        --export=ALL \
        --wrap="bash '${SELF_DIR}/merge_csv_parts.sh' '${OUT}'")"
    echo "merge job: ${MERGE_ID} (waits for the array)"
fi
echo
echo "watch:  squeue -j ${ARRAY_ID}"
echo "then:   column -s, -t ${OUT}/roi_tsnr.csv | less -S"
