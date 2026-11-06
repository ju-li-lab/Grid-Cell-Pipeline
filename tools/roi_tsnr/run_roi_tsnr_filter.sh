#!/usr/bin/env bash
# =============================================================================
#  Filter ROI voxels by tSNR - one job per subject x session
# =============================================================================
#  Look before you cut:
#    bash tools/roi_tsnr/run_roi_tsnr_filter.sh --dry-run --local
#
#  Then write the masks (left and right can differ - dropout is asymmetric):
#    bash tools/roi_tsnr/run_roi_tsnr_filter.sh --keep 60,40
#    bash tools/roi_tsnr/run_roi_tsnr_filter.sh --min-tsnr 45,55
#
#  Submits a SLURM array, ONE TASK PER SUBJECT x SESSION, then a merge job.
#  --local runs them one after another instead, which is what you want in an
#  interactive session:  srun --pty -c 2 --mem 8G bash   then --local.
#
#  Options:
#    --keep L,R       percent of each hemisphere to keep, by tSNR rank.
#                     Default 50,50. One value sets both. RANK RULE: n stays
#                     constant across subjects, which is what a group test
#                     wants.
#    --min-tsnr L,R   absolute tSNR cut instead. Overrides --keep. Lets n vary
#                     between subjects, which weights them unequally.
#    --dry-run        print the sweep table (100/75/50/25%) and write nothing.
#    --deriv DIR      derivatives tree with the functionals.
#                     Default OUTPUT_ROOT/derivatives/<DERIV_PREPROC>.
#    --rois DIR       flat dir of r<sub>_<ses>_hemi-*_label-*_mask.nii.
#                     Default OUTPUT_ROOT/GLM_runauto/ROI.
#    --out DIR        default <rois>_tsnr_L<pct>_R<pct>, so the folder name is
#                     the record. Filenames are unchanged, so config.sh needs
#                     no edit - just point the ROI dir at the new folder.
#    --prefix P       functional prefix, default 'u' (unwarped, UNSMOOTHED).
#                     Use 'ua' where slice timing ran first. NEVER 'su'.
#    --combine WHICH  mean (default) or min across the three runs. min is
#                     stricter and arguably right: SPM's implicit mask is an
#                     intersection across runs, so the worst run already
#                     decides whether a voxel survives.
#    --sessions LIST  comma-separated, DEFAULT ses-01. 2 and 3 are post-
#                     intervention: fix the voxel set on ses-01 and carry it
#                     forward rather than re-selecting, or a tSNR shift between
#                     sessions reads as an intervention effect.
#    --all-sessions   every session instead.
#    --sub SUB        limit to one subject (repeatable).
#    --tasks LIST     comma-separated task numbers, default 1,2,3.
#    --hpf SEC        high-pass, default 128, matching GLM.HPF_cutoff.
#    --min-voxels N   refuse to write a mask below this, default 5.
#    --local          run here, serially, no SLURM.
#    -h, --help       this text.
#
#  THIS WHOLE DIRECTORY IS OPTIONAL. Nothing in the pipeline calls it.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
CONFIG_FILE="${ROOT_DIR}/config.sh"

TRIAL="$(date +%Y%m%d_%H%M%S)"
DERIV_IN=""; ROIS=""; OUT=""
KEEP="50,50"; MINT=""; PREFIX="u"; COMBINE="mean"
SESSIONS="ses-01"; TASKS_IN="1,2,3"; HPF=128; MINVOX=5
DRY="false"; LOCAL="false"
SUBS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep)         KEEP="$2";      shift 2 ;;
        --min-tsnr)     MINT="$2";      shift 2 ;;
        --dry-run)      DRY="true";     shift ;;
        --deriv)        DERIV_IN="$2";  shift 2 ;;
        --rois)         ROIS="$2";      shift 2 ;;
        --out)          OUT="$2";       shift 2 ;;
        --prefix)       PREFIX="$2";    shift 2 ;;
        --combine)      COMBINE="$2";   shift 2 ;;
        --sessions)     SESSIONS="$2";  shift 2 ;;
        --all-sessions) SESSIONS="";    shift ;;
        --sub)          SUBS+=("$2");   shift 2 ;;
        --tasks)        TASKS_IN="$2";  shift 2 ;;
        --hpf)          HPF="$2";       shift 2 ;;
        --min-voxels)   MINVOX="$2";    shift 2 ;;
        --local)        LOCAL="true";   shift ;;
        -h|--help)
            sed -n '2,50p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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
        echo "       data is a distance-to-edge map, not a voxel property:" >&2
        echo "       5 mm FWHM on these voxels multiplies tSNR by 5.4, which is" >&2
        echo "       an average over ~30 voxels, wider than the ROI itself." >&2
        echo "       Use 'u' or 'ua'. See the README if you really mean it." >&2
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

# two numbers, or one meaning both
to_pair() {
    local v="$1" a b
    a="${v%%,*}"; b="${v##*,}"
    printf '[%s %s]' "$a" "$b"
}
KEEP_PCT="$(to_pair "$KEEP")"
MIN_TSNR="[]"
if [[ -n "$MINT" ]]; then
    MIN_TSNR="$(to_pair "$MINT")"
    [[ -z "$OUT" ]] && OUT="${ROIS}_tsnr_Labs${MINT%%,*}_Rabs${MINT##*,}"
fi
[[ -z "$OUT" ]] && OUT="${ROIS}_tsnr_L${KEEP%%,*}_R${KEEP##*,}"
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
if [[ -n "$MINT" ]]; then
    echo "  Rule:        tSNR >= ${MINT%%,*} (left), >= ${MINT##*,} (right)"
else
    echo "  Rule:        keep top ${KEEP%%,*}% (left), ${KEEP##*,}% (right)"
fi
echo "  Runs:        tasks ${TASKS_IN}, combined by ${COMBINE}"
echo "  Jobs:        ${#WORK[@]} (subject x session)"
[[ "$DRY" == "true" ]] && echo "  DRY RUN - nothing will be written"
echo "  Out:         ${OUT}"
echo "============================================"

export DERIV_ROOT_IN="$DERIV_IN" ROI_DIR="$ROIS" OUT_DIR="$OUT" WORK_LIST
export KEEP_PCT MIN_TSNR RTF_TASKS TRIAL
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
echo "then:   column -s, -t ${OUT}/roi_tsnr_sweep.csv | less -S"
