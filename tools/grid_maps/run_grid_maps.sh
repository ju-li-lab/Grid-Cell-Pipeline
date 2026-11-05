#!/usr/bin/env bash
# =============================================================================
#  Voxel-wise grid maps and orientation clusters, in parallel on the cluster
# =============================================================================
#  Turns what GLM1 already wrote into maps that say WHERE to look. Fits no
#  model, reads no functionals, and changes nothing in the GLM output.
#
#  Usage:
#    bash tools/grid_maps/run_grid_maps.sh --trial uw1_sm0
#    bash tools/grid_maps/run_grid_maps.sh --trial x --local         # here, serial
#    bash tools/grid_maps/run_grid_maps.sh --trial x --what clusters
#    bash tools/grid_maps/run_grid_maps.sh --trial x --sub sub-01s17
#
#  Submits a SLURM job array, ONE TASK PER SUBJECT, then a small merge job
#  that waits for the array and concatenates the per-subject CSVs. Pass
#  --local to skip SLURM and run the subjects one after another here.
#
#  Options:
#    --trial NAME     label written into the CSVs so trials can be compared.
#                     Defaults to the current date and time.
#    --glm DIR        GLM output directory holding the *_GLM1 folders.
#                     Defaults to OUTPUT_ROOT/GLM_output[_RUN_VARIANT].
#    --out DIR        where the maps and CSVs go. Default <glm>/grid_maps_out.
#    --masks DIR      ROI masks, as r<sub>_<ses>_hemi-*_mask.nii.
#                     Defaults to OUTPUT_ROOT/GLM_runauto/ROI.
#    --sessions LIST  comma-separated, DEFAULT ses-01 ONLY. Sessions 2 and 3
#                     are post-intervention, so a cross-session number there
#                     mixes measurement reliability with a real effect.
#    --all-sessions   every session instead.
#    --sub SUB        limit to one subject (repeatable).
#    --what WHICH     maps | clusters | both (default both).
#    --anat DIR       preprocessed derivatives tree, to also build the
#                     cross-session R map. Needs 2+ sessions, so it does
#                     nothing under the default --sessions ses-01.
#    --radius MM      sphere radius for coh-local/z-local, default 5.
#    --angle DEG      cluster agreement tolerance in orientation degrees,
#                     default 5. Read the note on percolation in
#                     grid_clusters.m before raising it.
#    --perms N        shuffles for the cluster null, default 500. 0 skips it.
#    --min-size N     smallest cluster written to the CSV, default 10.
#    --run N          cluster a single run instead of the across-run mean.
#    --local          run here, serially, no SLURM.
#    -h, --help       this text.
#
#  THIS WHOLE DIRECTORY IS OPTIONAL. Nothing in the pipeline calls it.
#  Delete tools/grid_maps/ and everything else keeps working.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
CONFIG_FILE="${ROOT_DIR}/config.sh"

TRIAL="$(date +%Y%m%d_%H%M%S)"
GLM_DIR=""; OUT_DIR=""; MASK_DIR=""; ANAT_ROOT=""
SESSIONS="ses-01"
WHAT="both"
RADIUS_MM=5; ANGLE_TOL=5; PERMS=500; MIN_SIZE=10; CLUST_RUN=0
LOCAL="false"
SUBS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --trial)        TRIAL="$2";     shift 2 ;;
        --glm)          GLM_DIR="$2";   shift 2 ;;
        --out)          OUT_DIR="$2";   shift 2 ;;
        --masks)        MASK_DIR="$2";  shift 2 ;;
        --anat)         ANAT_ROOT="$2"; shift 2 ;;
        --sessions)     SESSIONS="$2";  shift 2 ;;
        --all-sessions) SESSIONS="";    shift ;;
        --sub)          SUBS+=("$2");   shift 2 ;;
        --what)         WHAT="$2";      shift 2 ;;
        --radius)       RADIUS_MM="$2"; shift 2 ;;
        --angle)        ANGLE_TOL="$2"; shift 2 ;;
        --perms)        PERMS="$2";     shift 2 ;;
        --min-size)     MIN_SIZE="$2";  shift 2 ;;
        --run)          CLUST_RUN="$2"; shift 2 ;;
        --local)        LOCAL="true";   shift ;;
        -h|--help)
            sed -n '2,50p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
    esac
done

case "$WHAT" in
    maps|clusters|both) ;;
    *) echo "ERROR: --what must be maps, clusters or both (got: $WHAT)" >&2; exit 1 ;;
esac

[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config.sh not found: $CONFIG_FILE" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"
# shellcheck disable=SC1091
source "${ROOT_DIR}/bids_common.sh"

: "${MATLAB_MODULE:=matlab}"
: "${OUTPUT_ROOT:?ERROR: OUTPUT_ROOT not set in config.sh}"

if [[ -z "$GLM_DIR" ]]; then
    variant="$(printf '%s' "${RUN_VARIANT:-}" | tr -c 'A-Za-z0-9._-' '_' | sed -E 's/_+/_/g; s/^_//; s/_$//')"
    GLM_DIR="${OUTPUT_ROOT}/GLM_output${variant:+_${variant}}"
fi
[[ -z "$MASK_DIR" ]] && MASK_DIR="${OUTPUT_ROOT}/GLM_runauto/ROI"
[[ -z "$OUT_DIR"  ]] && OUT_DIR="${GLM_DIR}/grid_maps_out"

if [[ ! -d "$GLM_DIR" ]]; then
    echo "ERROR: no GLM output at: $GLM_DIR" >&2
    echo "       Run 03_gridcat.sh first, or pass --glm DIR." >&2
    exit 1
fi

# ---- subjects: one array task each ------------------------------------
if [[ ${#SUBS[@]} -eq 0 ]]; then
    # not mapfile: bash 3.2 (stock macOS) has no such builtin
    while IFS= read -r _sub; do
        [[ -n "$_sub" ]] && SUBS+=("$_sub")
    done < <(
        for d in "$GLM_DIR"/sub-*_GLM1; do
            [[ -d "$d" ]] || continue
            basename "$d" | sed -E 's/^(sub-[^_]+)_ses-.*$/\1/'
        done | sort -u
    )
fi
[[ ${#SUBS[@]} -gt 0 ]] || { echo "ERROR: no sub-*_GLM1 directories under $GLM_DIR" >&2; exit 1; }

mkdir -p "$OUT_DIR" "${OUT_DIR}/csv_parts"
SUBJECT_LIST="${OUT_DIR}/subjects_${TRIAL}.txt"
printf '%s\n' "${SUBS[@]}" > "$SUBJECT_LIST"

# MATLAB cell array of sessions
ses_arg="{}"
if [[ -n "$SESSIONS" ]]; then
    IFS=',' read -ra _s <<< "$SESSIONS"
    out="{"; for x in "${_s[@]}"; do out="${out}'${x}',"; done
    ses_arg="${out%,}}"
fi

echo "============================================"
echo "  Grid maps and orientation clusters"
echo "  Trial:    ${TRIAL}"
echo "  GLM1 out: ${GLM_DIR}"
echo "  Masks:    ${MASK_DIR}"
echo "  Sessions: ${SESSIONS:-<all>}"
echo "  Subjects: ${#SUBS[@]} (${SUBS[*]})"
echo "  What:     ${WHAT}"
echo "  Out:      ${OUT_DIR}"
echo "============================================"

export GLM_DIR OUT_DIR MASK_DIR ANAT_ROOT TRIAL WHAT \
       RADIUS_MM ANGLE_TOL PERMS MIN_SIZE CLUST_RUN SUBJECT_LIST
export SESSIONS_ARG="$ses_arg"
export PIPELINE_SCRIPT_DIR="$ROOT_DIR"
export GRID_MAPS_DIR="$SELF_DIR"

# ---- local, serial ----------------------------------------------------
if [[ "$LOCAL" == "true" ]]; then
    load_modules "${MATLAB_MODULE}"
    require_command matlab "running grid_maps"
    for i in "${!SUBS[@]}"; do
        SLURM_ARRAY_TASK_ID=$((i+1)) bash "${SELF_DIR}/grid_maps.sbatch"
    done
    echo; echo "merging:"
    bash "${SELF_DIR}/merge_csv_parts.sh" "$OUT_DIR"
    exit 0
fi

# ---- SLURM: array of subjects, then a dependent merge -----------------
require_command sbatch "submitting the job array"
mkdir -p "${ROOT_DIR}/logs"

cd "$ROOT_DIR"
ARRAY_ID="$(sbatch --parsable \
    --array="1-${#SUBS[@]}" \
    --partition="${SLURM_PARTITION:-compute}" \
    --export=ALL \
    "${SELF_DIR}/grid_maps.sbatch")"
echo "array job: ${ARRAY_ID} (${#SUBS[@]} tasks)"

MERGE_ID="$(sbatch --parsable \
    --job-name=gridmaps_merge \
    --partition="${SLURM_PARTITION:-compute}" \
    --dependency="afterany:${ARRAY_ID}" \
    --time=00:10:00 --mem=2G \
    --output="logs/%x_%j.out" --error="logs/%x_%j.err" \
    --export=ALL \
    --wrap="bash '${SELF_DIR}/merge_csv_parts.sh' '${OUT_DIR}'")"
echo "merge job: ${MERGE_ID} (waits for the array)"
echo
echo "watch:  squeue -j ${ARRAY_ID},${MERGE_ID}"
echo "then:   column -s, -t ${OUT_DIR}/clusters_summary.csv | less -S"
