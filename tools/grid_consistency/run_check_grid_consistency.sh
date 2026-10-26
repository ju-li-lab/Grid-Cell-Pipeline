#!/usr/bin/env bash
# =============================================================================
#  Cross-run grid orientation consistency — preprocessing QC
# =============================================================================
#  Asks one question: does a voxel report the same grid orientation in run 1
#  as in run 2? That is a preprocessing question, separate from whether
#  aligned beats misaligned, which is an analysis question.
#
#  Reads only the per-run images GLM1 writes, never allRunsAvg, so nothing
#  here depends on GLM2_averageMeanGridOriAcrossRuns. Change that setting
#  freely; these numbers stay comparable.
#
#  Usage:
#    bash tools/grid_consistency/run_check_grid_consistency.sh --trial sep4_fwhm5
#    bash tools/grid_consistency/run_check_grid_consistency.sh --trial x --sub sub-05s17
#    bash tools/grid_consistency/run_check_grid_consistency.sh --trial x --all-sessions
#
#  Options:
#    --trial NAME     label written into the CSV so trials can be compared.
#                     Defaults to the current date and time.
#    --glm DIR        GLM output directory holding the *_GLM1 folders.
#                     Defaults to OUTPUT_ROOT/GLM_output[_RUN_VARIANT].
#    --masks DIR      ROI masks. Defaults to OUTPUT_ROOT/GLM_runauto/ROI.
#    --sessions LIST  comma-separated, default ses-01 only. Sessions after an
#                     intervention can differ for real reasons, which would
#                     confound a preprocessing comparison.
#    --all-sessions   measure every session instead.
#    --sub SUB        limit to one subject (repeatable).
#    --csv FILE       CSV to append to (default <pipeline root>/grid_consistency.csv)
#    --no-shift-test  skip the mask-displacement scan.
#    -h, --help       this text.
#
#  THIS WHOLE DIRECTORY IS OPTIONAL. Nothing in the pipeline calls it.
#  Delete tools/grid_consistency/ and everything else keeps working.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# this tool sits two levels down, not one
ROOT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
CONFIG_FILE="${ROOT_DIR}/config.sh"

TRIAL="$(date +%Y%m%d_%H%M%S)"
GLM_DIR=""
MASK_DIR=""
SESSIONS="ses-01"
CSV=""
SHIFT="true"
SUBS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --trial)        TRIAL="$2";    shift 2 ;;
        --glm)          GLM_DIR="$2";  shift 2 ;;
        --masks)        MASK_DIR="$2"; shift 2 ;;
        --sessions)     SESSIONS="$2"; shift 2 ;;
        --all-sessions) SESSIONS="";   shift ;;
        --sub)          SUBS+=("$2");  shift 2 ;;
        --csv)          CSV="$2";      shift 2 ;;
        --no-shift-test) SHIFT="false"; shift ;;
        -h|--help)
            sed -n '2,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
    esac
done

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
[[ -z "$CSV"      ]] && CSV="${ROOT_DIR}/grid_consistency.csv"

if [[ ! -d "$GLM_DIR" ]]; then
    echo "ERROR: no GLM output at: $GLM_DIR" >&2
    echo "       Run 03_gridcat.sh first, or pass --glm DIR." >&2
    exit 1
fi

# MATLAB cell arrays
to_cell() { local out="{"; for x in "$@"; do out="${out}'${x}',"; done; printf '%s}' "${out%,}"; }
sub_arg="{}"; [[ ${#SUBS[@]} -gt 0 ]] && sub_arg="$(to_cell "${SUBS[@]}")"
ses_arg="{}"
if [[ -n "$SESSIONS" ]]; then
    IFS=',' read -ra _s <<< "$SESSIONS"
    ses_arg="$(to_cell "${_s[@]}")"
fi

echo "============================================"
echo "  Grid orientation consistency"
echo "  Trial:    ${TRIAL}"
echo "  GLM1 out: ${GLM_DIR}"
echo "  Masks:    ${MASK_DIR}"
echo "  Sessions: ${SESSIONS:-<all>}"
echo "  CSV:      ${CSV}"
echo "============================================"

load_modules "${MATLAB_MODULE}"
require_command matlab "running check_grid_consistency"

matlab -nodisplay -nosplash -r "\
    addpath('${SELF_DIR}'); \
    try, \
        check_grid_consistency('${GLM_DIR}', \
            'Trial', '${TRIAL}', \
            'MaskDir', '${MASK_DIR}', \
            'Sessions', ${ses_arg}, \
            'Subjects', ${sub_arg}, \
            'ShiftTest', ${SHIFT}, \
            'CSV', '${CSV}'); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
        exit(1); \
    end; \
    exit(0);"
