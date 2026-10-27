#!/usr/bin/env bash
# =============================================================================
#  Per-run data quality
# =============================================================================
#  Which runs are worth keeping? Reports temporal SNR (whole brain and inside
#  the ROI), head motion, signal spikes, and the fit of the task model - then
#  measures the SPLIT-HALF RELIABILITY of the task fit on your own data, so
#  you can see which numbers you are entitled to act on.
#
#  Usage:
#    bash tools/run_quality/run_check_run_quality.sh --trial sep4_fwhm5
#    bash tools/run_quality/run_check_run_quality.sh --trial x --sessions ses-01
#
#  Options:
#    --trial NAME     label written into the CSV. Default: date and time.
#    --deriv DIR      preprocessing derivatives. Default: from config.sh.
#    --glm DIR        GLM output with the *_GLM1 folders, needed for task_r2.
#                     Default: OUTPUT_ROOT/GLM_output[_RUN_VARIANT].
#    --events DIR     EventData. Default: OUTPUT_ROOT/GLM_runauto/EventData.
#    --roi DIR        ROI masks.  Default: OUTPUT_ROOT/GLM_runauto/ROI.
#    --sessions LIST  comma-separated, default all sessions.
#    --sub SUB        limit to one subject (repeatable).
#    --csv FILE       CSV to append to. Default <pipeline root>/run_quality.csv
#    -h, --help       this text.
#
#  READ THE RELIABILITY LINE BEFORE ACTING. On the dataset this was built
#  against, task_r2 had a split-half reliability of 0.10 - the two halves of
#  one run disagreed almost completely - so it could not be used to drop
#  runs. tsnr_roi and mean_fd were fine. Your data may differ, which is why
#  the tool measures it every time instead of assuming.
#
#  THIS WHOLE DIRECTORY IS OPTIONAL. Nothing in the pipeline calls it.
#  Delete tools/run_quality/ and everything else keeps working.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
CONFIG_FILE="${ROOT_DIR}/config.sh"

TRIAL="$(date +%Y%m%d_%H%M%S)"
DERIV=""; GLM_DIR=""; EV_DIR=""; ROI_DIR=""; SESSIONS=""; CSV=""; SUBS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --trial)    TRIAL="$2";    shift 2 ;;
        --deriv)    DERIV="$2";    shift 2 ;;
        --glm)      GLM_DIR="$2";  shift 2 ;;
        --events)   EV_DIR="$2";   shift 2 ;;
        --roi)      ROI_DIR="$2";  shift 2 ;;
        --sessions) SESSIONS="$2"; shift 2 ;;
        --sub)      SUBS+=("$2");  shift 2 ;;
        --csv)      CSV="$2";      shift 2 ;;
        -h|--help)  sed -n '2,37p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
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

[[ -z "$DERIV" ]] && DERIV="$(deriv_preproc)"
if [[ -z "$GLM_DIR" ]]; then
    variant="$(printf '%s' "${RUN_VARIANT:-}" | tr -c 'A-Za-z0-9._-' '_' | sed -E 's/_+/_/g; s/^_//; s/_$//')"
    GLM_DIR="${OUTPUT_ROOT}/GLM_output${variant:+_${variant}}"
fi
[[ -z "$EV_DIR"  ]] && EV_DIR="${OUTPUT_ROOT}/GLM_runauto/EventData"
[[ -z "$ROI_DIR" ]] && ROI_DIR="${OUTPUT_ROOT}/GLM_runauto/ROI"
[[ -z "$CSV"     ]] && CSV="${ROOT_DIR}/run_quality.csv"

[[ -d "$DERIV" ]] || { echo "ERROR: no derivatives at: $DERIV" >&2; exit 1; }
[[ -d "$GLM_DIR" ]] || { echo "NOTE: no GLM output at $GLM_DIR - task_r2 will be skipped."; GLM_DIR=""; }

to_cell() { local out="{"; for x in "$@"; do out="${out}'${x}',"; done; printf '%s}' "${out%,}"; }
sub_arg="{}"; [[ ${#SUBS[@]} -gt 0 ]] && sub_arg="$(to_cell "${SUBS[@]}")"
ses_arg="{}"
if [[ -n "$SESSIONS" ]]; then IFS=',' read -ra _s <<< "$SESSIONS"; ses_arg="$(to_cell "${_s[@]}")"; fi

echo "============================================"
echo "  Per-run data quality"
echo "  Trial:       ${TRIAL}"
echo "  Derivatives: ${DERIV}"
echo "  GLM1:        ${GLM_DIR:-<none, task_r2 skipped>}"
echo "  Sessions:    ${SESSIONS:-<all>}"
echo "  CSV:         ${CSV}"
echo "============================================"

load_modules "${MATLAB_MODULE}"
require_command matlab "running check_run_quality"

matlab -nodisplay -nosplash -r "\
    addpath('${SELF_DIR}'); \
    try, \
        check_run_quality('${DERIV}', \
            'GLMRoot','${GLM_DIR}', 'EventDir','${EV_DIR}', 'ROIDir','${ROI_DIR}', \
            'Sessions', ${ses_arg}, 'Subjects', ${sub_arg}, \
            'Trial','${TRIAL}', 'CSV','${CSV}'); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME,'extended')); exit(1); \
    end; exit(0);"
