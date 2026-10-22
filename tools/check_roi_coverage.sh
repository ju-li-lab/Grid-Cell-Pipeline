#!/usr/bin/env bash
# =============================================================================
#  How much of each ROI mask survives SPM's implicit mask
# =============================================================================
#  Run after preprocessing, before committing to a GridCAT run. Run it again
#  after changing a preprocessing setting, with a different --trial name, to
#  see whether the change helped.
#
#  Usage:
#    bash tools/check_roi_coverage.sh --trial sep4_fwhm5
#    bash tools/check_roi_coverage.sh --trial sep15_fwhm1 --deriv /other/spm-preproc
#    bash tools/check_roi_coverage.sh --trial x --sub sub-05s17
#
#  Options:
#    --trial NAME   label for this run, written into the CSV so trials can be
#                   compared. Defaults to the current date and time.
#    --deriv DIR    preprocessing derivatives to measure. Defaults to the one
#                   config.sh points at, so normally you can leave it out.
#    --csv FILE     CSV to append to (default: <pipeline root>/roi_coverage.csv).
#                   One row per subject x session x hemisphere, per trial.
#    --sub SUB      limit to one subject (repeatable).
#    -h, --help     this text.
#
#  Reading the output:
#    no_data   voxels with essentially no signal. These are dropout or a
#              misplaced mask; NO masking threshold recovers them.
#    n_<thr>   voxels that would carry a value in the GLM at that threshold.
#
#  Set C.GLM.maskingThreshold in spm_options.m to the lowest threshold that is
#  still on the flat part of the curve.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# config.sh and bids_common.sh live one level up, at the pipeline root.
ROOT_DIR="$(cd "${SELF_DIR}/.." && pwd)"
CONFIG_FILE="${ROOT_DIR}/config.sh"

TRIAL="$(date +%Y%m%d_%H%M%S)"
DERIV=""
CSV=""
SUBS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --trial) TRIAL="$2"; shift 2 ;;
        --deriv) DERIV="$2"; shift 2 ;;
        --csv)   CSV="$2";   shift 2 ;;
        --sub)   SUBS+=("$2"); shift 2 ;;
        -h|--help)
            sed -n '2,33p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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
[[ -z "$DERIV" ]] && DERIV="$(deriv_preproc)"
[[ -z "$CSV"   ]] && CSV="${ROOT_DIR}/roi_coverage.csv"

if [[ ! -d "$DERIV" ]]; then
    echo "ERROR: no preprocessing derivatives at: $DERIV" >&2
    echo "       Run the preprocessing first, or pass --deriv DIR." >&2
    exit 1
fi

# MATLAB cell array for the Subjects option, e.g. {'sub-01s17','sub-05s17'}
sub_arg="{}"
if [[ ${#SUBS[@]} -gt 0 ]]; then
    sub_arg="{"
    for s in "${SUBS[@]}"; do sub_arg="${sub_arg}'${s}',"; done
    sub_arg="${sub_arg%,}}"
fi

echo "============================================"
echo "  ROI coverage check"
echo "  Trial:       ${TRIAL}"
echo "  Derivatives: ${DERIV}"
echo "  CSV:         ${CSV}"
echo "============================================"

load_modules "${MATLAB_MODULE}"
require_command matlab "running check_roi_coverage"

matlab -nodisplay -nosplash -r "\
    addpath(genpath('${ROOT_DIR}/matlab')); \
    try, \
        check_roi_coverage('${DERIV}', \
            'Trial', '${TRIAL}', \
            'CSV', '${CSV}', \
            'Subjects', ${sub_arg}); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
        exit(1); \
    end; \
    exit(0);"
