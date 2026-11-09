#!/usr/bin/env bash
# =============================================================================
#  Does event time 0 land on the first saved volume?
# =============================================================================
#  Slides the rotate/stand/translate design against the functional data and
#  reports the shift that fits best. Deliberately independent of the grid
#  analysis: it uses the task's visual structure, which is there whether or
#  not the entorhinal cortex shows a grid code.
#
#  The answer is delta_hat, in seconds. Correct onsets as  onset - delta_hat.
#    ~0        event time 0 is the first saved volume. Nothing to do.
#    ~+2.06    one dummy volume was discarded after the trigger fired.
#    ~+4.12    two were.
#
#  Usage:
#    bash tools/onset_check/run_onset_check.sh
#    bash tools/onset_check/run_onset_check.sh --sub sub-01s17 --sessions ses-01
#    bash tools/onset_check/run_onset_check.sh --per-slice --json /path/to_bold.json
#
#  Options:
#    --events DIR    event tables. Default: searched for under
#                    OUTPUT_ROOT/GLM_runauto (EventFiles/Event_tables, EventData,
#                    then anywhere below it)
#    --deriv DIR     preprocessed tree. Default from config.sh
#    --prefix P      series to read. Default from func_prefix() in bids_common.sh
#    --sub SUB       limit to one subject (repeatable)
#    --sessions LIST comma-separated, default all
#    --tasks LIST    comma-separated, default 1,2,3
#    --t0 N          microtime onset. Default 9 (exactly 0.5 x TR at 16 bins)
#    --tr SEC        repetition time. Default read from spm_options.m
#    --range N       sweep +/- N TR. Default 5
#    --keep PCT      share of voxels kept, by residual SD. Default 2
#    --per-slice     fit one curve per slice on a series with NO slice timing
#                    correction, to confirm the offset without assuming the
#                    SliceTiming vector. Needs --json.
#    --json FILE     a BOLD JSON sidecar, for --per-slice
#    --out DIR       output directory. Default OUTPUT_ROOT/onset_check
#    -h, --help      this text
#
#  THIS WHOLE DIRECTORY IS OPTIONAL. Nothing in the pipeline calls it.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
CONFIG_FILE="${ROOT_DIR}/config.sh"

EVENTS=""; DERIV=""; PREFIX=""; SESSIONS=""; TASKS="1,2,3"
T0=9; RANGE=5; KEEP=2; PERSLICE="false"; JSON=""; OUT=""; TR_VAL=""
SUBS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --events)    EVENTS="$2";   shift 2 ;;
        --deriv)     DERIV="$2";    shift 2 ;;
        --prefix)    PREFIX="$2";   shift 2 ;;
        --sub)       SUBS+=("$2");  shift 2 ;;
        --sessions)  SESSIONS="$2"; shift 2 ;;
        --tasks)     TASKS="$2";    shift 2 ;;
        --t0)        T0="$2";       shift 2 ;;
        --tr)        TR_VAL="$2";   shift 2 ;;
        --range)     RANGE="$2";    shift 2 ;;
        --keep)      KEEP="$2";     shift 2 ;;
        --per-slice) PERSLICE="true"; shift ;;
        --json)      JSON="$2";     shift 2 ;;
        --out)       OUT="$2";      shift 2 ;;
        -h|--help)
            sed -n '2,46p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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
: "${SPM_DIR:?ERROR: SPM_DIR not set in config.sh}"

[[ -z "$DERIV"  ]] && DERIV="$(deriv_preproc)"
# Event tables: the pipeline itself tries a known layout and then falls back to
# a recursive scan (run_gridcat_analysis.m), so do the same rather than
# hard-coding one directory.
if [[ -z "$EVENTS" ]]; then
    GLM_RUNAUTO="${OUTPUT_ROOT}/GLM_runauto"
    for cand in "${GLM_RUNAUTO}/EventFiles/Event_tables" \
                "${GLM_RUNAUTO}/EventData" \
                "${GLM_RUNAUTO}"; do
        if compgen -G "${cand}/*EventData.txt" > /dev/null 2>&1; then
            EVENTS="$cand"; break
        fi
    done
    if [[ -z "$EVENTS" && -d "$GLM_RUNAUTO" ]]; then
        found="$(find "$GLM_RUNAUTO" -name '*EventData.txt' ! -name '._*' 2>/dev/null | head -1 || true)"
        [[ -n "$found" ]] && EVENTS="$(dirname "$found")"
    fi
fi
[[ -z "$OUT"    ]] && OUT="${OUTPUT_ROOT}/onset_check"

# TR lives in spm_options.m, not config.sh.
if [[ -z "$TR_VAL" ]]; then
    TR_VAL="$(grep -oE 'C\.TR[[:space:]]*=[[:space:]]*[0-9.]+' "${ROOT_DIR}/spm_options.m" \
              | head -1 | grep -oE '[0-9.]+$')"
fi
[[ -z "$TR_VAL" ]] && { echo "ERROR: could not read C.TR from spm_options.m; pass --tr." >&2; exit 1; }

# --per-slice needs the series BEFORE slice timing, which is the staged raw.
if [[ -z "$PREFIX" ]]; then
    if [[ "$PERSLICE" == "true" ]]; then PREFIX=""; else PREFIX="$(func_prefix)"; fi
fi

if [[ "$PERSLICE" == "true" && -z "$JSON" ]]; then
    echo "ERROR: --per-slice needs --json pointing at a BOLD JSON sidecar." >&2
    exit 1
fi
[[ ! -d "$DERIV"  ]] && { echo "ERROR: no preprocessed tree at: $DERIV" >&2; exit 1; }
if [[ -z "$EVENTS" || ! -d "$EVENTS" ]]; then
    echo "ERROR: could not find a directory holding *EventData.txt." >&2
    echo "       Looked under ${OUTPUT_ROOT}/GLM_runauto. Pass --events DIR." >&2
    exit 1
fi
if ! compgen -G "${EVENTS}/*.txt" > /dev/null 2>&1; then
    echo "ERROR: ${EVENTS} holds no .txt event tables." >&2
    exit 1
fi


to_cell() { local out="{"; for x in "$@"; do out="${out}'${x}',"; done; printf '%s}' "${out%,}"; }
sub_arg="{}"; [[ ${#SUBS[@]} -gt 0 ]] && sub_arg="$(to_cell "${SUBS[@]}")"
ses_arg="{}"
if [[ -n "$SESSIONS" ]]; then IFS=',' read -ra _s <<< "$SESSIONS"; ses_arg="$(to_cell "${_s[@]}")"; fi

echo "============================================"
echo "  onset_check"
echo "  Tree:     ${DERIV}"
echo "  Events:   ${EVENTS}"
echo "            $(find "$EVENTS" -maxdepth 1 -name '*.txt' ! -name '._*' 2>/dev/null | wc -l | tr -d ' ') event tables"
echo "  Prefix:   \"${PREFIX}\"$([[ "$PERSLICE" == "true" ]] && echo '  (no slice timing)')"
echo "  TR:       ${TR_VAL}"
echo "  Tasks:    ${TASKS}"
echo "  Sessions: ${SESSIONS:-<all>}"
echo "  Out:      ${OUT}"
echo "============================================"

load_modules "${MATLAB_MODULE}"
require_command matlab "running onset_check"

st_arg="[]"; [[ -n "$JSON" ]] && st_arg="'${JSON}'"

matlab -nodisplay -nosplash -r "\
    addpath('${SPM_DIR}'); addpath('${SELF_DIR}'); addpath('${ROOT_DIR}'); \
    try, \
        onset_check('${DERIV}', '${EVENTS}', \
            'Subjects', ${sub_arg}, \
            'Sessions', ${ses_arg}, \
            'Tasks', [${TASKS//,/ }], \
            'Prefix', '${PREFIX}', \
            'RawPrefix', '${PREFIX}', \
            'TR', ${TR_VAL}, \
            'MicrotimeOnset', ${T0}, \
            'RangeTR', ${RANGE}, \
            'KeepPct', ${KEEP}, \
            'PerSlice', ${PERSLICE}, \
            'SliceTiming', ${st_arg}, \
            'OutDir', '${OUT}'); \
    catch ME, \
        fprintf(2, '\n=== ERROR ===\n%s\n', getReport(ME, 'extended')); \
        exit(1); \
    end; \
    exit(0);"
