#!/bin/bash
# ======================================================================
# SCAN_MULTIRUN.SH
# ======================================================================
# Finds every subject/session that holds more than one run of the same
# scan — T1w, T2w, a task BOLD, a phasediff, a magnitude, a ready-made
# fieldmap, or the reverse phase-encode EPI — and writes run_selection.tsv
# so you can say which run the pipeline should use.
#
# It also reads AcquisitionTime out of the JSON sidecars and groups the
# session's scans into blocks: a subject who climbs out of the scanner and
# comes back leaves a gap of many minutes, and every run counter restarts
# independently on the way back in. That is why run-2 of a T2w cannot be
# assumed to belong with run-2 of a task, and why the suggestion in the file
# is made per block rather than per run number.
#
# Usage:
#   bash scan_multirun.sh                      # scan everything in BIDS_ROOT
#   bash scan_multirun.sh --list subses_list.txt
#   bash scan_multirun.sh --fresh              # discard previous selections
#   bash scan_multirun.sh path/to/pipeline_config.cfg
#
# Options:
#   --list FILE     Only scan the subject-sessions in this list.
#   --output FILE   Where to write (default: RUN_SELECTION_FILE, or
#                   <script dir>/run_selection.tsv).
#   --gap-minutes N Scans further apart than this are different blocks
#                   (default: RUN_MATCH_GAP_MIN from the config).
#   --fresh         Do not carry over the selections already in the file.
#   -h, --help      Show this help.
#
# Output:
#   run_selection.tsv, with columns
#     subject  session  type  available_runs  selected_run
#     functional_block  runs_detail  notes
#
#   selected_run is pre-filled with a suggestion. Check it: runs_detail shows
#   every run as run-N@HH:MM:SS[block]#series, so a suggestion from the wrong
#   block is visible at a glance.
# ======================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG_FILE=""
LIST=""
OUTPUT=""
GAP_MIN=""
KEEP_EXISTING=1

show_help() {
    awk 'NR>1 { if ($0 !~ /^#/) exit; print }' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --list)        LIST="$2";     shift 2 ;;
        --output)      OUTPUT="$2";   shift 2 ;;
        --gap-minutes) GAP_MIN="$2";  shift 2 ;;
        --fresh)       KEEP_EXISTING=0; shift ;;
        -h|--help)     show_help; exit 0 ;;
        -*)
            echo "ERROR: Unknown option: $1" >&2
            echo "       Run with --help to see the accepted options." >&2
            exit 1 ;;
        *)
            CONFIG_FILE="$1"; shift ;;
    esac
done

[[ -z "$CONFIG_FILE" ]] && CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/bids_common.sh"

if [[ -z "${BIDS_ROOT:-}" ]] || [[ ! -d "$BIDS_ROOT" ]]; then
    echo "ERROR: BIDS_ROOT not set or does not exist: ${BIDS_ROOT:-<not set>}" >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 not found — needed to read the JSON sidecars." >&2
    exit 1
fi

[[ -z "$OUTPUT" ]]  && OUTPUT="$(run_selection_file)"
[[ -z "$GAP_MIN" ]] && GAP_MIN="${RUN_MATCH_GAP_MIN:-20}"

BLUE='\033[0;34m'; NC='\033[0m'

printf "${BLUE}========================================${NC}\n"
printf "${BLUE}  Multi-Run Scanner${NC}\n"
printf "${BLUE}========================================${NC}\n"
printf "  BIDS root:  %s\n" "$BIDS_ROOT"
printf "  Output:     %s\n" "$OUTPUT"
printf "  Block gap:  %s minutes\n" "$GAP_MIN"
[[ -n "$LIST" ]] && printf "  Restricted: %s\n" "$LIST"
printf "\n"

args=(
    --bids-root "$BIDS_ROOT"
    --output "$OUTPUT"
    --gap-minutes "$GAP_MIN"
    --reverse-pattern "${TOPUP_REVERSE_PE_PATTERN:-}"
    --fieldmap-pattern "${FIELDMAP_PATTERN:-_fieldmap}"
    --magnitude-pattern "${FIELDMAP_MAGNITUDE_PATTERN:-_magnitude}"
)
[[ -n "$LIST" ]] && args+=(--list "$LIST")
[[ "$KEEP_EXISTING" -eq 1 ]] && args+=(--keep-existing)

python3 "${SCRIPT_DIR}/scan_multirun.py" "${args[@]}"
