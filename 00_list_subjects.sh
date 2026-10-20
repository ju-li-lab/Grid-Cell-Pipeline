#!/usr/bin/env bash
# =============================================================================
#  Discover subject-session pairs, and apply the filters from config.sh
# =============================================================================
#  Walks BIDS_ROOT for sub-XX/ses-YY pairs that actually have BOLD data, then
#  narrows that list by INCLUDE_/EXCLUDE_SUBJECTS, INCLUDE_/EXCLUDE_SESSIONS
#  and REQUIRE_ROI.
#
#  Writes:
#    subses_list.txt       the list every later stage reads, one "sub ses" per line
#    subses_list_full.txt  the unfiltered list, for reference
#
#  Usage:
#    bash 00_list_subjects.sh              discover and filter
#    bash 00_list_subjects.sh --no-filter  discover only
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SELF_DIR}/config.sh"

APPLY_FILTERS=true
[[ "${1:-}" == "--no-filter" ]] && APPLY_FILTERS=false

[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config.sh not found: $CONFIG_FILE" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"
# shellcheck disable=SC1091
source "${SELF_DIR}/bids_common.sh"

: "${BIDS_ROOT:?ERROR: BIDS_ROOT not set in config.sh}"
[[ -d "$BIDS_ROOT" ]] || { echo "ERROR: BIDS_ROOT does not exist: $BIDS_ROOT" >&2; exit 1; }

SUBSES_LIST="${SELF_DIR}/subses_list.txt"
FULL_LIST="${SELF_DIR}/subses_list_full.txt"

WORK=$(mktemp)
trap 'rm -f "$WORK" "$WORK".tmp' EXIT

# -----------------------------------------------------------------------------
#  Discover
# -----------------------------------------------------------------------------
#  A pair counts only if its func/ holds raw BOLD named for that subject.
#  Matching on the raw name means output an older version of this pipeline
#  left inside BIDS_ROOT is not mistaken for input.

echo "Scanning ${BIDS_ROOT}"

# Task labels collected as plain text, not an associative array, so this also
# runs under bash 3 on macOS where "declare -A" does not exist.
found_tasks=
n_files=0

# bold_files <dir> <sub>  -> NUL-separated raw BOLD paths
bold_files() {
    find "$1" -maxdepth 1 \( -name "$2_*_bold.nii" -o -name "$2_*_bold.nii.gz" \) -print0 2>/dev/null
}

# record <func-dir> <sub> <label-to-write>
record() {
    local func_dir="$1" sub="$2" label="$3" count
    [[ -d "$func_dir" ]] || return 0
    count=$(bold_files "$func_dir" "$sub" | tr -dc '\0' | wc -c)
    [[ "$count" -gt 0 ]] || return 0

    echo "$label" >> "$WORK"
    n_files=$((n_files + count))

    local f stem
    while IFS= read -r -d '' f; do
        stem=$(basename "$f"); stem=${stem%.gz}; stem=${stem%.nii}
        [[ $stem =~ task-([^_]+) ]] && found_tasks="${found_tasks}${BASH_REMATCH[1]}"$'\n'
    done < <(bold_files "$func_dir" "$sub")
}

while IFS= read -r -d '' sub_dir; do
    sub=$(basename "$sub_dir")
    if compgen -G "$sub_dir/ses-*" > /dev/null; then
        while IFS= read -r -d '' ses_dir; do
            record "$ses_dir/func" "$sub" "$sub $(basename "$ses_dir")"
        done < <(find "$sub_dir" -maxdepth 1 -type d -name "ses-*" -print0)
    else
        record "$sub_dir/func" "$sub" "$sub"      # dataset without sessions
    fi
done < <(find "$BIDS_ROOT" -maxdepth 1 -type d -name "sub-*" -print0)

sort -o "$WORK" "$WORK"
cp "$WORK" "$FULL_LIST"
n_found=$(grep -cve '^[[:space:]]*$' "$WORK" || true)

echo "  Found ${n_found} subject-session pairs, ${n_files} BOLD files"
if [[ -n "$found_tasks" ]]; then
    echo "  Tasks present: $(printf '%s' "$found_tasks" | sort -u | sed '/^$/d' | paste -sd' ' -)"
    echo "  TASKS in config.sh: ${TASKS:-<not set>}"
fi

if [[ "$n_found" -eq 0 ]]; then
    echo "ERROR: no subject-session pair in ${BIDS_ROOT} has BOLD data in func/." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
#  Filter
# -----------------------------------------------------------------------------

if [[ "$APPLY_FILTERS" != "true" ]]; then
    cp "$WORK" "$SUBSES_LIST"
    echo "Wrote ${SUBSES_LIST} (${n_found} pairs, --no-filter)"
    exit 0
fi

# keep_or_drop <keep|drop> <comma list> <subject|session> <label>
#   subject matches at the start of the line, session at the end.
keep_or_drop() {
    local mode="$1" list="$2" field="$3" label="$4"
    [[ -n "$list" ]] || return 0

    local pattern= item
    local IFS=','
    for item in $list; do
        item="${item//[[:space:]]/}"
        [[ -n "$item" ]] || continue
        [[ -n "$pattern" ]] && pattern="${pattern}|"
        if [[ "$field" == "subject" ]]; then
            pattern="${pattern}^${item}( |\$)"
        else
            pattern="${pattern} ${item}\$"
        fi
    done
    unset IFS
    [[ -n "$pattern" ]] || return 0

    if [[ "$mode" == "keep" ]]; then
        grep -E "$pattern" "$WORK" > "$WORK.tmp" || true
    else
        grep -vE "$pattern" "$WORK" > "$WORK.tmp" || true
    fi
    mv "$WORK.tmp" "$WORK"
    printf '  %-18s %-28s -> %s remaining\n' "$label" "$list" "$(grep -cve '^[[:space:]]*$' "$WORK" || true)"
}

echo "Filtering"
keep_or_drop keep "${INCLUDE_SUBJECTS:-}" subject INCLUDE_SUBJECTS
keep_or_drop keep "${INCLUDE_SESSIONS:-}" session INCLUDE_SESSIONS
keep_or_drop drop "${EXCLUDE_SUBJECTS:-}" subject EXCLUDE_SUBJECTS
keep_or_drop drop "${EXCLUDE_SESSIONS:-}" session EXCLUDE_SESSIONS

if [[ "${REQUIRE_ROI:-false}" == "true" ]]; then
    if [[ -z "${ROI_PATTERN_LEFT:-}${ROI_PATTERN_RIGHT:-}" ]]; then
        echo "  REQUIRE_ROI=true but no ROI_PATTERN_LEFT/RIGHT set - skipping this filter."
    else
        roi_deriv="$(deriv_rois 2>/dev/null || true)"
        : > "$WORK.tmp"
        dropped=0
        while read -r sub ses _; do
            found=false
            # move_rois.sh imports masks into the derivatives; older datasets
            # may still have them beside the raw anatomy, so check both.
            for anat in ${roi_deriv:+"$(deriv_session "$roi_deriv" "$sub" "$ses" anat)"} \
                        "${BIDS_ROOT}/${sub}/${ses}/anat"; do
                [[ -d "$anat" ]] || continue
                for pat in "${ROI_PATTERN_LEFT:-}" "${ROI_PATTERN_RIGHT:-}"; do
                    [[ -n "$pat" ]] || continue
                    if compgen -G "$anat/*$pat*" > /dev/null; then found=true; break; fi
                done
                [[ "$found" == "true" ]] && break
            done
            if [[ "$found" == "true" ]]; then
                echo "$sub $ses" >> "$WORK.tmp"
            else
                dropped=$((dropped + 1))
            fi
        done < "$WORK"
        mv "$WORK.tmp" "$WORK"
        printf '  %-18s %-28s -> %s remaining\n' REQUIRE_ROI "dropped $dropped without masks" \
               "$(grep -cve '^[[:space:]]*$' "$WORK" || true)"
    fi
fi

cp "$WORK" "$SUBSES_LIST"
n_kept=$(grep -cve '^[[:space:]]*$' "$SUBSES_LIST" || true)

echo
echo "Wrote ${SUBSES_LIST}: ${n_kept} of ${n_found} pairs kept"
echo "      ${FULL_LIST}: the unfiltered list"
if [[ "$n_kept" -eq 0 ]]; then
    echo "ERROR: every pair was filtered out. Check the INCLUDE_/EXCLUDE_ settings." >&2
    exit 1
fi
head -5 "$SUBSES_LIST" | sed 's/^/      /'
[[ "$n_kept" -gt 5 ]] && echo "      ... and $((n_kept - 5)) more"
exit 0
