#!/bin/bash
# =============================================================================
#  Filter subses_list.txt based on config settings
# =============================================================================
#
#  Reads the full subses_list.txt (from step0) and applies subject/session
#  filters defined in pipeline_config.cfg. The original list is backed up
#  as subses_list_full.txt; the filtered result overwrites subses_list.txt.
#
#  Filters (all optional, leave empty to skip):
#    INCLUDE_SUBJECTS  — keep only these subjects (whitelist)
#    INCLUDE_SESSIONS  — keep only these sessions (whitelist)
#    EXCLUDE_SUBJECTS  — remove these subjects (blacklist)
#    EXCLUDE_SESSIONS  — remove these sessions (blacklist)
#    REQUIRE_ROI       — remove pairs with no ROI mask, looking in the ROI
#                        derivatives dataset first and the raw anat/ second
#
#  Usage:
#    bash filter_subses_list.sh
#
#  Or use run_pipeline.sh, which calls this automatically after step 0.
# =============================================================================

set -euo pipefail

# ============================================================================
# AUTO-DETECT SCRIPT DIRECTORY AND SOURCE CONFIG
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

source "$CONFIG_FILE"
source "${SCRIPT_DIR}/bids_common.sh"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================

SUBSES_LIST="${SCRIPT_DIR}/subses_list.txt"

if [[ ! -f "$SUBSES_LIST" ]]; then
    echo "ERROR: subses_list.txt not found: $SUBSES_LIST"
    echo "       Run step0_make_subses_list.sh first."
    exit 1
fi

echo "=========================================="
echo "Filtering Subject-Session List"
echo "=========================================="
echo ""

# ============================================================================
# BACKUP ORIGINAL LIST
# ============================================================================

FULL_LIST="${SCRIPT_DIR}/subses_list_full.txt"
cp "$SUBSES_LIST" "$FULL_LIST"

TOTAL_BEFORE=$(wc -l < "$SUBSES_LIST")
echo "Input: $TOTAL_BEFORE subject-session pairs"
echo ""

# ============================================================================
# APPLY FILTERS
# ============================================================================

TEMP_LIST=$(mktemp)
trap "rm -f $TEMP_LIST" EXIT

cp "$SUBSES_LIST" "$TEMP_LIST"

# ---- INCLUDE_SUBJECTS (whitelist) ----
if [[ -n "${INCLUDE_SUBJECTS:-}" ]]; then
    echo "INCLUDE_SUBJECTS: ${INCLUDE_SUBJECTS}"
    # Build grep pattern: sub-01|sub-02|sub-03
    IFS=',' read -ra SUBS <<< "$INCLUDE_SUBJECTS"
    PATTERN=""
    for s in "${SUBS[@]}"; do
        s=$(echo "$s" | xargs)  # trim whitespace
        if [[ -n "$PATTERN" ]]; then
            PATTERN="${PATTERN}|"
        fi
        PATTERN="${PATTERN}^${s} "
    done
    grep -E "$PATTERN" "$TEMP_LIST" > "${TEMP_LIST}.tmp" || true
    mv "${TEMP_LIST}.tmp" "$TEMP_LIST"
    echo "  -> $(wc -l < "$TEMP_LIST") remaining"
fi

# ---- INCLUDE_SESSIONS (whitelist) ----
if [[ -n "${INCLUDE_SESSIONS:-}" ]]; then
    echo "INCLUDE_SESSIONS: ${INCLUDE_SESSIONS}"
    IFS=',' read -ra SESS <<< "$INCLUDE_SESSIONS"
    PATTERN=""
    for s in "${SESS[@]}"; do
        s=$(echo "$s" | xargs)
        if [[ -n "$PATTERN" ]]; then
            PATTERN="${PATTERN}|"
        fi
        PATTERN="${PATTERN} ${s}$"
    done
    grep -E "$PATTERN" "$TEMP_LIST" > "${TEMP_LIST}.tmp" || true
    mv "${TEMP_LIST}.tmp" "$TEMP_LIST"
    echo "  -> $(wc -l < "$TEMP_LIST") remaining"
fi

# ---- EXCLUDE_SUBJECTS (blacklist) ----
if [[ -n "${EXCLUDE_SUBJECTS:-}" ]]; then
    echo "EXCLUDE_SUBJECTS: ${EXCLUDE_SUBJECTS}"
    IFS=',' read -ra SUBS <<< "$EXCLUDE_SUBJECTS"
    PATTERN=""
    for s in "${SUBS[@]}"; do
        s=$(echo "$s" | xargs)
        if [[ -n "$PATTERN" ]]; then
            PATTERN="${PATTERN}|"
        fi
        PATTERN="${PATTERN}^${s} "
    done
    grep -vE "$PATTERN" "$TEMP_LIST" > "${TEMP_LIST}.tmp" || true
    mv "${TEMP_LIST}.tmp" "$TEMP_LIST"
    echo "  -> $(wc -l < "$TEMP_LIST") remaining"
fi

# ---- EXCLUDE_SESSIONS (blacklist) ----
if [[ -n "${EXCLUDE_SESSIONS:-}" ]]; then
    echo "EXCLUDE_SESSIONS: ${EXCLUDE_SESSIONS}"
    IFS=',' read -ra SESS <<< "$EXCLUDE_SESSIONS"
    PATTERN=""
    for s in "${SESS[@]}"; do
        s=$(echo "$s" | xargs)
        if [[ -n "$PATTERN" ]]; then
            PATTERN="${PATTERN}|"
        fi
        PATTERN="${PATTERN} ${s}$"
    done
    grep -vE "$PATTERN" "$TEMP_LIST" > "${TEMP_LIST}.tmp" || true
    mv "${TEMP_LIST}.tmp" "$TEMP_LIST"
    echo "  -> $(wc -l < "$TEMP_LIST") remaining"
fi

# ---- REQUIRE_ROI ----
if [[ "${REQUIRE_ROI:-false}" == "true" ]]; then
    echo "REQUIRE_ROI: checking for ROI masks..."

    ROI_PAT_L="${ROI_PATTERN_LEFT:-}"
    ROI_PAT_R="${ROI_PATTERN_RIGHT:-}"

    if [[ -z "$ROI_PAT_L" && -z "$ROI_PAT_R" ]]; then
        echo "  WARNING: REQUIRE_ROI=true but no ROI_PATTERN_LEFT/RIGHT set. Skipping."
    else
        ROI_DERIV="$(deriv_rois 2>/dev/null || echo '')"
        [[ -n "$ROI_DERIV" ]] && echo "  Looking in: ${ROI_DERIV}/<sub>/<ses>/anat"
        echo "          and ${BIDS_ROOT}/<sub>/<ses>/anat"

        > "${TEMP_LIST}.tmp"
        removed_roi=0
        while read -r sub ses rest; do
            found=false

            # move_rois.sh imports masks into the derivatives; older datasets
            # may still have them beside the raw anatomy, so check both.
            for anat_dir in ${ROI_DERIV:+"$(deriv_session "$ROI_DERIV" "$sub" "$ses" anat)"} \
                            "${BIDS_ROOT}/${sub}/${ses}/anat"; do
                [[ -d "$anat_dir" ]] || continue

                if [[ -n "$ROI_PAT_L" ]] && ls "$anat_dir"/*"$ROI_PAT_L"* &>/dev/null; then
                    found=true
                fi
                if [[ -n "$ROI_PAT_R" ]] && [[ "$found" == "false" ]] && \
                   ls "$anat_dir"/*"$ROI_PAT_R"* &>/dev/null; then
                    found=true
                fi
                [[ "$found" == "true" ]] && break
            done

            if [[ "$found" == "true" ]]; then
                echo "$sub $ses" >> "${TEMP_LIST}.tmp"
            else
                ((removed_roi++)) || true
            fi
        done < "$TEMP_LIST"
        mv "${TEMP_LIST}.tmp" "$TEMP_LIST"
        echo "  Removed $removed_roi pairs without ROI masks"
        echo "  -> $(wc -l < "$TEMP_LIST") remaining"
    fi
fi

# ============================================================================
# WRITE FILTERED LIST
# ============================================================================

TOTAL_AFTER=$(wc -l < "$TEMP_LIST")
REMOVED=$((TOTAL_BEFORE - TOTAL_AFTER))

cp "$TEMP_LIST" "$SUBSES_LIST"

echo ""
echo "=========================================="
echo "Result: $TOTAL_AFTER / $TOTAL_BEFORE kept ($REMOVED removed)"
echo "=========================================="
echo ""

if [[ $TOTAL_AFTER -eq 0 ]]; then
    echo "WARNING: All subject-session pairs were filtered out!"
    echo "  Check your INCLUDE/EXCLUDE settings in pipeline_config.cfg."
    echo "  Full list preserved in: $FULL_LIST"
    exit 1
fi

echo "Filtered list written to: $SUBSES_LIST"
echo "Full (unfiltered) backup: $FULL_LIST"
echo ""

if [[ $TOTAL_AFTER -le 10 ]]; then
    echo "Contents:"
    cat "$SUBSES_LIST"
else
    echo "First 10 entries:"
    head -10 "$SUBSES_LIST"
    echo "  ... and $((TOTAL_AFTER - 10)) more"
fi

echo ""
echo "Filtering complete."
