#!/bin/bash
# =============================================================================
#  STEP 0: Discover subjects and sessions from BIDS directory
# =============================================================================
#
#  This script:
#    1. Sources pipeline_config.cfg (auto-detects location)
#    2. Validates BIDS_ROOT exists
#    3. Finds all sub-XX/ses-YY pairs with func/ directories
#    4. Checks for .nii or .nii.gz functional files (raw names only — output
#       an older version left in the BIDS directory is not counted)
#    5. Discovers available tasks dynamically (no hardcoding)
#    6. Writes subses_list.txt with format: sub-XX ses-YY
#    7. Prints clear summary of what was found
#

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

# Source config
source "$CONFIG_FILE"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================

if [[ -z "${BIDS_ROOT:-}" ]]; then
    echo "ERROR: BIDS_ROOT not set in config file"
    exit 1
fi

if [[ ! -d "$BIDS_ROOT" ]]; then
    echo "ERROR: BIDS_ROOT does not exist: $BIDS_ROOT"
    exit 1
fi

# ============================================================================
# DISCOVER SUBJECTS AND SESSIONS
# ============================================================================

echo "=========================================="
echo "Step 0: Discovering Subjects and Sessions"
echo "=========================================="
echo "BIDS_ROOT: $BIDS_ROOT"
echo ""

SUBSES_LIST="${SCRIPT_DIR}/subses_list.txt"

# Create temporary file for output
TEMP_LIST=$(mktemp)
trap "rm -f $TEMP_LIST" EXIT

# Counter for statistics
declare -i total_subses=0
declare -i total_files=0
# Task labels, one per line. A plain string rather than an associative array so
# the script also runs on bash 3 (macOS), where "declare -A" does not exist.
found_tasks=

# Find all subject directories (sub-XX)
while IFS= read -r -d '' sub_dir; do
    sub_id=$(basename "$sub_dir")

    # Find all session directories (ses-YY) within this subject
    if [[ -d "$sub_dir" ]]; then
        # Check if this subject has sessions
        if ls -d "$sub_dir"/ses-* &> /dev/null; then
            # Has sessions
            while IFS= read -r -d '' ses_dir; do
                ses_id=$(basename "$ses_dir")

                # Check if func directory exists
                if [[ -d "$ses_dir/func" ]]; then
                    # Check for functional files (.nii or .nii.gz)
                    nii_count=$(find "$ses_dir/func" -maxdepth 1 \( -name "${sub_id}_*_bold.nii" -o -name "${sub_id}_*_bold.nii.gz" \) 2>/dev/null | wc -l)

                    if [[ $nii_count -gt 0 ]]; then
                        # This is a valid sub-ses pair
                        echo "$sub_id $ses_id" >> "$TEMP_LIST"
                        ((total_subses++)) || true
                        ((total_files+=nii_count)) || true

                        # Discover task labels from BOLD filenames
                        while IFS= read -r -d '' bold_file; do
                            filename=$(basename "$bold_file" .nii.gz)
                            filename=$(basename "$filename" .nii)

                            # Extract task label using pattern matching
                            if [[ $filename =~ task-([^_]+) ]]; then
                                found_tasks="${found_tasks}${BASH_REMATCH[1]}"$'\n'
                            fi
                        done < <(find "$ses_dir/func" -maxdepth 1 \( -name "${sub_id}_*_bold.nii" -o -name "${sub_id}_*_bold.nii.gz" \) -print0 2>/dev/null)
                    fi
                fi
            done < <(find "$sub_dir" -maxdepth 1 -type d -name "ses-*" -print0)
        else
            # No sessions, check subject-level func directory
            if [[ -d "$sub_dir/func" ]]; then
                nii_count=$(find "$sub_dir/func" -maxdepth 1 \( -name "${sub_id}_*_bold.nii" -o -name "${sub_id}_*_bold.nii.gz" \) 2>/dev/null | wc -l)

                if [[ $nii_count -gt 0 ]]; then
                    # Single-session subject
                    echo "$sub_id" >> "$TEMP_LIST"
                    ((total_subses++)) || true
                    ((total_files+=nii_count)) || true

                    # Discover task labels
                    while IFS= read -r -d '' bold_file; do
                        filename=$(basename "$bold_file" .nii.gz)
                        filename=$(basename "$filename" .nii)

                        if [[ $filename =~ task-([^_]+) ]]; then
                            found_tasks="${found_tasks}${BASH_REMATCH[1]}"$'\n'
                        fi
                    done < <(find "$sub_dir/func" -maxdepth 1 \( -name "${sub_id}_*_bold.nii" -o -name "${sub_id}_*_bold.nii.gz" \) -print0 2>/dev/null)
                fi
            fi
        fi
    fi
done < <(find "$BIDS_ROOT" -maxdepth 1 -type d -name "sub-*" -print0)

# ============================================================================
# WRITE OUTPUT AND SUMMARY
# ============================================================================

# Sort and write final list
sort "$TEMP_LIST" > "$SUBSES_LIST"

echo "Results:"
echo "--------"
echo "  Subject-session pairs found: $total_subses"
echo "  Total BOLD files: $total_files"
echo ""

if [[ -n "$found_tasks" ]]; then
    echo "  Tasks discovered:"
    printf '%s' "$found_tasks" | sort -u | sed '/^$/d; s/^/    - /'
    echo ""
fi

echo "Config file settings:"
echo "  TASKS (from config): ${TASKS:-not set}"
echo "  EXCLUDE_TASKS (from config): ${EXCLUDE_TASKS:-none}"
echo ""

if [[ $total_subses -eq 0 ]]; then
    echo "WARNING: No valid subject-session pairs found!"
    echo "  Check BIDS_ROOT path and BIDS structure."
    exit 1
fi

echo "Output written to: $SUBSES_LIST"
echo ""
echo "First 5 entries:"
head -5 "$SUBSES_LIST"

if [[ $total_subses -gt 5 ]]; then
    echo "  ... and $((total_subses - 5)) more"
fi

echo ""
echo "Step 0 complete. Ready for preprocessing."
echo ""

exit 0
