#!/bin/bash
# ======================================================================
# SCAN_MULTIRUN.SH
# ======================================================================
# Scans the BIDS directory for subjects/sessions that have:
#   - Multiple T2w runs (run-1_T2w, run-2_T2w, ...)
#   - Multiple runs of the same task BOLD (task-run2_run-1_bold, ...)
#   - Multiple reverse-PE runs (task-reverse_run-1_bold, ...)
#
# Produces run_selection.tsv — a tab-separated file listing every
# subject-session that has duplicate runs.  The user edits this file
# to choose which run to use, then the pipeline reads it.
#
# Usage:
#   bash scan_multirun.sh [path/to/pipeline_config.cfg]
#
# Output:
#   <SCRIPT_DIR>/run_selection.tsv
# ======================================================================

set -euo pipefail

# ======================================================================
# Setup
# ======================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -ge 1 ]]; then
    CONFIG_FILE="$1"
else
    CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

source "$CONFIG_FILE"

if [[ -z "${BIDS_ROOT:-}" ]] || [[ ! -d "$BIDS_ROOT" ]]; then
    echo "ERROR: BIDS_ROOT not set or does not exist: ${BIDS_ROOT:-<not set>}"
    exit 1
fi

OUTPUT_FILE="${SCRIPT_DIR}/run_selection.tsv"

# ======================================================================
# Colors
# ======================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# ======================================================================
# Scan BIDS directory
# ======================================================================
echo "========================================"
echo "  Multi-Run Scanner"
echo "========================================"
echo "BIDS_ROOT: $BIDS_ROOT"
echo ""

# Track statistics
declare -i multi_t2w_count=0
declare -i multi_bold_count=0
declare -i multi_reverse_count=0
declare -i multi_fieldmap_count=0
declare -i total_entries=0

# Start building TSV content
# Header
TSV_HEADER="subject\tsession\ttype\tavailable_runs\tselected_run"
TSV_LINES=()

for sub_dir in "$BIDS_ROOT"/sub-*/; do
    [[ -d "$sub_dir" ]] || continue
    SUB=$(basename "$sub_dir")

    for ses_dir in "$sub_dir"/ses-*/; do
        [[ -d "$ses_dir" ]] || continue
        SES=$(basename "$ses_dir")
        anat_dir="$ses_dir/anat"
        func_dir="$ses_dir/func"
        fmap_dir="$ses_dir/fmap"

        # ----- Check for multi-run T2w -----
        if [[ -d "$anat_dir" ]]; then
            # Find all T2w NIfTI files
            T2W_FILES=()
            while IFS= read -r -d '' f; do
                T2W_FILES+=("$(basename "$f")")
            done < <(find "$anat_dir" -maxdepth 1 \( -name "*_T2w.nii" -o -name "*_T2w.nii.gz" \) -print0 2>/dev/null)

            if [[ ${#T2W_FILES[@]} -gt 1 ]]; then
                # Multiple T2w runs found
                ((multi_t2w_count++)) || true
                ((total_entries++)) || true

                # Extract run labels
                runs=""
                for f in "${T2W_FILES[@]}"; do
                    # Extract run-N from filename
                    if [[ "$f" =~ run-([0-9]+) ]]; then
                        run_label="run-${BASH_REMATCH[1]}"
                    else
                        run_label="no-run-label"
                    fi
                    if [[ -n "$runs" ]]; then runs+=","; fi
                    runs+="$run_label"
                done
                # Sort runs
                runs=$(echo "$runs" | tr ',' '\n' | sort | tr '\n' ',' | sed 's/,$//')

                TSV_LINES+=("${SUB}\t${SES}\tT2w\t${runs}\t")
            fi
        fi

        # ----- Check for multi-run task BOLD -----
        if [[ -d "$func_dir" ]]; then
            # Get unique task labels from BOLD files
            declare -A task_files
            while IFS= read -r -d '' f; do
                fname=$(basename "$f")
                # Extract task label: everything between task- and the next _bold or _run-
                if [[ "$fname" =~ task-([^_]+)(_run-[0-9]+)?_bold ]]; then
                    task_label="${BASH_REMATCH[1]}"
                    run_part="${BASH_REMATCH[2]:-}"

                    # Skip if this is a reverse-PE EPI (handled separately)
                    if [[ "$task_label" == "reverse" ]]; then
                        continue
                    fi

                    # Track: key = task label, value = list of full filenames
                    if [[ -v task_files["$task_label"] ]]; then
                        task_files["$task_label"]+=",$fname"
                    else
                        task_files["$task_label"]="$fname"
                    fi
                fi
            done < <(find "$func_dir" -maxdepth 1 \( -name "*_bold.nii" -o -name "*_bold.nii.gz" \) -print0 2>/dev/null)

            # For each task, check if there are multiple runs
            for task_label in "${!task_files[@]}"; do
                IFS=',' read -ra files <<< "${task_files[$task_label]}"
                if [[ ${#files[@]} -gt 1 ]]; then
                    ((multi_bold_count++)) || true
                    ((total_entries++)) || true

                    # Extract run labels
                    runs=""
                    for f in "${files[@]}"; do
                        if [[ "$f" =~ _run-([0-9]+)_bold ]]; then
                            run_label="run-${BASH_REMATCH[1]}"
                        else
                            run_label="no-run-label"
                        fi
                        if [[ -n "$runs" ]]; then runs+=","; fi
                        runs+="$run_label"
                    done
                    runs=$(echo "$runs" | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')

                    TSV_LINES+=("${SUB}\t${SES}\ttask-${task_label}\t${runs}\t")
                fi
            done
            unset task_files

            # ----- Check for multi-run reverse-PE -----
            REVERSE_FILES=()
            while IFS= read -r -d '' f; do
                REVERSE_FILES+=("$(basename "$f")")
            done < <(find "$func_dir" -maxdepth 1 \( -name "*task-reverse*_bold.nii" -o -name "*task-reverse*_bold.nii.gz" \) -print0 2>/dev/null)

            if [[ ${#REVERSE_FILES[@]} -gt 1 ]]; then
                ((multi_reverse_count++)) || true
                ((total_entries++)) || true

                runs=""
                for f in "${REVERSE_FILES[@]}"; do
                    if [[ "$f" =~ _run-([0-9]+) ]]; then
                        run_label="run-${BASH_REMATCH[1]}"
                    else
                        run_label="no-run-label"
                    fi
                    if [[ -n "$runs" ]]; then runs+=","; fi
                    runs+="$run_label"
                done
                runs=$(echo "$runs" | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')

                TSV_LINES+=("${SUB}\t${SES}\ttask-reverse\t${runs}\t")
            fi
        fi

        # ----- Check for multi-run fieldmaps (PREPROC_MODE=precalc_fieldmap) -----
        # These are what make_fieldmaps.sh writes. The name must start with
        # <sub>_<ses> so leftovers such as vdm5_*_fieldmap.nii are ignored.
        if [[ -d "$fmap_dir" ]]; then
            FIELDMAP_FILES=()
            while IFS= read -r -d '' f; do
                FIELDMAP_FILES+=("$(basename "$f")")
            done < <(find "$fmap_dir" -maxdepth 1 \
                        \( -name "${SUB}_${SES}*_fieldmap.nii" -o -name "${SUB}_${SES}*_fieldmap.nii.gz" \) \
                        -print0 2>/dev/null)

            if [[ ${#FIELDMAP_FILES[@]} -gt 1 ]]; then
                ((multi_fieldmap_count++)) || true
                ((total_entries++)) || true

                runs=""
                for f in "${FIELDMAP_FILES[@]}"; do
                    if [[ "$f" =~ _run-([0-9]+) ]]; then
                        run_label="run-${BASH_REMATCH[1]}"
                    else
                        run_label="no-run-label"
                    fi
                    if [[ -n "$runs" ]]; then runs+=","; fi
                    runs+="$run_label"
                done
                runs=$(echo "$runs" | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')

                TSV_LINES+=("${SUB}\t${SES}\tfieldmap\t${runs}\t")
            fi
        fi
    done
done

# ======================================================================
# Write output
# ======================================================================

if [[ ${#TSV_LINES[@]} -eq 0 ]]; then
    echo -e "${GREEN}No multi-run cases found. All subjects/sessions have single runs.${NC}"
    echo ""
    echo "No run_selection.tsv needed."
    exit 0
fi

# Write the TSV
{
    echo -e "$TSV_HEADER"
    for line in "${TSV_LINES[@]}"; do
        echo -e "$line"
    done
} > "$OUTPUT_FILE"

# ======================================================================
# Summary
# ======================================================================
echo -e "${YELLOW}Multi-run cases found:${NC}"
echo "  T2w multi-run:       $multi_t2w_count subject-sessions"
echo "  Task BOLD multi-run: $multi_bold_count subject-session-tasks"
echo "  Reverse multi-run:   $multi_reverse_count subject-sessions"
echo "  Fieldmap multi-run:  $multi_fieldmap_count subject-sessions"
echo ""
echo "  Total entries:       $total_entries"
echo ""
echo -e "${BLUE}Output written to:${NC} $OUTPUT_FILE"
echo ""
echo "========================================"
echo "  WHAT TO DO NEXT"
echo "========================================"
echo ""
echo "1. Open run_selection.tsv in a text editor or spreadsheet program"
echo "2. For each row, fill in the 'selected_run' column with the run you want"
echo "   (e.g., 'run-1' or 'run-2')"
echo "3. Save the file"
echo "4. The pipeline will use your selections when processing"
echo ""
echo "If 'selected_run' is left empty for an entry, the pipeline will:"
echo "  - For T2w: use the LAST run (highest run number)"
echo "  - For task BOLD: use the LAST run (highest run number)"
echo "  - For reverse-PE: use the FIRST run (lowest run number)"
echo "  - For fieldmap:   use the LAST run (highest run number)"
echo ""
echo -e "${YELLOW}Preview of run_selection.tsv:${NC}"
echo "---"
column -t -s $'\t' "$OUTPUT_FILE" | head -20
if [[ $total_entries -gt 19 ]]; then
    echo "  ... ($((total_entries - 19)) more rows)"
fi
echo "---"
