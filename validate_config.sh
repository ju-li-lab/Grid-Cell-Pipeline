#!/bin/bash
# ======================================================================
# VALIDATE_CONFIG.SH
# ======================================================================
# Validates the central pipeline configuration file.
# Checks paths, acquisition parameters, BIDS structure, and HPC settings.
#
# Usage:
#   bash validate_config.sh [path/to/pipeline_config.cfg]
#
# If no path given, looks for pipeline_config.cfg in the same directory.
#
# Output:
#   - Colored summary (green ✓, red ✗, yellow ⚠)
#   - Overall PASS/FAIL status at the end
#
# Exit codes:
#   0 = PASS (all critical checks succeeded)
#   1 = FAIL (one or more critical checks failed)
# ======================================================================

set -euo pipefail

# ======================================================================
# ANSI Color codes for output
# ======================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'  # No Color

CHECKMARK='✓'
CROSS='✗'
WARNING='⚠'

# ======================================================================
# Configuration
# ======================================================================

# Get script directory (where validate_config.sh lives)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Determine config file path
if [[ $# -ge 1 ]]; then
    CONFIG_FILE="$1"
else
    CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"
fi

# Tracking results
CHECKS_PASSED=0
CHECKS_FAILED=0
CHECKS_WARNED=0

# ======================================================================
# Helper Functions
# ======================================================================

print_header() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_check_pass() {
    echo -e "${GREEN}${CHECKMARK}${NC} $1"
    ((CHECKS_PASSED++)) || true
}

print_check_fail() {
    echo -e "${RED}${CROSS}${NC} $1"
    ((CHECKS_FAILED++)) || true
}

print_check_warn() {
    echo -e "${YELLOW}${WARNING}${NC} $1"
    ((CHECKS_WARNED++)) || true
}

print_detail() {
    echo "  → $1"
}

# ======================================================================
# Main Validation Script
# ======================================================================

echo ""
echo -e "${BLUE}GridCAT Pipeline Configuration Validator${NC}"
echo "Config file: ${CONFIG_FILE}"

# Check if config file exists
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo ""
    echo -e "${RED}ERROR: Configuration file not found!${NC}"
    echo "  Expected: $CONFIG_FILE"
    echo ""
    echo "Please ensure pipeline_config.cfg exists in:"
    echo "  - The directory containing this script, OR"
    echo "  - Pass the path as an argument: bash validate_config.sh /path/to/pipeline_config.cfg"
    exit 1
fi

# Source the config file
if ! source "$CONFIG_FILE" 2>/dev/null; then
    echo ""
    echo -e "${RED}ERROR: Cannot source configuration file!${NC}"
    echo "The file may contain syntax errors."
    exit 1
fi

print_header "1. CONFIGURATION FILE"

# Check if config file is readable
if [[ -r "$CONFIG_FILE" ]]; then
    print_check_pass "Configuration file is readable"
else
    print_check_fail "Configuration file is not readable (check permissions)"
fi

# Count configuration entries
CFG_ENTRIES=$(grep -c "^[^#].*=" "$CONFIG_FILE" || echo 0)
print_detail "Found $CFG_ENTRIES configuration entries"

# ======================================================================
# PATH VALIDATION
# ======================================================================

print_header "2. PATHS"

# Required path variables
declare -a REQUIRED_PATHS=(
    "BIDS_ROOT"
    "OUTPUT_ROOT"
    "SCRIPT_DIR"
    "SPM_DIR"
    "GRIDCAT_DIR"
)

for path_var in "${REQUIRED_PATHS[@]}"; do
    if [[ -z "${!path_var:-}" ]]; then
        print_check_fail "$path_var is not defined"
        continue
    fi

    path_value="${!path_var}"

    if [[ -d "$path_value" ]]; then
        print_check_pass "$path_var exists"
        print_detail "$path_value"
    else
        print_check_fail "$path_var does not exist or is not accessible"
        print_detail "$path_value"
    fi
done

# Check OUTPUT_ROOT is writable (or parent directory is)
if [[ -z "${OUTPUT_ROOT:-}" ]]; then
    print_check_warn "OUTPUT_ROOT not defined, skipping write check"
elif [[ -d "$OUTPUT_ROOT" ]]; then
    if [[ -w "$OUTPUT_ROOT" ]]; then
        print_check_pass "OUTPUT_ROOT is writable"
    else
        print_check_fail "OUTPUT_ROOT exists but is not writable"
    fi
else
    # Check if parent directory is writable
    OUTPUT_PARENT=$(dirname "$OUTPUT_ROOT")
    if [[ -d "$OUTPUT_PARENT" ]] && [[ -w "$OUTPUT_PARENT" ]]; then
        print_check_pass "OUTPUT_ROOT parent directory is writable (will be created)"
    else
        print_check_fail "OUTPUT_ROOT parent is not writable or doesn't exist"
        print_detail "$OUTPUT_PARENT"
    fi
fi

# ======================================================================
# ACQUISITION PARAMETERS
# ======================================================================

print_header "3. ACQUISITION PARAMETERS"

# Validate TOTAL_READOUT_MS
if [[ -z "${TOTAL_READOUT_MS:-}" ]]; then
    print_check_fail "TOTAL_READOUT_MS is not defined"
else
    if [[ "$TOTAL_READOUT_MS" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        if (( $(echo "$TOTAL_READOUT_MS > 0" | bc -l) )); then
            print_check_pass "TOTAL_READOUT_MS is valid"
            print_detail "$TOTAL_READOUT_MS ms"
        else
            print_check_fail "TOTAL_READOUT_MS must be positive (got $TOTAL_READOUT_MS)"
        fi
    else
        print_check_fail "TOTAL_READOUT_MS is not a valid number (got $TOTAL_READOUT_MS)"
    fi
fi

# Validate echo times
if [[ -z "${TE_SHORT_MS:-}" ]] || [[ -z "${TE_LONG_MS:-}" ]]; then
    print_check_fail "TE_SHORT_MS and/or TE_LONG_MS not defined"
else
    if [[ "$TE_SHORT_MS" =~ ^[0-9]+(\.[0-9]+)?$ ]] && [[ "$TE_LONG_MS" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        if (( $(echo "$TE_SHORT_MS < $TE_LONG_MS" | bc -l) )); then
            print_check_pass "Echo times are valid"
            print_detail "TE1=$TE_SHORT_MS ms, TE2=$TE_LONG_MS ms"
        else
            print_check_fail "TE_SHORT_MS should be less than TE_LONG_MS"
        fi
    else
        print_check_fail "Echo times are not valid numbers"
    fi
fi

# Validate BLIP_DIRECTION
if [[ -z "${BLIP_DIRECTION:-}" ]]; then
    print_check_fail "BLIP_DIRECTION is not defined"
else
    if [[ "$BLIP_DIRECTION" == "1" ]] || [[ "$BLIP_DIRECTION" == "-1" ]]; then
        print_check_pass "BLIP_DIRECTION is valid"
        print_detail "Value: $BLIP_DIRECTION"
    else
        print_check_fail "BLIP_DIRECTION must be 1 or -1 (got $BLIP_DIRECTION)"
    fi
fi

# Validate TR
if [[ -z "${TR:-}" ]]; then
    print_check_fail "TR is not defined"
else
    if [[ "$TR" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        if (( $(echo "$TR > 0" | bc -l) )); then
            print_check_pass "TR is valid"
            print_detail "$TR seconds"
        else
            print_check_fail "TR must be positive (got $TR)"
        fi
    else
        print_check_fail "TR is not a valid number (got $TR)"
    fi
fi

# ======================================================================
# STUDY DESIGN
# ======================================================================

print_header "4. STUDY DESIGN"

# Validate TASKS (accepts numeric: 1,2,3 or string labels: run1,run2)
if [[ -z "${TASKS:-}" ]]; then
    print_check_warn "TASKS is not defined (will process all tasks found)"
else
    if [[ "$TASKS" =~ ^[a-zA-Z0-9]+([,][a-zA-Z0-9]+)*$ ]]; then
        print_check_pass "TASKS is valid"
        print_detail "Tasks: $TASKS"
    else
        print_check_fail "TASKS format invalid (use comma-separated labels like: 1,2,3 or run1,run2)"
        print_detail "Got: $TASKS"
    fi
fi

# Validate PREPROC_MODE
if [[ -z "${PREPROC_MODE:-}" ]]; then
    print_check_fail "PREPROC_MODE is not defined"
else
    if [[ "$PREPROC_MODE" =~ ^(realign_unwarp|topup|realign_only)$ ]]; then
        print_check_pass "PREPROC_MODE is valid"
        print_detail "Mode: $PREPROC_MODE"
    else
        print_check_fail "PREPROC_MODE must be 'realign_unwarp', 'topup', or 'realign_only' (got $PREPROC_MODE)"
    fi
fi

# Validate ROI_MODE
if [[ -z "${ROI_MODE:-}" ]]; then
    print_check_warn "ROI_MODE is not defined (will default to 'both')"
else
    if [[ "$ROI_MODE" =~ ^(both|bilat_only|lr_only)$ ]]; then
        print_check_pass "ROI_MODE is valid"
        print_detail "Mode: $ROI_MODE"
    else
        print_check_fail "ROI_MODE must be 'both', 'bilat_only', or 'lr_only' (got $ROI_MODE)"
    fi
fi

# Validate ROI source settings (for move_rois.sh)
if [[ -n "${ROI_SOURCE_DIR:-}" ]]; then
    if [[ -d "$ROI_SOURCE_DIR" ]]; then
        print_check_pass "ROI_SOURCE_DIR exists"
        ROI_COUNT=$(find "$ROI_SOURCE_DIR" -type f -name "*.nii*" 2>/dev/null | wc -l)
        print_detail "$ROI_SOURCE_DIR ($ROI_COUNT NIfTI files)"
    else
        print_check_warn "ROI_SOURCE_DIR does not exist: $ROI_SOURCE_DIR"
    fi
else
    print_check_warn "ROI_SOURCE_DIR not set (move_rois.sh won't work without it)"
fi

if [[ -n "${ROI_SOURCE_LABEL:-}" ]]; then
    print_check_pass "ROI_SOURCE_LABEL is set"
    print_detail "Label: $ROI_SOURCE_LABEL"
else
    print_check_warn "ROI_SOURCE_LABEL not set"
fi

# Validate run selection file
if [[ -n "${RUN_SELECTION_FILE:-}" ]]; then
    if [[ -f "$RUN_SELECTION_FILE" ]]; then
        SEL_ROWS=$(tail -n +2 "$RUN_SELECTION_FILE" | grep -cve '^\s*$' || echo 0)
        print_check_pass "RUN_SELECTION_FILE exists"
        print_detail "$RUN_SELECTION_FILE ($SEL_ROWS entries)"
        # Check for unfilled entries
        UNFILLED=$(tail -n +2 "$RUN_SELECTION_FILE" | awk -F'\t' '{if ($5 == "" || $5 ~ /^\s*$/) print}' | wc -l)
        if [[ "$UNFILLED" -gt 0 ]]; then
            print_check_warn "$UNFILLED entries in run_selection.tsv have no selected_run (will use defaults)"
        fi
    else
        print_check_warn "RUN_SELECTION_FILE specified but not found: $RUN_SELECTION_FILE"
        print_detail "Run: bash scan_multirun.sh to generate it"
    fi
else
    print_check_pass "RUN_SELECTION_FILE not set (auto-select mode for multi-run cases)"
fi

# ======================================================================
# BIDS STRUCTURE
# ======================================================================

print_header "5. BIDS STRUCTURE"

if [[ -z "${BIDS_ROOT:-}" ]]; then
    print_check_fail "BIDS_ROOT not defined, cannot check BIDS structure"
elif [[ ! -d "$BIDS_ROOT" ]]; then
    print_check_fail "BIDS_ROOT directory does not exist"
    print_detail "$BIDS_ROOT"
else
    print_check_pass "BIDS_ROOT directory exists"

    # Check for sub-XX directories
    SUBJ_DIRS=$(find "$BIDS_ROOT" -maxdepth 1 -type d -name "sub-*" 2>/dev/null | wc -l)
    if [[ $SUBJ_DIRS -gt 0 ]]; then
        print_check_pass "Found BIDS subject directories"
        print_detail "$SUBJ_DIRS subject(s)"

        # Sample first subject for detailed checks
        FIRST_SUBJ=$(find "$BIDS_ROOT" -maxdepth 1 -type d -name "sub-*" 2>/dev/null | sort | head -1)
        if [[ -n "$FIRST_SUBJ" ]]; then
            SUBJ_NAME=$(basename "$FIRST_SUBJ")

            # Check for session directories
            SESS_DIRS=$(find "$FIRST_SUBJ" -maxdepth 1 -type d -name "ses-*" 2>/dev/null | wc -l)
            if [[ $SESS_DIRS -gt 0 ]]; then
                print_check_pass "Found session directories in first subject"
                print_detail "$SUBJ_NAME has $SESS_DIRS session(s)"

                FIRST_SESS=$(find "$FIRST_SUBJ" -maxdepth 1 -type d -name "ses-*" 2>/dev/null | sort | head -1)
                if [[ -n "$FIRST_SESS" ]]; then
                    SESS_NAME=$(basename "$FIRST_SESS")

                    # Check for func directory
                    if [[ -d "$FIRST_SESS/func" ]]; then
                        print_check_pass "Found 'func' directory"
                        FUNC_FILES=$(find "$FIRST_SESS/func" -name "*_bold.nii*" 2>/dev/null | wc -l)
                        print_detail "$FUNC_FILES functional file(s)"
                    else
                        print_check_warn "No 'func' directory found in $SUBJ_NAME/$SESS_NAME"
                    fi

                    # Check for anat directory
                    if [[ -d "$FIRST_SESS/anat" ]]; then
                        print_check_pass "Found 'anat' directory"
                        ANAT_FILES=$(find "$FIRST_SESS/anat" -name "*.nii*" 2>/dev/null | wc -l)
                        print_detail "$ANAT_FILES anatomical file(s)"

                        # Check T2w files
                        T2W_COUNT=$(find "$FIRST_SESS/anat" -name "*_T2w.nii*" 2>/dev/null | wc -l)
                        if [[ $T2W_COUNT -eq 1 ]]; then
                            print_check_pass "Found 1 T2w image"
                            print_detail "$(basename "$(find "$FIRST_SESS/anat" -name "*_T2w.nii*" 2>/dev/null | head -1)")"
                        elif [[ $T2W_COUNT -gt 1 ]]; then
                            print_check_warn "Found $T2W_COUNT T2w runs (multi-run detected)"
                            find "$FIRST_SESS/anat" -name "*_T2w.nii*" 2>/dev/null | while read -r tf; do
                                print_detail "  $(basename "$tf")"
                            done
                            if [[ -n "${RUN_SELECTION_FILE:-}" ]] && [[ -f "${RUN_SELECTION_FILE:-}" ]]; then
                                print_detail "Run selection file will determine which T2w to use"
                            else
                                print_detail "No run selection file — pipeline will use last T2w run"
                                print_detail "Run: bash scan_multirun.sh to generate run_selection.tsv"
                            fi
                        elif [[ $T2W_COUNT -eq 0 ]]; then
                            print_check_warn "No T2w images found (needed for ROI coregistration)"
                        fi

                        # Check for multi-run task BOLD
                        if [[ -d "$FIRST_SESS/func" ]]; then
                            MULTI_BOLD=$(find "$FIRST_SESS/func" -name "*_bold.nii*" 2>/dev/null | \
                                sed -n 's/.*task-\([^_]*\)_run-\([0-9]*\).*/\1/p' | sort -u)
                            if [[ -n "$MULTI_BOLD" ]]; then
                                print_check_warn "Multi-run task BOLD detected for: $MULTI_BOLD"
                                if [[ -z "${RUN_SELECTION_FILE:-}" ]] || [[ ! -f "${RUN_SELECTION_FILE:-}" ]]; then
                                    print_detail "No run selection file — pipeline will use last BOLD run"
                                fi
                            fi
                        fi
                    else
                        print_check_warn "No 'anat' directory found in $SUBJ_NAME/$SESS_NAME"
                    fi

                    # Check for fmap directory
                    if [[ -d "$FIRST_SESS/fmap" ]]; then
                        print_check_pass "Found 'fmap' directory (for distortion correction)"
                        FMAP_FILES=$(find "$FIRST_SESS/fmap" -name "*.nii*" 2>/dev/null | wc -l)
                        print_detail "$FMAP_FILES fieldmap file(s)"

                        # Check for reverse-PE EPI if topup mode
                        if [[ "$PREPROC_MODE" == "topup" ]] && [[ -n "${TOPUP_REVERSE_PE_PATTERN:-}" ]]; then
                            REVERSE_PE=""
                            SEARCH_DIR="${TOPUP_REVERSE_PE_DIR:-auto}"
                            # Determine which directories to search
                            if [[ "$SEARCH_DIR" == "fmap" ]] || [[ "$SEARCH_DIR" == "auto" ]]; then
                                REVERSE_PE=$(find "$FIRST_SESS/fmap" -name "*${TOPUP_REVERSE_PE_PATTERN}*" 2>/dev/null | head -1)
                            fi
                            if [[ -z "$REVERSE_PE" ]] && { [[ "$SEARCH_DIR" == "func" ]] || [[ "$SEARCH_DIR" == "auto" ]]; }; then
                                REVERSE_PE=$(find "$FIRST_SESS/func" -name "*${TOPUP_REVERSE_PE_PATTERN}*" 2>/dev/null | head -1)
                            fi
                            if [[ -f "$REVERSE_PE" ]]; then
                                print_check_pass "Found reverse-PE EPI for topup"
                                print_detail "$(basename "$REVERSE_PE") (in $(basename "$(dirname "$REVERSE_PE")")/)"
                            else
                                SEARCHED="fmap/"
                                [[ "$SEARCH_DIR" == "func" ]] && SEARCHED="func/"
                                [[ "$SEARCH_DIR" == "auto" ]] && SEARCHED="fmap/ and func/"
                                print_check_fail "Reverse-PE EPI not found for topup (pattern: *${TOPUP_REVERSE_PE_PATTERN}* in $SEARCHED)"
                            fi
                        fi
                    else
                        if [[ "$PREPROC_MODE" == "realign_unwarp" ]]; then
                            print_check_warn "No 'fmap' directory found, but PREPROC_MODE=$PREPROC_MODE"
                        elif [[ "$PREPROC_MODE" == "topup" ]]; then
                            SEARCH_DIR="${TOPUP_REVERSE_PE_DIR:-auto}"
                            if [[ "$SEARCH_DIR" != "func" ]]; then
                                print_check_warn "No 'fmap' directory found, but PREPROC_MODE=topup and TOPUP_REVERSE_PE_DIR=$SEARCH_DIR"
                            fi
                        fi
                    fi
                fi
            else
                print_check_warn "No session directories in $SUBJ_NAME (may be ses-level studies)"
            fi
        fi
    else
        print_check_fail "No BIDS subject directories (sub-*) found"
        print_detail "BIDS_ROOT: $BIDS_ROOT"
    fi

    # Check for required BIDS files
    if [[ -n "$FIRST_SUBJ" ]] && [[ -n "$FIRST_SESS" ]]; then
        # Check for JSON sidecars
        FUNC_JSON=$(find "$FIRST_SESS/func" -name "*_bold.json" 2>/dev/null | head -1)
        if [[ -f "$FUNC_JSON" ]]; then
            print_check_pass "Found functional JSON sidecar"
            print_detail "$(basename "$FUNC_JSON")"
        else
            print_check_warn "No functional JSON sidecar found (needed for TotalReadoutTime)"
        fi

        if [[ -d "$FIRST_SESS/fmap" ]]; then
            FMAP_JSON=$(find "$FIRST_SESS/fmap" -name "*_phasediff.json" 2>/dev/null | head -1)
            if [[ -f "$FMAP_JSON" ]]; then
                print_check_pass "Found fieldmap JSON sidecar"
                print_detail "$(basename "$FMAP_JSON")"
            else
                print_check_warn "No fieldmap JSON sidecar found"
            fi
        fi
    fi
fi

# ======================================================================
# SOFTWARE DIRECTORIES
# ======================================================================

print_header "6. SOFTWARE REQUIREMENTS"

# Check SPM
if [[ -z "${SPM_DIR:-}" ]]; then
    print_check_fail "SPM_DIR is not defined"
else
    if [[ -d "$SPM_DIR" ]]; then
        if [[ -f "$SPM_DIR/spm.m" ]]; then
            print_check_pass "SPM installation found"
            print_detail "$SPM_DIR"
        else
            print_check_fail "SPM_DIR exists but spm.m not found"
            print_detail "Expected: $SPM_DIR/spm.m"
        fi
    else
        print_check_fail "SPM_DIR does not exist or is not accessible"
        print_detail "$SPM_DIR"
    fi
fi

# Check GridCAT
if [[ -z "${GRIDCAT_DIR:-}" ]]; then
    print_check_fail "GRIDCAT_DIR is not defined"
else
    if [[ -d "$GRIDCAT_DIR" ]]; then
        if [[ -f "$GRIDCAT_DIR/specifyGLM.m" ]]; then
            print_check_pass "GridCAT installation found"
            print_detail "$GRIDCAT_DIR"
        else
            print_check_fail "GRIDCAT_DIR exists but specifyGLM.m not found"
            print_detail "Expected: $GRIDCAT_DIR/specifyGLM.m"
        fi
    else
        print_check_warn "GRIDCAT_DIR does not exist (may be installed elsewhere)"
        print_detail "$GRIDCAT_DIR"
    fi
fi

# ======================================================================
# TOPUP SETTINGS (only when PREPROC_MODE=topup)
# ======================================================================

if [[ "${PREPROC_MODE:-}" == "topup" ]]; then
    print_header "6b. FSL TOPUP SETTINGS"

    # Validate FSL module
    if [[ -z "${FSL_MODULE:-}" ]]; then
        print_check_warn "FSL_MODULE not defined (will try default: fsl)"
    else
        print_check_pass "FSL_MODULE is set"
        print_detail "Module: $FSL_MODULE"
    fi

    # Validate PE directions
    VALID_PE_DIRS="^(x|x-|y|y-|z|z-)$"

    if [[ -z "${TOPUP_PE_DIR_BOLD:-}" ]]; then
        print_check_fail "TOPUP_PE_DIR_BOLD is not defined"
    else
        if [[ "$TOPUP_PE_DIR_BOLD" =~ $VALID_PE_DIRS ]]; then
            print_check_pass "TOPUP_PE_DIR_BOLD is valid"
            print_detail "Direction: $TOPUP_PE_DIR_BOLD"
        else
            print_check_fail "TOPUP_PE_DIR_BOLD must be x, x-, y, y-, z, or z- (got $TOPUP_PE_DIR_BOLD)"
        fi
    fi

    if [[ -z "${TOPUP_PE_DIR_REVERSE:-}" ]]; then
        print_check_fail "TOPUP_PE_DIR_REVERSE is not defined"
    else
        if [[ "$TOPUP_PE_DIR_REVERSE" =~ $VALID_PE_DIRS ]]; then
            print_check_pass "TOPUP_PE_DIR_REVERSE is valid"
            print_detail "Direction: $TOPUP_PE_DIR_REVERSE"
        else
            print_check_fail "TOPUP_PE_DIR_REVERSE must be x, x-, y, y-, z, or z- (got $TOPUP_PE_DIR_REVERSE)"
        fi
    fi

    # Warn if forward and reverse PE directions don't look opposite
    if [[ -n "${TOPUP_PE_DIR_BOLD:-}" ]] && [[ -n "${TOPUP_PE_DIR_REVERSE:-}" ]]; then
        BOLD_BASE="${TOPUP_PE_DIR_BOLD//-/}"
        REV_BASE="${TOPUP_PE_DIR_REVERSE//-/}"
        if [[ "$BOLD_BASE" != "$REV_BASE" ]]; then
            print_check_warn "PE directions do not share the same axis ($TOPUP_PE_DIR_BOLD vs $TOPUP_PE_DIR_REVERSE)"
        elif [[ "$TOPUP_PE_DIR_BOLD" == "$TOPUP_PE_DIR_REVERSE" ]]; then
            print_check_fail "Forward and reverse PE directions are identical ($TOPUP_PE_DIR_BOLD)"
        else
            print_check_pass "Forward and reverse PE directions are opposite"
        fi
    fi

    # Validate readout time
    if [[ -z "${TOPUP_READOUT_SEC:-}" ]]; then
        print_check_fail "TOPUP_READOUT_SEC is not defined"
    else
        if [[ "$TOPUP_READOUT_SEC" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
            if (( $(echo "$TOPUP_READOUT_SEC > 0" | bc -l) )); then
                print_check_pass "TOPUP_READOUT_SEC is valid"
                print_detail "$TOPUP_READOUT_SEC seconds"
            else
                print_check_fail "TOPUP_READOUT_SEC must be positive"
            fi
        else
            print_check_fail "TOPUP_READOUT_SEC is not a valid number (got $TOPUP_READOUT_SEC)"
        fi
    fi

    # Validate reverse-PE search directory
    if [[ -z "${TOPUP_REVERSE_PE_DIR:-}" ]]; then
        print_check_warn "TOPUP_REVERSE_PE_DIR not set (defaults to auto)"
    else
        if [[ "$TOPUP_REVERSE_PE_DIR" =~ ^(fmap|func|auto)$ ]]; then
            print_check_pass "TOPUP_REVERSE_PE_DIR is valid"
            print_detail "Search directory: $TOPUP_REVERSE_PE_DIR"
        else
            print_check_fail "TOPUP_REVERSE_PE_DIR must be 'fmap', 'func', or 'auto' (got $TOPUP_REVERSE_PE_DIR)"
        fi
    fi

    # Validate reverse-PE pattern
    if [[ -z "${TOPUP_REVERSE_PE_PATTERN:-}" ]]; then
        print_check_fail "TOPUP_REVERSE_PE_PATTERN is not defined"
    else
        print_check_pass "TOPUP_REVERSE_PE_PATTERN is set"
        print_detail "Pattern: $TOPUP_REVERSE_PE_PATTERN (substring match)"
    fi

    # Validate topup config file (optional)
    if [[ -n "${TOPUP_CONFIG:-}" ]]; then
        if [[ -f "$TOPUP_CONFIG" ]]; then
            print_check_pass "Custom topup config file found"
            print_detail "$TOPUP_CONFIG"
        else
            print_check_fail "TOPUP_CONFIG file not found: $TOPUP_CONFIG"
        fi
    else
        print_check_pass "Using FSL default topup config (b02b0.cnf)"
    fi

    # Validate apply method
    if [[ -z "${TOPUP_APPLY_METHOD:-}" ]]; then
        print_check_warn "TOPUP_APPLY_METHOD not defined (will default to applytopup)"
    else
        if [[ "$TOPUP_APPLY_METHOD" =~ ^(applytopup|vdm)$ ]]; then
            print_check_pass "TOPUP_APPLY_METHOD is valid"
            print_detail "Method: $TOPUP_APPLY_METHOD"
        else
            print_check_fail "TOPUP_APPLY_METHOD must be 'applytopup' or 'vdm' (got $TOPUP_APPLY_METHOD)"
        fi
    fi

    # Validate interpolation
    if [[ -n "${TOPUP_INTERP:-}" ]]; then
        if [[ "$TOPUP_INTERP" =~ ^(spline|trilinear)$ ]]; then
            print_check_pass "TOPUP_INTERP is valid"
            print_detail "Interpolation: $TOPUP_INTERP"
        else
            print_check_fail "TOPUP_INTERP must be 'spline' or 'trilinear' (got $TOPUP_INTERP)"
        fi
    fi
fi

# ======================================================================
# GLM / GridCAT SETTINGS
# ======================================================================

print_header "7. GLM & GRIDCAT PARAMETERS"

# Validate X_FOLD_SYMMETRY
if [[ -z "${X_FOLD_SYMMETRY:-}" ]]; then
    print_check_warn "X_FOLD_SYMMETRY not defined (will use default: 6)"
else
    if [[ "$X_FOLD_SYMMETRY" =~ ^[0-9]+$ ]] && [[ $X_FOLD_SYMMETRY -gt 0 ]]; then
        print_check_pass "X_FOLD_SYMMETRY is valid"
        print_detail "Value: $X_FOLD_SYMMETRY"
    else
        print_check_fail "X_FOLD_SYMMETRY must be a positive integer (got $X_FOLD_SYMMETRY)"
    fi
fi

# Validate MASKING_THRESHOLD
if [[ -z "${MASKING_THRESHOLD:-}" ]]; then
    print_check_warn "MASKING_THRESHOLD not defined"
else
    if [[ "$MASKING_THRESHOLD" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        if (( $(echo "$MASKING_THRESHOLD >= 0 && $MASKING_THRESHOLD <= 1" | bc -l) )); then
            print_check_pass "MASKING_THRESHOLD is valid"
            print_detail "Value: $MASKING_THRESHOLD"
        else
            print_check_fail "MASKING_THRESHOLD must be between 0 and 1 (got $MASKING_THRESHOLD)"
        fi
    else
        print_check_fail "MASKING_THRESHOLD is not a valid number"
    fi
fi

# Validate HPF_CUTOFF
if [[ -z "${HPF_CUTOFF:-}" ]]; then
    print_check_warn "HPF_CUTOFF not defined (will use default)"
else
    if [[ "$HPF_CUTOFF" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        if (( $(echo "$HPF_CUTOFF > 0" | bc -l) )); then
            print_check_pass "HPF_CUTOFF is valid"
            print_detail "$HPF_CUTOFF seconds"
        else
            print_check_fail "HPF_CUTOFF must be positive"
        fi
    else
        print_check_fail "HPF_CUTOFF is not a valid number"
    fi
fi

# ======================================================================
# SLURM / HPC SETTINGS
# ======================================================================

print_header "8. HPC & SLURM CONFIGURATION"

# Validate SLURM_PARTITION
if [[ -z "${SLURM_PARTITION:-}" ]]; then
    print_check_warn "SLURM_PARTITION not defined"
else
    print_check_pass "SLURM_PARTITION is set"
    print_detail "Partition: $SLURM_PARTITION"
fi

# Validate time limits (HH:MM:SS format)
validate_time_format() {
    local time_str="$1"
    local time_name="$2"
    if [[ $time_str =~ ^[0-9]{1,2}:[0-9]{2}:[0-9]{2}$ ]]; then
        return 0
    else
        return 1
    fi
}

if [[ -z "${PREPROC_TIME:-}" ]]; then
    print_check_warn "PREPROC_TIME not defined"
else
    if validate_time_format "$PREPROC_TIME" "PREPROC_TIME"; then
        print_check_pass "PREPROC_TIME format is valid"
        print_detail "$PREPROC_TIME"
    else
        print_check_fail "PREPROC_TIME format invalid (use HH:MM:SS)"
        print_detail "Got: $PREPROC_TIME"
    fi
fi

if [[ -z "${GRIDCAT_TIME:-}" ]]; then
    print_check_warn "GRIDCAT_TIME not defined"
else
    if validate_time_format "$GRIDCAT_TIME" "GRIDCAT_TIME"; then
        print_check_pass "GRIDCAT_TIME format is valid"
        print_detail "$GRIDCAT_TIME"
    else
        print_check_fail "GRIDCAT_TIME format invalid (use HH:MM:SS)"
        print_detail "Got: $GRIDCAT_TIME"
    fi
fi

# Validate CPU counts
validate_cpu() {
    local cpu_str="$1"
    if [[ $cpu_str =~ ^[0-9]+$ ]] && [[ $cpu_str -gt 0 ]]; then
        return 0
    else
        return 1
    fi
}

if [[ -z "${PREPROC_CPUS:-}" ]]; then
    print_check_warn "PREPROC_CPUS not defined"
else
    if validate_cpu "$PREPROC_CPUS"; then
        print_check_pass "PREPROC_CPUS is valid"
        print_detail "$PREPROC_CPUS CPUs"
    else
        print_check_fail "PREPROC_CPUS must be a positive integer (got $PREPROC_CPUS)"
    fi
fi

if [[ -z "${GRIDCAT_CPUS:-}" ]]; then
    print_check_warn "GRIDCAT_CPUS not defined"
else
    if validate_cpu "$GRIDCAT_CPUS"; then
        print_check_pass "GRIDCAT_CPUS is valid"
        print_detail "$GRIDCAT_CPUS CPUs"
    else
        print_check_fail "GRIDCAT_CPUS must be a positive integer (got $GRIDCAT_CPUS)"
    fi
fi

# Validate MAX_WORKERS
if [[ -z "${MAX_WORKERS:-}" ]]; then
    print_check_warn "MAX_WORKERS not defined"
else
    if validate_cpu "$MAX_WORKERS"; then
        print_check_pass "MAX_WORKERS is valid"
        print_detail "$MAX_WORKERS workers"

        # Warn if MAX_WORKERS >= GRIDCAT_CPUS
        if [[ -n "${GRIDCAT_CPUS:-}" ]] && [[ $MAX_WORKERS -ge $GRIDCAT_CPUS ]]; then
            print_check_warn "MAX_WORKERS should be less than GRIDCAT_CPUS (currently $MAX_WORKERS >= $GRIDCAT_CPUS)"
        fi
    else
        print_check_fail "MAX_WORKERS must be a positive integer"
    fi
fi

# Validate memory format
if [[ -z "${PREPROC_MEM_PER_CPU:-}" ]]; then
    print_check_warn "PREPROC_MEM_PER_CPU not defined"
else
    if [[ "$PREPROC_MEM_PER_CPU" =~ ^[0-9]+[GMK]$ ]]; then
        print_check_pass "PREPROC_MEM_PER_CPU format is valid"
        print_detail "$PREPROC_MEM_PER_CPU"
    else
        print_check_fail "PREPROC_MEM_PER_CPU format invalid (use: 12G, 1024M, etc.)"
        print_detail "Got: $PREPROC_MEM_PER_CPU"
    fi
fi

if [[ -z "${GRIDCAT_MEM:-}" ]]; then
    print_check_warn "GRIDCAT_MEM not defined"
else
    if [[ "$GRIDCAT_MEM" =~ ^[0-9]+[GMK]$ ]]; then
        print_check_pass "GRIDCAT_MEM format is valid"
        print_detail "$GRIDCAT_MEM"
    else
        print_check_fail "GRIDCAT_MEM format invalid (use: 32G, 1024M, etc.)"
        print_detail "Got: $GRIDCAT_MEM"
    fi
fi

# Validate MATLAB_MODULE
if [[ -z "${MATLAB_MODULE:-}" ]]; then
    print_check_warn "MATLAB_MODULE not defined (will try default: matlab)"
else
    print_check_pass "MATLAB_MODULE is set"
    print_detail "Module: $MATLAB_MODULE"
fi

# ======================================================================
# DEBUG / FLAGS
# ======================================================================

print_header "9. DEBUG & FLAGS"

# Validate FAIL_ON_MISSING
if [[ -z "${FAIL_ON_MISSING:-}" ]]; then
    print_check_warn "FAIL_ON_MISSING not defined (will default to true)"
else
    if [[ "$FAIL_ON_MISSING" =~ ^(true|false)$ ]]; then
        print_check_pass "FAIL_ON_MISSING is valid"
        print_detail "Value: $FAIL_ON_MISSING"
    else
        print_check_fail "FAIL_ON_MISSING must be 'true' or 'false' (got $FAIL_ON_MISSING)"
    fi
fi

# Validate DRY_RUN
if [[ -z "${DRY_RUN:-}" ]]; then
    print_check_warn "DRY_RUN not defined (will default to false)"
else
    if [[ "$DRY_RUN" =~ ^(true|false)$ ]]; then
        print_check_pass "DRY_RUN is valid"
        print_detail "Value: $DRY_RUN"
    else
        print_check_fail "DRY_RUN must be 'true' or 'false' (got $DRY_RUN)"
    fi
fi

# Check COPY_MODE
if [[ -z "${COPY_MODE:-}" ]]; then
    print_check_warn "COPY_MODE not defined (will default to copy)"
else
    if [[ "$COPY_MODE" =~ ^(copy|symlink)$ ]]; then
        print_check_pass "COPY_MODE is valid"
        print_detail "Mode: $COPY_MODE"
    else
        print_check_fail "COPY_MODE must be 'copy' or 'symlink' (got $COPY_MODE)"
    fi
fi

# ======================================================================
# SUMMARY
# ======================================================================

print_header "VALIDATION SUMMARY"

TOTAL_CHECKS=$((CHECKS_PASSED + CHECKS_FAILED + CHECKS_WARNED))

echo ""
echo -e "Passed:  ${GREEN}$CHECKS_PASSED${NC}"
echo -e "Failed:  ${RED}$CHECKS_FAILED${NC}"
echo -e "Warned:  ${YELLOW}$CHECKS_WARNED${NC}"
echo -e "Total:   $TOTAL_CHECKS"
echo ""

# Final decision
if [[ $CHECKS_FAILED -eq 0 ]]; then
    echo -e "${GREEN}VALIDATION: PASS${NC}"
    echo ""
    echo "Your configuration is ready! You can now run the pipeline."
    echo ""
    if [[ $CHECKS_WARNED -gt 0 ]]; then
        echo -e "${YELLOW}Note: There are $CHECKS_WARNED warning(s). Review them above.${NC}"
    fi
    echo ""
    exit 0
else
    echo -e "${RED}VALIDATION: FAIL${NC}"
    echo ""
    echo "Please fix the $CHECKS_FAILED error(s) above before running the pipeline."
    echo ""
    exit 1
fi
