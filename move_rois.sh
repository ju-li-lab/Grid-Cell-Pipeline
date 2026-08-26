#!/usr/bin/env bash
# ======================================================================
# MOVE_ROIS.SH
# ======================================================================
# Copies ROI masks from a source directory into the correct BIDS
# anat/ folders, renaming them to match the pipeline's ROI_PATTERN.
#
# Source files are expected to follow the naming convention:
#   sub-XX_ses-YY_[run-N_]<LABEL>_<hemi>.nii[.gz]
# e.g.:
#   sub-01s13_ses-01_ErC_left.nii.gz
#   sub-04s13_ses-01_run-2_ErC_left.nii.gz
#
# Target names use the ROI_PATTERN_LEFT / ROI_PATTERN_RIGHT from config:
#   sub-XX_ses-YY_<ROI_PATTERN_LEFT>
# e.g.:
#   sub-01s13_ses-01_hemi-left_label-ErC_mask.nii
#
# The script can also decompress .nii.gz to .nii if the pipeline
# expects uncompressed NIfTI files.
#
# Usage:
#   bash move_rois.sh              # dry-run (preview)
#   bash move_rois.sh --execute    # actually copy and rename
#
# All settings come from pipeline_config.cfg.
# ======================================================================

set -euo pipefail

# ======================================================================
# Colors
# ======================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# ======================================================================
# Setup
# ======================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo -e "${RED}ERROR: Config file not found: $CONFIG_FILE${NC}"
    exit 1
fi

source "$CONFIG_FILE"

# Check required config
for var in BIDS_ROOT ROI_SOURCE_DIR ROI_SOURCE_LABEL ROI_PATTERN_LEFT ROI_PATTERN_RIGHT; do
    if [[ -z "${!var:-}" ]]; then
        echo -e "${RED}ERROR: $var not set in $CONFIG_FILE${NC}"
        exit 1
    fi
done

if [[ ! -d "$ROI_SOURCE_DIR" ]]; then
    echo -e "${RED}ERROR: ROI source directory not found: $ROI_SOURCE_DIR${NC}"
    exit 1
fi

if [[ ! -d "$BIDS_ROOT" ]]; then
    echo -e "${RED}ERROR: BIDS_ROOT not found: $BIDS_ROOT${NC}"
    exit 1
fi

# Parse arguments
DRY_RUN=1
DECOMPRESS=0
if [[ "${1:-}" == "--execute" ]] || [[ "${1:-}" == "-x" ]]; then
    DRY_RUN=0
fi
# Check if ROI patterns end in .nii (not .nii.gz) — need decompression
if [[ "$ROI_PATTERN_LEFT" == *.nii ]] && [[ "$ROI_PATTERN_LEFT" != *.nii.gz ]]; then
    DECOMPRESS=1
fi

echo "========================================"
echo "  ROI Mask Mover"
echo "========================================"
echo "Source:       $ROI_SOURCE_DIR"
echo "BIDS target:  $BIDS_ROOT"
echo "ROI label:    $ROI_SOURCE_LABEL"
echo "Pattern L:    $ROI_PATTERN_LEFT"
echo "Pattern R:    $ROI_PATTERN_RIGHT"
echo "Decompress:   $( [[ $DECOMPRESS -eq 1 ]] && echo 'yes (.nii.gz → .nii)' || echo 'no' )"
if [[ $DRY_RUN -eq 1 ]]; then
    echo -e "${YELLOW}Mode:         DRY RUN (preview only)${NC}"
    echo "  Use --execute or -x to actually copy files"
else
    echo -e "${GREEN}Mode:         EXECUTE (will copy files)${NC}"
fi
echo "========================================"
echo ""

# ======================================================================
# Counters
# ======================================================================
declare -i copied=0
declare -i skipped=0
declare -i errors=0

# ======================================================================
# Process each ROI file
# ======================================================================
process_roi() {
    local src_file="$1"
    local base
    base="$(basename "$src_file")"

    # Parse filename: sub-XX_ses-YY_[run-N_]<LABEL>_<hemi>.nii[.gz]
    local sub ses run_part label hemi ext

    # Get extension
    if [[ "$base" == *.nii.gz ]]; then
        ext=".nii.gz"
        local base_noext="${base%.nii.gz}"
    elif [[ "$base" == *.nii ]]; then
        ext=".nii"
        local base_noext="${base%.nii}"
    else
        echo -e "  ${YELLOW}[SKIP]${NC} Not a NIfTI file: $base"
        ((skipped++)) || true
        return
    fi

    # Extract sub-XX
    if [[ "$base_noext" =~ ^(sub-[^_]+) ]]; then
        sub="${BASH_REMATCH[1]}"
    else
        echo -e "  ${YELLOW}[SKIP]${NC} Cannot parse subject: $base"
        ((skipped++)) || true
        return
    fi

    # Extract ses-YY
    if [[ "$base_noext" =~ _(ses-[^_]+) ]]; then
        ses="${BASH_REMATCH[1]}"
    else
        echo -e "  ${YELLOW}[SKIP]${NC} Cannot parse session: $base"
        ((skipped++)) || true
        return
    fi

    # Extract hemisphere (last component before extension)
    if [[ "$base_noext" =~ _(left|right|lh|rh|L|R)$ ]]; then
        hemi="${BASH_REMATCH[1]}"
    else
        echo -e "  ${YELLOW}[SKIP]${NC} Cannot parse hemisphere: $base"
        ((skipped++)) || true
        return
    fi

    # Normalise hemisphere
    case "$hemi" in
        left|lh|L) hemi="left" ;;
        right|rh|R) hemi="right" ;;
    esac

    # Check target anat directory
    local target_dir="$BIDS_ROOT/$sub/$ses/anat"
    if [[ ! -d "$target_dir" ]]; then
        echo -e "  ${YELLOW}[SKIP]${NC} No anat/ dir: $target_dir (from $base)"
        ((skipped++)) || true
        return
    fi

    # Build target filename
    local target_pattern
    if [[ "$hemi" == "left" ]]; then
        target_pattern="$ROI_PATTERN_LEFT"
    else
        target_pattern="$ROI_PATTERN_RIGHT"
    fi

    local target_name="${sub}_${ses}${target_pattern}"
    local target_path="$target_dir/$target_name"

    # Check if target already exists
    if [[ -f "$target_path" ]]; then
        echo -e "  ${BLUE}[EXISTS]${NC} $target_path"
        ((skipped++)) || true
        return
    fi

    # Execute
    if [[ $DRY_RUN -eq 1 ]]; then
        if [[ $DECOMPRESS -eq 1 ]] && [[ "$ext" == ".nii.gz" ]]; then
            echo -e "  ${BLUE}[DRY]${NC} gunzip + copy: $base → $target_name"
        else
            echo -e "  ${BLUE}[DRY]${NC} copy: $base → $target_name"
        fi
    else
        if [[ $DECOMPRESS -eq 1 ]] && [[ "$ext" == ".nii.gz" ]]; then
            # Decompress to target
            gunzip -c "$src_file" > "$target_path"
            echo -e "  ${GREEN}[OK]${NC} Decompressed: $base → $target_name"
        else
            cp -n "$src_file" "$target_path"
            echo -e "  ${GREEN}[OK]${NC} Copied: $base → $target_name"
        fi
    fi
    ((copied++)) || true
}

# Find all NIfTI files in source directory that match the ROI label
# Patterns: sub-*_ses-*_<LABEL>_left.nii*, sub-*_ses-*_run-*_<LABEL>_left.nii*, etc.
found_files=0
while IFS= read -r -d '' f; do
    ((found_files++)) || true
    process_roi "$f"
done < <(find "$ROI_SOURCE_DIR" -type f \( \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_left.nii" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_left.nii.gz" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_right.nii" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_right.nii.gz" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_lh.nii" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_lh.nii.gz" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_rh.nii" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_rh.nii.gz" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_L.nii" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_L.nii.gz" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_R.nii" -o \
    -name "sub-*_ses-*_*${ROI_SOURCE_LABEL}_R.nii.gz" \
\) -print0 2>/dev/null | sort -z)

# ======================================================================
# Summary
# ======================================================================
echo ""
echo "========================================"
echo "  Summary"
echo "========================================"
echo "  Files found:   $found_files"
echo "  Copied:        $copied"
echo "  Skipped:       $skipped"
echo "  Errors:        $errors"
echo ""

if [[ $found_files -eq 0 ]]; then
    echo -e "${RED}No ROI files found matching *${ROI_SOURCE_LABEL}_{left,right}.nii* in:${NC}"
    echo "  $ROI_SOURCE_DIR"
    echo ""
    echo "Check ROI_SOURCE_DIR and ROI_SOURCE_LABEL in pipeline_config.cfg"
    exit 1
fi

if [[ $DRY_RUN -eq 1 ]]; then
    echo -e "${YELLOW}This was a dry run. To actually copy files, run:${NC}"
    echo "  bash $0 --execute"
fi
echo ""
