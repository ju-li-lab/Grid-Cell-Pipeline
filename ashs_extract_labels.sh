#!/usr/bin/env bash
# ==============================================================================
# ASHS Label Extraction Script
# ==============================================================================
# Extracts specific region labels from ASHS segmentation output into
# individual binary mask files.
#
# Usage:
#   bash ashs_extract_labels.sh                     # uses default config
#   bash ashs_extract_labels.sh my_config.sh        # uses custom config
#
# Requirements:
#   - Python 3 with nibabel and numpy (standard on neuroimaging HPCs)
#   - ASHS output directory with final/ segmentation files
#
# Author: Generated for GridCell Analysis pipeline
# Date:   2026-02-06
# ==============================================================================
set -uo pipefail

# ---- Color output helpers ----------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

info()    { echo -e "${CYAN}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERROR]${NC} $*" >&2; }
header()  { echo -e "\n${BOLD}━━━ $* ━━━${NC}"; }

# ---- Default label name lookup -----------------------------------------------
declare -A DEFAULT_LABEL_NAMES=(
    [0]="Background"
    [1]="CA1"
    [2]="CA2"
    [3]="DG"
    [4]="CA3"
    [5]="Tail"
    [6]="Hippocampal_sulcus"
    [7]="MISC"
    [8]="SUB"
    [9]="ERC"
    [10]="BA35"
    [11]="BA36"
    [12]="PHC"
    [13]="MISC2"
    [14]="CS"
)

# ---- Load configuration -----------------------------------------------------
CONFIG_FILE="${1:-$(dirname "$0")/ashs_extract_config.sh}"

if [[ ! -f "$CONFIG_FILE" ]]; then
    error "Config file not found: $CONFIG_FILE"
    echo "  Copy ashs_extract_config.sh next to this script and edit it."
    exit 1
fi

info "Loading config: $CONFIG_FILE"
source "$CONFIG_FILE"

# ---- Validate configuration --------------------------------------------------
header "Validating configuration"

# Check Python + nibabel
if ! python3 -c "import nibabel, numpy" 2>/dev/null; then
    error "Python 3 with nibabel and numpy is required."
    echo "  On HPC try:  module load python  or  pip install --user nibabel numpy"
    exit 1
fi
success "Python 3 + nibabel + numpy available"

# Check ASHS_DIR
if [[ -z "${ASHS_DIR:-}" ]] || [[ ! -d "$ASHS_DIR" ]]; then
    error "ASHS_DIR does not exist: ${ASHS_DIR:-<not set>}"
    exit 1
fi
success "ASHS directory: $ASHS_DIR"

# Check SEG_TYPE
case "${SEG_TYPE:-}" in
    corr_usegray|corr_nogray|heur) ;;
    *)
        error "Invalid SEG_TYPE: '${SEG_TYPE:-}'. Must be: corr_usegray, corr_nogray, or heur"
        exit 1
        ;;
esac
success "Segmentation type: $SEG_TYPE"

# Check LABELS
if [[ -z "${LABELS_TO_EXTRACT:-}" ]]; then
    error "LABELS_TO_EXTRACT is empty. Set at least one label number."
    exit 1
fi
# Convert to array
read -ra LABEL_ARRAY <<< "$LABELS_TO_EXTRACT"
success "Labels to extract: ${LABEL_ARRAY[*]}"

# Parse label names
read -ra NAME_ARRAY <<< "${LABEL_NAMES:-}"
if [[ ${#NAME_ARRAY[@]} -gt 0 ]] && [[ ${#NAME_ARRAY[@]} -ne ${#LABEL_ARRAY[@]} ]]; then
    error "LABEL_NAMES count (${#NAME_ARRAY[@]}) does not match LABELS_TO_EXTRACT count (${#LABEL_ARRAY[@]})"
    exit 1
fi

# Build label-to-name mapping
declare -A EXTRACT_NAMES
for i in "${!LABEL_ARRAY[@]}"; do
    lbl="${LABEL_ARRAY[$i]}"
    if [[ ${#NAME_ARRAY[@]} -gt 0 ]]; then
        EXTRACT_NAMES[$lbl]="${NAME_ARRAY[$i]}"
    else
        EXTRACT_NAMES[$lbl]="${DEFAULT_LABEL_NAMES[$lbl]:-label${lbl}}"
    fi
done

info "Label name mapping:"
for lbl in "${LABEL_ARRAY[@]}"; do
    echo "     Label $lbl → ${EXTRACT_NAMES[$lbl]}"
done

# Check hemispheres
case "${HEMISPHERES:-both}" in
    left)  HEMI_LIST=("left") ;;
    right) HEMI_LIST=("right") ;;
    both)  HEMI_LIST=("left" "right") ;;
    *)
        error "Invalid HEMISPHERES: '${HEMISPHERES:-}'. Must be: left, right, or both"
        exit 1
        ;;
esac
success "Hemispheres: ${HEMI_LIST[*]}"

# Check output mode
case "${OUTPUT_MODE:-flat}" in
    flat)
        if [[ -z "${FLAT_OUTPUT_DIR:-}" ]]; then
            error "OUTPUT_MODE is 'flat' but FLAT_OUTPUT_DIR is not set."
            exit 1
        fi
        ;;
    bids)
        if [[ -z "${BIDS_OUTPUT_DIR:-}" ]]; then
            error "OUTPUT_MODE is 'bids' but BIDS_OUTPUT_DIR is not set."
            exit 1
        fi
        if [[ ! -d "$BIDS_OUTPUT_DIR" ]]; then
            error "BIDS_OUTPUT_DIR does not exist: $BIDS_OUTPUT_DIR"
            exit 1
        fi
        ;;
    *)
        error "Invalid OUTPUT_MODE: '${OUTPUT_MODE:-}'. Must be: flat or bids"
        exit 1
        ;;
esac
success "Output mode: $OUTPUT_MODE"

# ---- Discover subjects and sessions -----------------------------------------
header "Discovering subjects and sessions"

# Build subject list
if [[ -n "${SUBJECTS:-}" ]]; then
    read -ra SUB_LIST <<< "$SUBJECTS"
else
    SUB_LIST=()
    for d in "$ASHS_DIR"/sub-*/; do
        [[ -d "$d" ]] && SUB_LIST+=("$(basename "$d")")
    done
fi

if [[ ${#SUB_LIST[@]} -eq 0 ]]; then
    error "No subjects found in $ASHS_DIR"
    exit 1
fi
info "Found ${#SUB_LIST[@]} subject(s): ${SUB_LIST[*]}"

# ---- Counters ----------------------------------------------------------------
TOTAL_EXTRACTED=0
TOTAL_SKIPPED=0
TOTAL_ERRORS=0
SKIPPED_FILES=()
ERROR_FILES=()

# ---- Python extraction helper (embedded) -------------------------------------
# This function writes a temporary Python script and runs it.
# It extracts a single label from a NIfTI file and saves as binary mask.
extract_label_python() {
    local input_file="$1"
    local output_file="$2"
    local label_value="$3"

    python3 << PYEOF
import nibabel as nib
import numpy as np
import sys

try:
    img = nib.load("${input_file}")
    data = img.get_fdata()

    # Create binary mask for the requested label
    mask = (data == ${label_value}).astype(np.uint8)

    voxel_count = int(np.sum(mask))
    if voxel_count == 0:
        print(f"WARNING: Label ${label_value} has 0 voxels in ${input_file}", file=sys.stderr)
        sys.exit(2)  # special exit code for empty label

    # Save with same header/affine as original
    out_img = nib.Nifti1Image(mask, img.affine, img.header)
    out_img.header.set_data_dtype(np.uint8)
    nib.save(out_img, "${output_file}")

    print(f"{voxel_count}")
    sys.exit(0)
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
}

# Combined multi-label extraction
extract_combined_python() {
    local input_file="$1"
    local output_file="$2"
    shift 2
    local labels=("$@")
    local labels_str=$(IFS=,; echo "${labels[*]}")

    python3 << PYEOF
import nibabel as nib
import numpy as np
import sys

try:
    img = nib.load("${input_file}")
    data = img.get_fdata()
    labels = [${labels_str}]

    # Keep original label values for requested labels, zero everything else
    mask = np.zeros_like(data, dtype=np.uint8)
    for lbl in labels:
        mask[data == lbl] = lbl

    voxel_count = int(np.sum(mask > 0))

    out_img = nib.Nifti1Image(mask, img.affine, img.header)
    out_img.header.set_data_dtype(np.uint8)
    nib.save(out_img, "${output_file}")

    print(f"{voxel_count}")
    sys.exit(0)
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
}

# ---- Helper: choose run directory --------------------------------------------
choose_run_dir() {
    local sub="$1"
    local ses="$2"
    local anat_dir="$ASHS_DIR/$sub/$ses/anat"

    if [[ ! -d "$anat_dir" ]]; then
        echo ""
        return
    fi

    # Find all run directories that contain a final/ folder
    local run_dirs=()
    for rd in "$anat_dir"/*/; do
        [[ -d "$rd/final" ]] && run_dirs+=("$(basename "$rd")")
    done

    if [[ ${#run_dirs[@]} -eq 0 ]]; then
        echo ""
        return
    fi

    if [[ ${#run_dirs[@]} -eq 1 ]]; then
        echo "${run_dirs[0]}"
        return
    fi

    # Multiple runs found — handle according to MULTI_RUN_MODE
    case "${MULTI_RUN_MODE:-ask}" in
        first)
            warn "$sub/$ses has ${#run_dirs[@]} runs. Auto-selecting first: ${run_dirs[0]}"
            echo "${run_dirs[0]}"
            ;;
        all)
            # Return all, pipe-separated
            local IFS='|'
            echo "${run_dirs[*]}"
            ;;
        ask)
            echo ""
            warn "$sub/$ses has ${#run_dirs[@]} run directories:"
            for i in "${!run_dirs[@]}"; do
                echo "     [$((i+1))] ${run_dirs[$i]}"
            done
            echo -n "  Enter choice (number), 'all', or 'skip': "
            read -r choice
            case "$choice" in
                skip|s|S)
                    echo "SKIP"
                    ;;
                all|a|A)
                    local IFS='|'
                    echo "${run_dirs[*]}"
                    ;;
                *)
                    if [[ "$choice" =~ ^[0-9]+$ ]] && (( choice >= 1 && choice <= ${#run_dirs[@]} )); then
                        echo "${run_dirs[$((choice-1))]}"
                    else
                        warn "Invalid choice. Skipping $sub/$ses."
                        echo "SKIP"
                    fi
                    ;;
            esac
            ;;
    esac
}

# ---- Helper: determine output path -------------------------------------------
get_output_path() {
    local sub="$1"
    local ses="$2"
    local hemi="$3"
    local label_name="$4"
    local run_tag="$5"  # empty string or e.g. "_run-t1norun_t2norun"

    local ext="${OUTPUT_FORMAT:-nii.gz}"
    local fname="${sub}_${ses}_hemi-${hemi}_label-${label_name}${run_tag}_mask.${ext}"

    case "${OUTPUT_MODE:-flat}" in
        flat)
            echo "${FLAT_OUTPUT_DIR}/${fname}"
            ;;
        bids)
            local bids_anat="${BIDS_OUTPUT_DIR}/${sub}/${ses}/anat"
            if [[ ! -d "$bids_anat" ]]; then
                warn "BIDS anat dir does not exist, creating: $bids_anat"
                if [[ "${DRY_RUN:-no}" != "yes" ]]; then
                    mkdir -p "$bids_anat"
                fi
            fi
            echo "${bids_anat}/${fname}"
            ;;
    esac
}

get_combined_output_path() {
    local sub="$1"
    local ses="$2"
    local hemi="$3"
    local run_tag="$4"

    local ext="${OUTPUT_FORMAT:-nii.gz}"
    # Build combined label name
    local combined_name=""
    for lbl in "${LABEL_ARRAY[@]}"; do
        combined_name+="${EXTRACT_NAMES[$lbl]}+"
    done
    combined_name="${combined_name%+}"  # remove trailing +

    local fname="${sub}_${ses}_hemi-${hemi}_label-${combined_name}${run_tag}_mask.${ext}"

    case "${OUTPUT_MODE:-flat}" in
        flat)
            echo "${FLAT_OUTPUT_DIR}/${fname}"
            ;;
        bids)
            local bids_anat="${BIDS_OUTPUT_DIR}/${sub}/${ses}/anat"
            mkdir -p "$bids_anat" 2>/dev/null || true
            echo "${bids_anat}/${fname}"
            ;;
    esac
}

# ---- Main processing loop ---------------------------------------------------
header "Starting label extraction"

# Create output directory if flat mode
if [[ "${OUTPUT_MODE}" == "flat" ]] && [[ "${DRY_RUN:-no}" != "yes" ]]; then
    mkdir -p "$FLAT_OUTPUT_DIR"
fi

for sub in "${SUB_LIST[@]}"; do
    # Verify subject directory exists
    if [[ ! -d "$ASHS_DIR/$sub" ]]; then
        warn "Subject directory not found: $ASHS_DIR/$sub — skipping"
        continue
    fi

    # Build session list for this subject
    if [[ -n "${SESSIONS:-}" ]]; then
        read -ra SES_LIST <<< "$SESSIONS"
    else
        SES_LIST=()
        for d in "$ASHS_DIR/$sub"/ses-*/; do
            [[ -d "$d" ]] && SES_LIST+=("$(basename "$d")")
        done
    fi

    if [[ ${#SES_LIST[@]} -eq 0 ]]; then
        warn "$sub: no sessions found — skipping"
        continue
    fi

    for ses in "${SES_LIST[@]}"; do
        header "$sub / $ses"

        # Check session directory exists
        if [[ ! -d "$ASHS_DIR/$sub/$ses" ]]; then
            warn "$sub/$ses does not exist — skipping"
            continue
        fi

        # Choose run directory
        run_result=$(choose_run_dir "$sub" "$ses")

        if [[ -z "$run_result" ]]; then
            warn "$sub/$ses: no ASHS final/ output found — skipping"
            TOTAL_SKIPPED=$((TOTAL_SKIPPED + 1))
            SKIPPED_FILES+=("$sub/$ses (no final/ directory)")
            continue
        fi

        if [[ "$run_result" == "SKIP" ]]; then
            info "$sub/$ses: skipped by user"
            TOTAL_SKIPPED=$((TOTAL_SKIPPED + 1))
            SKIPPED_FILES+=("$sub/$ses (user skipped)")
            continue
        fi

        # Split pipe-separated runs if mode=all
        IFS='|' read -ra RUN_LIST <<< "$run_result"

        for run_dir in "${RUN_LIST[@]}"; do
            local_final="$ASHS_DIR/$sub/$ses/anat/$run_dir/final"

            # Build run tag for filename (only if multiple runs)
            if [[ ${#RUN_LIST[@]} -gt 1 ]]; then
                # Convert t1-run01_t2-run01 -> run-t1run01t2run01
                run_tag="_run-$(echo "$run_dir" | tr -d '-' | tr '_' '-')"
            else
                run_tag=""
            fi

            info "Processing: $sub/$ses  run=$run_dir"

            for hemi in "${HEMI_LIST[@]}"; do
                # Determine input filename
                case "$SEG_TYPE" in
                    corr_usegray) seg_suffix="lfseg_corr_usegray" ;;
                    corr_nogray)  seg_suffix="lfseg_corr_nogray" ;;
                    heur)         seg_suffix="lfseg_heur" ;;
                esac

                # Try multiple filename conventions in order:
                #   1) sub-XX_ses-YY_<hemi>_<seg>.nii.gz  (session in filename)
                #   2) sub-XX_ses-YY_<hemi>_<seg>.nii
                #   3) sub-XX_<hemi>_<seg>.nii.gz         (session only in path)
                #   4) sub-XX_<hemi>_<seg>.nii
                input_file=""
                for candidate in \
                    "${local_final}/${sub}_${ses}_${hemi}_${seg_suffix}.nii.gz" \
                    "${local_final}/${sub}_${ses}_${hemi}_${seg_suffix}.nii" \
                    "${local_final}/${sub}_${hemi}_${seg_suffix}.nii.gz" \
                    "${local_final}/${sub}_${hemi}_${seg_suffix}.nii"
                do
                    if [[ -f "$candidate" ]]; then
                        input_file="$candidate"
                        break
                    fi
                done

                if [[ -z "$input_file" ]]; then
                    warn "  Segmentation not found: ${sub}_(${ses}_)?${hemi}_${seg_suffix}.nii(.gz) in $local_final"
                    TOTAL_ERRORS=$((TOTAL_ERRORS + 1))
                    ERROR_FILES+=("$sub/$ses/$hemi — file not found")
                    continue
                fi

                # Extract each label
                for lbl in "${LABEL_ARRAY[@]}"; do
                    label_name="${EXTRACT_NAMES[$lbl]}"
                    output_file=$(get_output_path "$sub" "$ses" "$hemi" "$label_name" "$run_tag")

                    # Check overwrite
                    if [[ -f "$output_file" ]] && [[ "${OVERWRITE:-no}" != "yes" ]]; then
                        warn "  Output exists, skipping (set OVERWRITE=yes to replace): $(basename "$output_file")"
                        TOTAL_SKIPPED=$((TOTAL_SKIPPED + 1))
                        continue
                    fi

                    if [[ "${DRY_RUN:-no}" == "yes" ]]; then
                        info "  [DRY RUN] Would extract label $lbl ($label_name) from $hemi → $(basename "$output_file")"
                        continue
                    fi

                    # Run extraction
                    voxels=$(extract_label_python "$input_file" "$output_file" "$lbl" 2>&1) || {
                        exit_code=$?
                        if [[ $exit_code -eq 2 ]]; then
                            warn "  Label $lbl ($label_name) $hemi: 0 voxels — file saved but empty"
                            # Still save the file (it was created by Python)
                            TOTAL_EXTRACTED=$((TOTAL_EXTRACTED + 1))
                        else
                            error "  Failed to extract label $lbl ($label_name) from $hemi: $voxels"
                            TOTAL_ERRORS=$((TOTAL_ERRORS + 1))
                            ERROR_FILES+=("$sub/$ses/$hemi label-$lbl — extraction failed")
                        fi
                        continue
                    }

                    success "  $hemi label $lbl ($label_name): $voxels voxels → $(basename "$output_file")"
                    TOTAL_EXTRACTED=$((TOTAL_EXTRACTED + 1))
                done

                # Combined mask if requested
                if [[ "${COMBINE_LABELS:-no}" == "yes" ]]; then
                    combined_output=$(get_combined_output_path "$sub" "$ses" "$hemi" "$run_tag")

                    if [[ -f "$combined_output" ]] && [[ "${OVERWRITE:-no}" != "yes" ]]; then
                        warn "  Combined output exists, skipping: $(basename "$combined_output")"
                    elif [[ "${DRY_RUN:-no}" == "yes" ]]; then
                        info "  [DRY RUN] Would create combined mask → $(basename "$combined_output")"
                    else
                        voxels=$(extract_combined_python "$input_file" "$combined_output" "${LABEL_ARRAY[@]}" 2>&1) || {
                            error "  Failed to create combined mask for $hemi: $voxels"
                            TOTAL_ERRORS=$((TOTAL_ERRORS + 1))
                            continue
                        }
                        success "  $hemi combined mask: $voxels voxels → $(basename "$combined_output")"
                        TOTAL_EXTRACTED=$((TOTAL_EXTRACTED + 1))
                    fi
                fi
            done
        done
    done
done

# ---- Summary -----------------------------------------------------------------
header "Extraction complete"
echo ""
info "Extracted:  $TOTAL_EXTRACTED files"
info "Skipped:    $TOTAL_SKIPPED"
info "Errors:     $TOTAL_ERRORS"

if [[ ${#SKIPPED_FILES[@]} -gt 0 ]]; then
    echo ""
    warn "Skipped items:"
    for s in "${SKIPPED_FILES[@]}"; do
        echo "     - $s"
    done
fi

if [[ ${#ERROR_FILES[@]} -gt 0 ]]; then
    echo ""
    error "Error items:"
    for e in "${ERROR_FILES[@]}"; do
        echo "     - $e"
    done
fi

echo ""
if [[ "${OUTPUT_MODE}" == "flat" ]]; then
    info "Output directory: $FLAT_OUTPUT_DIR"
else
    info "Output BIDS directory: $BIDS_OUTPUT_DIR"
fi

# Write extraction log
LOG_DIR="${FLAT_OUTPUT_DIR:-${BIDS_OUTPUT_DIR}}"
LOG_FILE="${LOG_DIR}/extraction_log_$(date +%Y%m%d_%H%M%S).txt"
if [[ "${DRY_RUN:-no}" != "yes" ]] && [[ -d "$LOG_DIR" ]]; then
    {
        echo "ASHS Label Extraction Log"
        echo "========================="
        echo "Date:           $(date)"
        echo "Config:         $CONFIG_FILE"
        echo "ASHS dir:       $ASHS_DIR"
        echo "Seg type:       $SEG_TYPE"
        echo "Labels:         ${LABELS_TO_EXTRACT}"
        echo "Hemispheres:    ${HEMISPHERES}"
        echo "Output mode:    ${OUTPUT_MODE}"
        echo "Output dir:     ${LOG_DIR}"
        echo ""
        echo "Results:"
        echo "  Extracted:    $TOTAL_EXTRACTED"
        echo "  Skipped:      $TOTAL_SKIPPED"
        echo "  Errors:       $TOTAL_ERRORS"
        if [[ ${#ERROR_FILES[@]} -gt 0 ]]; then
            echo ""
            echo "Errors:"
            for e in "${ERROR_FILES[@]}"; do
                echo "  - $e"
            done
        fi
    } > "$LOG_FILE"
    info "Log saved: $LOG_FILE"
fi

echo ""
success "Done."
