#!/usr/bin/env bash
set -euo pipefail

# ===================== USER SETTINGS =====================
ROI_SRC="/sc-projects/sc-proj-cc02-brainspace/sirius_wip/ROI"          # ROI files
BIDS_ROOT="/sc-projects/sc-proj-cc02-brainspace/sirius_wip/b2_bids"   # your BIDS root, containing sub-*/ses-*/anat

DRY_RUN=0                 # set to 0 to actually copy
RENAME_TO_BIDS_STYLE=0    # set to 1 to rename in target
# =========================================================

echo "ROI_SRC   = $ROI_SRC"
echo "BIDS_ROOT = $BIDS_ROOT"
echo "DRY_RUN   = $DRY_RUN"
echo "RENAME    = $RENAME_TO_BIDS_STYLE"
echo

copy_one () {
  local f="$1"
  local base
  base="$(basename "$f")"

  # Extract subject/session from filename
  local sub ses
  sub="$(echo "$base" | sed -nE 's/^(sub-[^_]+)_.*$/\1/p')"
  ses="$(echo "$base" | sed -nE 's/^sub-[^_]+_(ses-[^_]+)_.*$/\1/p')"

  if [[ -z "${sub:-}" || -z "${ses:-}" ]]; then
    echo "[SKIP] Could not parse sub/ses from: $base"
    return 0
  fi

  local target_dir="$BIDS_ROOT/$sub/$ses/anat"
  if [[ ! -d "$target_dir" ]]; then
    echo "[SKIP] Target anat folder missing: $target_dir  (from $base)"
    return 0
  fi

  local target_name="$base"

  if [[ "$RENAME_TO_BIDS_STYLE" -eq 1 ]]; then
    local ext base_noext
    if [[ "$base" == *.nii.gz ]]; then
      ext=".nii.gz"
      base_noext="${base%.nii.gz}"
    else
      ext=".nii"
      base_noext="${base%.nii}"
    fi

    local roi_label hemi
    roi_label="$(echo "$base_noext" | sed -nE 's/^sub-[^_]+_ses-[^_]+_T2w_(.*)_(left|right)$/\1/p')"
    hemi="$(echo "$base_noext" | sed -nE 's/^sub-[^_]+_ses-[^_]+_T2w_.*_(left|right)$/\1/p')"

    if [[ -n "${roi_label:-}" && -n "${hemi:-}" ]]; then
      target_name="${sub}_${ses}_space-T2w_desc-${roi_label}_${hemi}_roi${ext}"
    fi
  fi

  local target_path="$target_dir/$target_name"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY] cp -n \"$f\" \"$target_path\""
  else
    cp -n "$f" "$target_path"
    echo "[OK ] Copied -> $target_path"
  fi
}

# Find and process files without mapfile
found_any=0
find "$ROI_SRC" -type f \( \
  -name "sub-*_ses-*_T2w_*_left.nii" -o -name "sub-*_ses-*_T2w_*_right.nii" -o \
  -name "sub-*_ses-*_T2w_*_left.nii.gz" -o -name "sub-*_ses-*_T2w_*_right.nii.gz" \
\) -print0 | while IFS= read -r -d '' f; do
  found_any=1
  copy_one "$f"
done

# Note: due to subshell behavior in some shells, found_any may not update reliably.
# So also do a quick check:
if ! find "$ROI_SRC" -type f \( \
  -name "sub-*_ses-*_T2w_*_left.nii" -o -name "sub-*_ses-*_T2w_*_right.nii" -o \
  -name "sub-*_ses-*_T2w_*_left.nii.gz" -o -name "sub-*_ses-*_T2w_*_right.nii.gz" \
\) | head -n 1 | grep -q .; then
  echo "No matching ROI files found in: $ROI_SRC"
  exit 1
fi

echo
echo "Done."
echo "If you want to actually copy, set: DRY_RUN=0"
