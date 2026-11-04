#!/usr/bin/env bash
# tsnr_report.sh <derivatives_root> [out.csv]   every 4D run in the tree -> CSV
# tsnr_report.sh <one_bold.nii>                 one run -> CSV on stdout
#
# Temporal SNR for every 4D run x ROI mask in a preprocessed derivatives tree.
# Needs FSL on the PATH and nothing else. For the cluster, and to run the runs
# in parallel, use tsnr_report.sbatch - it calls this script in its one-run
# form, which is why that form exists.
#
#   tools/tsnr_report.sh /path/to/derivatives/spm-preproc-oap-viarest tsnr.csv
#
# Assumes only <root>/sub-*/ses-*/{func,anat}/. Image and mask names are globs,
# so a dataset that names its ROIs differently needs no edit here:
#
#   FUNC_GLOB='ua*_bold.nii*' tools/tsnr_report.sh <root>      # unsmoothed
#   MASK_GLOB='*_roi.nii*'    tools/tsnr_report.sh <root>      # other ROIs
#   HPF_SEC=0                 tools/tsnr_report.sh <root>      # skip detrend
#
# The mask column holds the file's basename, unparsed - ROI naming varies
# between studies and whatever reads the CSV can split it.
#
# COLUMNS
#   roi_tsnr    mean over mask voxels of mean/SD over time. Zeros are kept:
#               a dropout voxel is a real 0, not a missing value.
#   brain_tsnr  the same over the whole brain, so roi_tsnr has a reference.
#   spm_global  SPM's global: mean of voxels above (image mean)/8.
#   roi_over_mthresh
#               roi_mean / (0.8 * spm_global). Below 1 means the ROI sits
#               under SPM's default implicit-mask threshold and the GLM will
#               drop those voxels - the same arithmetic as maskingThreshold
#               in spm_options.m.
#
# By default the drift is removed first (high-pass at HPF_SEC, matching the
# GLM's HPF_cutoff) so the SD is the noise the GLM actually sees. Raw tSNR
# counts scanner drift as noise and reads lower.

target=${1:?usage: tsnr_report.sh <derivatives_root|one_bold.nii> [out.csv]}
out=${2:-tsnr.csv}
FUNC_GLOB=${FUNC_GLOB:-sua*_bold.nii*}
MASK_GLOB=${MASK_GLOB:-r*_mask.nii*}
HPF_SEC=${HPF_SEC:-128}

HEADER="subject,session,run,image,mask,n_vox,roi_tsnr,roi_mean,brain_tsnr,spm_global,roi_over_mthresh"

command -v fslmaths >/dev/null 2>&1 || { echo "FSL is not on the PATH." >&2; exit 1; }

tmp=$(mktemp -d) || exit 1
trap 'rm -rf "$tmp"' EXIT

# fslstats prints a TRAILING SPACE after every value, which printf %f rejects
# with "invalid number" and then writes 0. Everything goes through this.
stat1() { fslstats "$@" | awk '{ print $1+0 }'; }

# one_run <4d_bold>  ->  one CSV row per ROI mask, on stdout.
# Returns 2 if the file is not 4D, so the tree walk can skip 3D means (which
# match the same glob) while a single-file call can complain about one.
one_run() {
    local func=$1
    local nvol sesdir ses sub base run img mall thr gl btsnr trs sigma wrote=0

    nvol=$(fslval "$func" dim4 2>/dev/null | tr -d '[:space:]')
    case $nvol in ''|*[!0-9]*) return 2 ;; esac
    [ "$nvol" -gt 1 ] || return 2

    sesdir=$(dirname "$(dirname "$func")")
    ses=$(basename "$sesdir")
    sub=$(basename "$(dirname "$sesdir")")
    base=$(basename "$func")
    run=${base#*task-}; run=${run%_bold*}
    img=${base%%sub-*}; img=${img:-raw}
    printf '  %s %s %-10s %s (%s vols)\n' "$sub" "$ses" "$run" "$img" "$nvol" >&2

    fslmaths "$func" -Tmean "$tmp/mean" || return 1
    if [ "${HPF_SEC%%.*}" -gt 0 ] 2>/dev/null; then
        trs=$(fslval "$func" pixdim4 | tr -d '[:space:]')
        sigma=$(awk -v c="$HPF_SEC" -v t="$trs" 'BEGIN{ printf "%.4f", (t>0)? c/(2*t) : 0 }')
        fslmaths "$func" -bptf "$sigma" -1 "$tmp/filt" || return 1
        fslmaths "$tmp/filt" -Tstd "$tmp/sd" || return 1
    else
        fslmaths "$func" -Tstd "$tmp/sd" || return 1
    fi
    fslmaths "$tmp/mean" -div "$tmp/sd" "$tmp/tsnr" || return 1

    # SPM's global, and a brain mask built the same way
    mall=$(stat1 "$tmp/mean" -m)
    thr=$(awk -v m="$mall" 'BEGIN{ printf "%.8f", m/8 }')
    gl=$(stat1 "$tmp/mean" -l "$thr" -M)
    fslmaths "$tmp/mean" -thr "$thr" -bin "$tmp/brain"
    btsnr=$(stat1 "$tmp/tsnr" -k "$tmp/brain" -m)

    local mask mname nvox rt rm_ ratio
    for mask in "$sesdir"/anat/$MASK_GLOB; do
        mname=$(basename "$mask"); mname=${mname%%.nii*}
        nvox=$(stat1 "$mask" -l 0.5 -V)
        rt=$(stat1 "$tmp/tsnr" -k "$mask" -m)
        rm_=$(stat1 "$tmp/mean" -k "$mask" -m)
        ratio=$(awk -v r="$rm_" -v g="$gl" 'BEGIN{ printf "%.4f", (g>0)? r/(0.8*g) : 0 }')
        printf '%s,%s,%s,%s,%s,%s,%.2f,%.2f,%.2f,%.2f,%s\n' \
            "$sub" "$ses" "$run" "$img" "$mname" "$nvox" \
            "$rt" "$rm_" "$btsnr" "$gl" "$ratio"
        wrote=1
    done
    if [ "$wrote" -eq 0 ]; then       # no masks here: keep the brain numbers
        printf '%s,%s,%s,%s,,,,,%.2f,%.2f,\n' \
            "$sub" "$ses" "$run" "$img" "$btsnr" "$gl"
    fi
    rm -f "$tmp"/mean.* "$tmp"/sd.* "$tmp"/filt.* "$tmp"/tsnr.* "$tmp"/brain.*
}

# ---- one run, for tsnr_report.sbatch and for looking at a single file ----
if [ -f "$target" ]; then
    nv=$(fslval "$target" dim4 2>/dev/null | tr -d '[:space:]')
    case $nv in ''|*[!0-9]*) nv=0 ;; esac
    [ "$nv" -gt 1 ] || { echo "Not a 4D image: $target" >&2; exit 1; }
    echo "$HEADER"
    one_run "$target"
    exit $?
fi

# ---- the whole tree ----
[ -d "$target" ] || { echo "No such file or directory: $target" >&2; exit 1; }
shopt -s nullglob
echo "$HEADER" > "$out"
nrun=0
for func in "$target"/sub-*/ses-*/func/$FUNC_GLOB; do
    if one_run "$func" >> "$out"; then
        nrun=$((nrun + 1))
    fi
done
echo "$nrun runs -> $out" >&2
