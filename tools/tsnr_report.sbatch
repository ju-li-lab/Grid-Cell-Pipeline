#!/usr/bin/env bash
# =============================================================================
#  tsnr_report.sbatch - tSNR for a whole derivatives tree, runs in parallel
#
#    mkdir -p logs                                   # once; SLURM needs it
#    sbatch tools/tsnr_report.sbatch <derivatives_root> [out.csv]
#    sbatch -c 16 tools/tsnr_report.sbatch <root>     # 16 runs at a time
#
#  Each run is independent - one mean, one SD, a few fslstats - so there is
#  nothing to coordinate. This is one job on one node running N runs at a
#  time, N = --cpus-per-task. 54 runs take ~22 min one at a time and about
#  22/N min here; at the default 8 that is under three minutes.
#
#  A job array would spread over several nodes and finish sooner still, but it
#  would need a second job to merge the parts. One node is quick enough that
#  the extra machinery is not worth it.
#
#  FUNC_GLOB / MASK_GLOB / HPF_SEC work exactly as in tsnr_report.sh:
#    sbatch --export=ALL,FUNC_GLOB='ua*_bold.nii*' \
#           tools/tsnr_report.sbatch <root> tsnr_unsmoothed.csv
# =============================================================================
#SBATCH --job-name=tsnr
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

ROOT=${1:?usage: sbatch tools/tsnr_report.sbatch <derivatives_root> [out.csv]}

# sbatch copies this script into the SLURM spool, so its own path says nothing
# about where the repo is. SLURM_SUBMIT_DIR is where sbatch was run from.
BOOT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
OUT=${2:-${BOOT_DIR}/tsnr.csv}

CONFIG_FILE="${PIPELINE_CONFIG:-${BOOT_DIR}/config.sh}"
# shellcheck disable=SC1090
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE"          # for FSL_MODULE

if [[ -f "${BOOT_DIR}/bids_common.sh" ]]; then
    # shellcheck disable=SC1090
    source "${BOOT_DIR}/bids_common.sh"
    load_modules "${FSL_MODULE:-fsl}"
elif command -v module >/dev/null 2>&1; then
    module purge
    module load "${FSL_MODULE:-fsl}"
else
    echo "NOTE: no 'module' command available - using fsl from PATH."
fi

command -v fslmaths >/dev/null 2>&1 || {
    echo "FSL is not available. Check FSL_MODULE in config.sh (now: ${FSL_MODULE:-unset})." >&2
    exit 1
}

TSNR="${BOOT_DIR}/tools/tsnr_report.sh"
[[ -x "$TSNR" ]] || { echo "Not executable: $TSNR" >&2; exit 1; }
[[ -d "$ROOT" ]]  || { echo "No such directory: $ROOT" >&2; exit 1; }

NCPU=${SLURM_CPUS_PER_TASK:-4}
export OMP_NUM_THREADS=1          # N fslmaths processes, one thread each

shopt -s nullglob
runs=( "$ROOT"/sub-*/ses-*/func/${FUNC_GLOB:-sua*_bold.nii*} )
(( ${#runs[@]} )) || { echo "No images matched under $ROOT" >&2; exit 1; }

echo "root     : $ROOT"
echo "images   : ${#runs[@]} matched ${FUNC_GLOB:-sua*_bold.nii*} (3D ones are skipped)"
echo "parallel : $NCPU at a time"
echo "out      : $OUT"
echo

PARTS=$(mktemp -d)
trap 'rm -rf "$PARTS"' EXIT
export PARTS TSNR
export FUNC_GLOB="${FUNC_GLOB:-sua*_bold.nii*}"
export MASK_GLOB="${MASK_GLOB:-r*_mask.nii*}"
export HPF_SEC="${HPF_SEC:-128}"

# One CSV fragment per run, named after the image, so the concatenation below
# is in a stable order and two runs can never interleave inside a line.
worker() { "$TSNR" "$1" > "${PARTS}/$(basename "$1").csv"; }
export -f worker

printf '%s\0' "${runs[@]}" \
    | xargs -0 -P "$NCPU" -I{} bash -c 'worker "$@"' _ {} \
    || echo "WARNING: at least one run failed - see above." >&2

mkdir -p "$(dirname "$OUT")"
awk 'NR==1 || !/^subject,/' "${PARTS}"/*.csv > "$OUT"

rows=$(( $(wc -l < "$OUT") - 1 ))
echo
echo "$rows rows from ${#runs[@]} matched images -> $OUT"
