#!/usr/bin/env bash
# =============================================================================
#  Remove pipeline output that earlier versions wrote into BIDS_ROOT
# =============================================================================
#  Before the derivatives layout, the preprocessing ran directly on the BIDS
#  directory: SPM wrote u*_bold.nii, rp_*.txt, voxel displacement maps, topup
#  results and resliced ROI masks next to the raw data, and every .nii.gz that
#  SPM had to read was left behind decompressed as well.
#
#  A dataset processed with one of those versions still holds all of it. It is
#  dead weight — nothing reads it any more — and while it is there you cannot
#  tell by looking which files came off the scanner.
#
#  This lists what it would remove and does nothing until you pass --execute.
#
#  WHAT IT REMOVES
#    func/  <u>sub-*_bold*.nii  <su>sub-*_bold*.nii  *_bold_dc.nii
#           rp_sub-*.txt  mean*sub-*.nii  vdm_*  vdm5_*  topup_*
#           fieldmap_hz.nii  preproc_provenance.json  preproc_stages.log
#    anat/  <r>sub-*.nii — the masks the coregistration resliced
#    any/   a .nii that has a .nii.gz of the same name beside it
#
#  Every pattern is anchored on the SPM prefix followed by "sub-", because a
#  raw BIDS filename starts with "sub-" itself: a pattern like su*_bold*.nii
#  matches sub-01_ses-01_task-run1_bold.nii, which is raw data.
#
#  WHAT IT NEVER TOUCHES
#    Any file named like raw data — sub-*_bold.nii, _T1w, _T2w, _epi,
#    _phasediff, _magnitude*, _fieldmap, _sbref, _events.tsv, a JSON sidecar,
#    or a mask that is not a resliced copy — enforced as a filter over
#    whatever the patterns matched, not just by the patterns themselves.
#
#    The one exception is the decompressed twins: a .nii that has a .nii.gz of
#    the same name beside it is a copy SPM made, and the .nii.gz original is
#    left alone. Pass --keep-uncompressed to skip those too.
#
#  Usage:
#    bash clean_bids_leftovers.sh              # list what would go
#    bash clean_bids_leftovers.sh --execute    # actually delete
#    bash clean_bids_leftovers.sh --sub sub-01s13 --ses ses-01
#
#  Options:
#    --execute        Delete. Without it nothing is written.
#    --sub SUB        Restrict to one subject...
#    --ses SES        ...and one session.
#    --keep-uncompressed  Do not remove .nii files that have a .nii.gz twin.
#    --config FILE    Config file (default: <script dir>/pipeline_config.cfg).
#    -h, --help       Show this help.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG_FILE="${SELF_DIR}/pipeline_config.cfg"
SUB_ARG=""
SES_ARG=""
EXECUTE=false
KEEP_UNCOMPRESSED=false

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

show_help() {
    awk 'NR>1 { if ($0 !~ /^#/) exit; print }' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --execute)           EXECUTE=true; shift ;;
        --sub)               SUB_ARG="$2"; shift 2 ;;
        --ses)               SES_ARG="$2"; shift 2 ;;
        --keep-uncompressed) KEEP_UNCOMPRESSED=true; shift ;;
        --config)            CONFIG_FILE="$2"; shift 2 ;;
        -h|--help)           show_help; exit 0 ;;
        *)
            printf "${RED}ERROR: Unknown option: %s${NC}\n" "$1" >&2
            printf "       Run with --help to see the accepted options.\n" >&2
            exit 1 ;;
    esac
done

if [[ ! -f "$CONFIG_FILE" ]]; then
    printf "${RED}ERROR: Config file not found: %s${NC}\n" "$CONFIG_FILE" >&2
    exit 1
fi
# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${BIDS_ROOT:?ERROR: BIDS_ROOT not set in pipeline_config.cfg}"

if [[ ! -d "$BIDS_ROOT" ]]; then
    printf "${RED}ERROR: BIDS_ROOT not found: %s${NC}\n" "$BIDS_ROOT" >&2
    exit 1
fi

SCOPE="$BIDS_ROOT"
if [[ -n "$SUB_ARG" ]]; then
    SCOPE="${BIDS_ROOT}/${SUB_ARG}"
    [[ -n "$SES_ARG" ]] && SCOPE="${SCOPE}/${SES_ARG}"
    if [[ ! -d "$SCOPE" ]]; then
        printf "${RED}ERROR: not found: %s${NC}\n" "$SCOPE" >&2
        exit 1
    fi
fi

printf "${BLUE}=============================================${NC}\n"
printf "${BLUE}  Clean pipeline leftovers out of BIDS_ROOT${NC}\n"
printf "${BLUE}=============================================${NC}\n"
printf "  Scanning: %s\n" "$SCOPE"
if [[ "$EXECUTE" == "true" ]]; then
    printf "  ${GREEN}Mode:     EXECUTE — files will be deleted${NC}\n"
else
    printf "  ${YELLOW}Mode:     DRY RUN — nothing will be deleted${NC}\n"
fi
printf "\n"

TMP="$(mktemp)"
trap 'rm -f "$TMP" "${TMP}.gz" "${TMP}.keep"' EXIT

# ---- SPM output, by name ----------------------------------------------------
# The prefixes come from the config, and every one of them is followed by
# "sub-". Without that anchor "su*_bold*.nii" matches the raw
# sub-01_ses-01_task-run1_bold.nii, because a BIDS name starts with "sub-".
U_PREFIX="${RESLICE_PREFIX:-u}"
S_PREFIX="${SMOOTH_PREFIX:-s}"
R_PREFIX="${ROI_PREFIX:-r}"

find "$SCOPE" -type d -name func -print0 2>/dev/null | while IFS= read -r -d '' d; do
    find "$d" -maxdepth 1 -type f \( \
        -name "${U_PREFIX}sub-*_bold*.nii"            -o \
        -name "${S_PREFIX}${U_PREFIX}sub-*_bold*.nii" -o \
        -name "*_bold_dc.nii"                         -o \
        -name "rp_sub-*.txt"                          -o \
        -name "mean${U_PREFIX}sub-*.nii"              -o \
        -name "meansub-*.nii"                         -o \
        -name "meanu_session.nii"                     -o \
        -name "vdm_*"    -o -name "vdm5_*"            -o \
        -name "topup_*"  -o -name "fieldmap_hz.nii"   -o \
        -name "preproc_provenance.json" -o -name "preproc_stages.log" \
    \) 2>/dev/null
done >> "$TMP"

# Masks the coregistration resliced: the ROI prefix, then the mask's own name.
if [[ -n "$R_PREFIX" ]]; then
    find "$SCOPE" -type d -name anat -print0 2>/dev/null | while IFS= read -r -d '' d; do
        find "$d" -maxdepth 1 -type f -name "${R_PREFIX}sub-*.nii" 2>/dev/null
    done >> "$TMP"
fi

# Stale VDMs the FieldMap toolbox dropped into fmap/
find "$SCOPE" -type d -name fmap -print0 2>/dev/null | while IFS= read -r -d '' d; do
    find "$d" -maxdepth 1 -type f \( -name "vdm_*" -o -name "vdm5_*" \) 2>/dev/null
done >> "$TMP"

# ---- Guard the pattern matches -----------------------------------------------
# Whatever the patterns above matched, a file named like raw data does not get
# deleted. This is the last line of defence between a mistyped prefix and the
# scanner's output.
grep -vE '/sub-[^/]*_(bold|sbref|T1w|T2w|T2star|FLAIR|PD|epi|phasediff|phase[12]|magnitude[12]?|fieldmap|dwi|events)\.(nii|nii\.gz|tsv|json)$' \
     "$TMP" > "${TMP}.keep" 2>/dev/null || true
mv "${TMP}.keep" "$TMP"

# ---- Decompressed twins ------------------------------------------------------
# ensure_nii used to gunzip in place, so the raw .nii.gz and a .nii of the same
# name both exist. BIDS never stores a scan both ways, so a .nii sitting beside
# a .nii.gz of the same name is a decompressed copy and the original stays.
# Added after the guard on purpose: these twins are named exactly like the raw
# files they came from, which is what the guard is there to protect.
if [[ "$KEEP_UNCOMPRESSED" != "true" ]]; then
    find "$SCOPE" -type f -name "*.nii.gz" 2>/dev/null > "${TMP}.gz" || true
    while IFS= read -r gz; do
        [[ -z "$gz" ]] && continue
        plain="${gz%.gz}"
        [[ -f "$plain" ]] && printf '%s\n' "$plain"
    done < "${TMP}.gz" >> "$TMP"
fi

sort -u "$TMP" -o "$TMP"

COUNT=$(grep -cve '^[[:space:]]*$' "$TMP" || true)
if [[ "$COUNT" -eq 0 ]]; then
    printf "${GREEN}Nothing to clean — BIDS_ROOT holds no pipeline output.${NC}\n\n"
    exit 0
fi

BYTES=$(while IFS= read -r f; do
            [[ -f "$f" ]] && wc -c < "$f"
        done < "$TMP" | awk '{s+=$1} END {print s+0}')
HUMAN=$(awk -v b="$BYTES" 'BEGIN {
    split("B KB MB GB TB", u, " "); i = 1
    while (b >= 1024 && i < 5) { b /= 1024; i++ }
    printf "%.1f %s", b, u[i]
}')

printf "%s file(s), %s\n\n" "$COUNT" "$HUMAN"
sed "s|^${BIDS_ROOT}/||" "$TMP" | head -40
[[ "$COUNT" -gt 40 ]] && printf "  ... and %s more\n" "$((COUNT - 40))"
printf "\n"

if [[ "$EXECUTE" != "true" ]]; then
    printf "${YELLOW}Dry run. To delete these, run:${NC}\n"
    printf "  bash %s --execute\n\n" "$0"
    exit 0
fi

DELETED=0
FAILED=0
while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    if rm -f "$f" 2>/dev/null; then
        DELETED=$((DELETED + 1))
    else
        FAILED=$((FAILED + 1))
        printf "${RED}could not remove: %s${NC}\n" "$f" >&2
    fi
done < "$TMP"

printf "${GREEN}Deleted %s file(s)${NC}" "$DELETED"
[[ "$FAILED" -gt 0 ]] && printf ", ${RED}%s failed${NC}" "$FAILED"
printf "\n\n"
printf "From here on the preprocessing writes to the derivatives instead, so\n"
printf "BIDS_ROOT stays as it is.\n\n"

[[ "$FAILED" -gt 0 ]] && exit 1
exit 0
