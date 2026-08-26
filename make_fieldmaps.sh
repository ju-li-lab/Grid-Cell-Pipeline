#!/usr/bin/env bash
# =============================================================================
#  Create B0 fieldmaps from a structural T1w + phasediff (FSL)
# =============================================================================
#  For sessions where the scanner produced a phase-difference image but no
#  usable magnitude image, this builds the magnitude that
#  fsl_prepare_fieldmap needs out of the T1w:
#
#    1. bet      — brain-extract the T1w
#    2. flirt    — resample the T1w into the phasediff grid (sform/qform based)
#    3. flirt    — resample the brain mask the same way (nearest neighbour)
#    4. fslmaths — apply the mask, giving a brain-only pseudo-magnitude
#    5. fsl_prepare_fieldmap — unwrap the phase and write a fieldmap in rad/s
#
#  Outputs land in the session's fmap/ directory with BIDS "Case 3" names, so
#  the preprocessing picks them up without any renaming:
#
#    fmap/<sub>_<ses>[_run-N]_fieldmap.nii[.gz]    the fieldmap, in rad/s
#    fmap/<sub>_<ses>[_run-N]_magnitude.nii[.gz]   the brain-only magnitude
#    fmap/<sub>_<ses>[_run-N]_fieldmap.json        units, IntendedFor, provenance
#
#  Then set PREPROC_MODE=precalc_fieldmap in pipeline_config.cfg and run the
#  preprocessing as usual.
#
#  Usage:
#    bash make_fieldmaps.sh                     # every session in subses_list.txt
#    bash make_fieldmaps.sh --dry-run           # show what would happen
#    bash make_fieldmaps.sh --sub sub-a01 --ses ses-02
#    bash make_fieldmaps.sh --force             # overwrite existing fieldmaps
#
#  Options:
#    --sub SUB          Process one subject...
#    --ses SES          ...and this session (both required together).
#    --list FILE        Subject-session list (default: <script dir>/subses_list.txt).
#    --config FILE      Config file (default: <script dir>/pipeline_config.cfg).
#    --magnitude-source structural|auto
#                       structural = always build the magnitude from the T1w
#                                    (default; keeps every session identical).
#                       auto       = use a real *_magnitude1 image when the
#                                    session has one, T1w otherwise.
#    --delta-te MS      Echo time difference in ms. Default: read per session
#                       from the phasediff JSON (EchoTime2 - EchoTime1).
#    --force            Rebuild fieldmaps that already exist.
#    --keep-work        Keep the intermediate bet/flirt files for inspection.
#    --dry-run          Print the commands without running them.
#    -h, --help         Show this help.
# =============================================================================

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG_FILE="${SELF_DIR}/pipeline_config.cfg"
SUB_ARG=""
SES_ARG=""
LIST=""
MAG_SOURCE=""
DELTA_TE_ARG=""
FORCE=false
KEEP_WORK=false
DRY_RUN_ARG=false

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

show_help() {
    awk 'NR>1 { if ($0 !~ /^#/) exit; print }' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

info()  { printf "${BLUE}ℹ${NC} %b\n" "$*"; }
ok()    { printf "${GREEN}✓${NC} %b\n" "$*"; }
warn()  { printf "${YELLOW}⚠${NC} %b\n" "$*"; }
err()   { printf "${RED}✗${NC} %b\n" "$*" >&2; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sub)              SUB_ARG="$2";       shift 2 ;;
        --ses)              SES_ARG="$2";       shift 2 ;;
        --list)             LIST="$2";          shift 2 ;;
        --config)           CONFIG_FILE="$2";   shift 2 ;;
        --magnitude-source) MAG_SOURCE="$2";    shift 2 ;;
        --delta-te)         DELTA_TE_ARG="$2";  shift 2 ;;
        --force)            FORCE=true;         shift ;;
        --keep-work)        KEEP_WORK=true;     shift ;;
        --dry-run)          DRY_RUN_ARG=true;   shift ;;
        -h|--help)          show_help; exit 0 ;;
        *)
            err "Unknown option: $1"
            printf "       Run with --help to see the accepted options.\n" >&2
            exit 1 ;;
    esac
done

if [[ ! -f "$CONFIG_FILE" ]]; then
    err "Config file not found: $CONFIG_FILE"
    exit 1
fi
source "$CONFIG_FILE"

: "${BIDS_ROOT:?ERROR: BIDS_ROOT not set in pipeline_config.cfg}"

[[ -z "$LIST" ]]       && LIST="${SCRIPT_DIR:-$SELF_DIR}/subses_list.txt"
[[ -z "$MAG_SOURCE" ]] && MAG_SOURCE="${FIELDMAP_MAGNITUDE_SOURCE:-structural}"
[[ -z "$DELTA_TE_ARG" ]] && DELTA_TE_ARG="${FIELDMAP_DELTA_TE:-}"

SCANNER="${FIELDMAP_SCANNER:-SIEMENS}"
BET_F="${FIELDMAP_BET_F:-0.5}"

# --dry-run on the command line, or DRY_RUN=true in the config
DRY_RUN_EFF=false
if [[ "$DRY_RUN_ARG" == "true" || "${DRY_RUN:-false}" == "true" ]]; then
    DRY_RUN_EFF=true
fi

if [[ ! "$MAG_SOURCE" =~ ^(structural|auto)$ ]]; then
    err "--magnitude-source must be 'structural' or 'auto' (got '${MAG_SOURCE}')"
    exit 1
fi

if [[ -n "$SUB_ARG" || -n "$SES_ARG" ]]; then
    if [[ -z "$SUB_ARG" || -z "$SES_ARG" ]]; then
        err "--sub and --ses must be given together."
        exit 1
    fi
fi

# ---- FSL ----
if ! command -v fsl_prepare_fieldmap >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load "${FSL_MODULE:-fsl}" 2>/dev/null || true
    fi
fi
if [[ "$DRY_RUN_EFF" != "true" ]]; then
    for tool in bet flirt fslmaths fsl_prepare_fieldmap; do
        command -v "$tool" >/dev/null 2>&1 || {
            err "FSL tool '$tool' not found on PATH."
            printf "       Load FSL first (module load %s), or run with --dry-run.\n" "${FSL_MODULE:-fsl}" >&2
            exit 1
        }
    done

    # bet and fsl_prepare_fieldmap are wrapper scripts that call their helpers
    # through $FSLDIR. Without it they fail with a bare "/bin/remove_ext: No
    # such file or directory", which is not obvious to diagnose.
    if [[ -z "${FSLDIR:-}" ]]; then
        _bet_path="$(command -v bet)"
        _guess="$(cd "$(dirname "$_bet_path")/.." && pwd)"
        if [[ -x "${_guess}/bin/remove_ext" ]]; then
            export FSLDIR="$_guess"
            warn "FSLDIR was not set — using ${FSLDIR}"
        else
            err "FSLDIR is not set and could not be guessed from $(dirname "$_bet_path")."
            printf "       Source FSL's environment first, e.g.\n" >&2
            printf "         export FSLDIR=/path/to/fsl && . \$FSLDIR/etc/fslconf/fsl.sh\n" >&2
            printf "       On the cluster: module load %s\n" "${FSL_MODULE:-fsl}" >&2
            exit 1
        fi
    fi
fi

# Match the compression of the dataset rather than imposing FSL's default, so
# the new files sit alongside the existing ones in the same format.
export FSLOUTPUTTYPE="${FSLOUTPUTTYPE:-NIFTI_GZ}"

run_cmd() {
    printf "    %s\n" "$*"
    [[ "$DRY_RUN_EFF" == "true" ]] && return 0

    local out rc
    out="$("$@" 2>&1)" && return 0
    rc=$?
    err "Command failed (exit ${rc}): $*"
    [[ -n "$out" ]] && printf "      %s\n" "$out" | head -20 >&2
    return "$rc"
}

# Echo the first path that exists, from the arguments given
first_existing() {
    local f
    for f in "$@"; do
        [[ -f "$f" ]] && { printf '%s' "$f"; return 0; }
    done
    return 1
}

# ---- Build the list of subject-sessions to process ----
declare -a TARGETS=()
if [[ -n "$SUB_ARG" ]]; then
    TARGETS+=("${SUB_ARG} ${SES_ARG}")
else
    if [[ ! -f "$LIST" ]]; then
        err "Subject list not found: $LIST"
        printf "       Run step0_make_subses_list.sh first, or pass --sub/--ses.\n" >&2
        exit 1
    fi
    while read -r s ss _rest; do
        [[ -z "${s:-}" ]] && continue
        TARGETS+=("${s} ${ss}")
    done < "$LIST"
fi

printf "\n"
printf "${BLUE}=============================================${NC}\n"
printf "${BLUE}  Create fieldmaps from structural + phasediff${NC}\n"
printf "${BLUE}=============================================${NC}\n"
printf "  BIDS root:        %s\n" "$BIDS_ROOT"
printf "  Sessions:         %d\n" "${#TARGETS[@]}"
printf "  Scanner:          %s\n" "$SCANNER"
printf "  Magnitude source: %s\n" "$MAG_SOURCE"
printf "  BET -f:           %s\n" "$BET_F"
if [[ -n "$DELTA_TE_ARG" ]]; then
    printf "  Delta TE:         %s ms (fixed)\n" "$DELTA_TE_ARG"
else
    printf "  Delta TE:         from each phasediff JSON (EchoTime2 - EchoTime1)\n"
fi
[[ "$DRY_RUN_EFF" == "true" ]] && printf "  ${YELLOW}DRY RUN — nothing will be written${NC}\n"
printf "\n"

N_MADE=0; N_SKIPPED=0; N_FAILED=0

for entry in "${TARGETS[@]}"; do
    read -r SUB SES <<< "$entry"
    SESDIR="${BIDS_ROOT}/${SUB}/${SES}"
    FMAPDIR="${SESDIR}/fmap"
    ANATDIR="${SESDIR}/anat"

    printf "${BLUE}── %s / %s ──${NC}\n" "$SUB" "$SES"

    if [[ ! -d "$FMAPDIR" ]]; then
        warn "no fmap/ directory — skipping"
        N_SKIPPED=$((N_SKIPPED+1)); printf "\n"; continue
    fi

    # ---- Every phasediff in this session (may carry run-N entities) ----
    declare -a PHASEDIFFS=()
    while IFS= read -r f; do
        [[ -n "$f" ]] && PHASEDIFFS+=("$f")
    done < <(find "$FMAPDIR" -maxdepth 1 \
                  \( -name "${SUB}_${SES}*phasediff.nii" -o -name "${SUB}_${SES}*phasediff.nii.gz" \) \
                  2>/dev/null | sort)

    if [[ ${#PHASEDIFFS[@]} -eq 0 ]]; then
        warn "no ${SUB}_${SES}*phasediff.nii[.gz] in fmap/ — skipping"
        N_SKIPPED=$((N_SKIPPED+1)); printf "\n"; continue
    fi

    for PHASE in "${PHASEDIFFS[@]}"; do
        # Entity stem: the filename minus the _phasediff suffix and extension.
        # sub-a01_ses-04_run-1_phasediff.nii -> sub-a01_ses-04_run-1
        base="$(basename "$PHASE")"
        stem="${base%_phasediff.nii.gz}"; stem="${stem%_phasediff.nii}"

        # Keep the dataset's own compression for the new files
        if [[ "$PHASE" == *.nii.gz ]]; then
            export FSLOUTPUTTYPE=NIFTI_GZ; OUT_EXT=".nii.gz"
        else
            export FSLOUTPUTTYPE=NIFTI;    OUT_EXT=".nii"
        fi

        OUT_FMAP="${FMAPDIR}/${stem}_fieldmap${OUT_EXT}"
        OUT_MAG="${FMAPDIR}/${stem}_magnitude${OUT_EXT}"
        OUT_JSON="${FMAPDIR}/${stem}_fieldmap.json"

        printf "  ${BLUE}%s${NC}\n" "$stem"

        existing="$(first_existing "${FMAPDIR}/${stem}_fieldmap.nii" "${FMAPDIR}/${stem}_fieldmap.nii.gz" || true)"
        if [[ -n "$existing" && "$FORCE" != "true" ]]; then
            info "fieldmap already exists — skipping (use --force to rebuild)"
            printf "      %s\n" "$existing"
            N_SKIPPED=$((N_SKIPPED+1)); continue
        fi

        # ---- Delta TE ----
        PHASE_JSON="${FMAPDIR}/${stem}_phasediff.json"
        if [[ -n "$DELTA_TE_ARG" ]]; then
            DELTA_TE="$DELTA_TE_ARG"
            info "delta TE: ${DELTA_TE} ms (from config/--delta-te)"
        elif [[ -f "$PHASE_JSON" ]]; then
            DELTA_TE=$(python3 -c "
import json,sys
try:
    j=json.load(open('${PHASE_JSON}'))
    te1=j.get('EchoTime1'); te2=j.get('EchoTime2')
    if te1 is None or te2 is None: sys.exit(1)
    print(round((te2-te1)*1000, 6))
except Exception:
    sys.exit(1)
" 2>/dev/null) || DELTA_TE=""
            if [[ -z "$DELTA_TE" ]]; then
                err "could not read EchoTime1/EchoTime2 from ${PHASE_JSON}"
                printf "      Set FIELDMAP_DELTA_TE in the config or pass --delta-te.\n" >&2
                N_FAILED=$((N_FAILED+1)); continue
            fi
            info "delta TE: ${DELTA_TE} ms (from $(basename "$PHASE_JSON"))"
        else
            err "no ${stem}_phasediff.json and no FIELDMAP_DELTA_TE set"
            printf "      Pass --delta-te MS, or set FIELDMAP_DELTA_TE in the config.\n" >&2
            N_FAILED=$((N_FAILED+1)); continue
        fi

        # ---- Work directory ----
        WORK="${FMAPDIR}/.work_${stem}"
        if [[ "$DRY_RUN_EFF" != "true" ]]; then
            rm -rf "$WORK"; mkdir -p "$WORK"
        fi

        MAG_KIND=""
        MAG_BRAIN="${WORK}/magnitude_brain"

        # ---- Real magnitude, if we are allowed to use one ----
        REAL_MAG=""
        if [[ "$MAG_SOURCE" == "auto" ]]; then
            REAL_MAG="$(first_existing \
                "${FMAPDIR}/${stem}_magnitude1.nii"    "${FMAPDIR}/${stem}_magnitude1.nii.gz" \
                "${FMAPDIR}/${SUB}_${SES}_magnitude1.nii" "${FMAPDIR}/${SUB}_${SES}_magnitude1.nii.gz" \
                || true)"
        fi

        if [[ -n "$REAL_MAG" ]]; then
            # The session has a real GRE magnitude — brain-extract it directly
            MAG_KIND="magnitude1"
            info "magnitude: $(basename "$REAL_MAG") (real GRE magnitude)"
            run_cmd bet "$REAL_MAG" "$MAG_BRAIN" -f "$BET_F" || { N_FAILED=$((N_FAILED+1)); continue; }
        else
            # ---- Build a pseudo-magnitude from the T1w ----
            MAG_KIND="T1w"

            # Prefer a T1w carrying the same run entity as this phasediff
            RUN_ENT=""
            [[ "$stem" =~ (_run-[A-Za-z0-9]+) ]] && RUN_ENT="${BASH_REMATCH[1]}"

            T1W=""
            if [[ -n "$RUN_ENT" ]]; then
                T1W="$(first_existing \
                    "${ANATDIR}/${SUB}_${SES}${RUN_ENT}_T1w.nii" \
                    "${ANATDIR}/${SUB}_${SES}${RUN_ENT}_T1w.nii.gz" || true)"
            fi
            if [[ -z "$T1W" ]]; then
                # Otherwise the last T1w in sorted order, matching how the
                # pipeline picks a T2w when several runs exist.
                T1W=$(find "$ANATDIR" -maxdepth 1 \( -name "${SUB}_${SES}*_T1w.nii" -o -name "${SUB}_${SES}*_T1w.nii.gz" \) 2>/dev/null | sort | tail -1)
            fi

            if [[ -z "$T1W" ]]; then
                err "no ${SUB}_${SES}*_T1w.nii[.gz] in ${ANATDIR}"
                N_FAILED=$((N_FAILED+1)); continue
            fi
            info "magnitude: built from $(basename "$T1W")"

            T1_BRAIN="${WORK}/T1w_brain"
            MAG_RAW="${WORK}/magnitude"
            MAG_MASK="${WORK}/magnitude_brain_mask"

            # 1. brain-extract the T1w (-m also writes <out>_mask)
            run_cmd bet "$T1W" "$T1_BRAIN" -m -f "$BET_F" || { N_FAILED=$((N_FAILED+1)); continue; }

            # 2. resample the T1w into the phasediff grid using the stored
            #    scanner coordinates (-usesqform), no registration
            run_cmd flirt -in "$T1W" -ref "$PHASE" -applyxfm -usesqform -out "$MAG_RAW" \
                || { N_FAILED=$((N_FAILED+1)); continue; }

            # 3. same for the brain mask, nearest neighbour to keep it binary
            T1_MASK="$(first_existing "${T1_BRAIN}_mask.nii" "${T1_BRAIN}_mask.nii.gz" || echo "${T1_BRAIN}_mask")"
            run_cmd flirt -in "$T1_MASK" -ref "$PHASE" -applyxfm -usesqform \
                    -out "$MAG_MASK" -interp nearestneighbour \
                || { N_FAILED=$((N_FAILED+1)); continue; }

            # 4. mask the resampled T1w -> brain-only pseudo-magnitude
            run_cmd fslmaths "$MAG_RAW" -mas "$MAG_MASK" "$MAG_BRAIN" \
                || { N_FAILED=$((N_FAILED+1)); continue; }
        fi

        # ---- 5. Unwrap the phase and write the fieldmap (rad/s) ----
        run_cmd fsl_prepare_fieldmap "$SCANNER" "$PHASE" "$MAG_BRAIN" "${WORK}/fieldmap" "$DELTA_TE" \
            || { N_FAILED=$((N_FAILED+1)); continue; }

        # ---- Move the results into fmap/ under BIDS names ----
        if [[ "$DRY_RUN_EFF" == "true" ]]; then
            printf "    would write %s\n" "$OUT_FMAP"
            printf "    would write %s\n" "$OUT_MAG"
            printf "    would write %s\n" "$OUT_JSON"
            N_MADE=$((N_MADE+1))
            continue
        fi

        FMAP_SRC="$(first_existing "${WORK}/fieldmap.nii" "${WORK}/fieldmap.nii.gz" || true)"
        MAGB_SRC="$(first_existing "${MAG_BRAIN}.nii" "${MAG_BRAIN}.nii.gz" || true)"
        if [[ -z "$FMAP_SRC" || -z "$MAGB_SRC" ]]; then
            err "fsl_prepare_fieldmap produced no output in ${WORK}"
            N_FAILED=$((N_FAILED+1)); continue
        fi

        # Remove any previous output in the other compression, so the session
        # never ends up with both a .nii and a .nii.gz fieldmap.
        rm -f "${FMAPDIR}/${stem}_fieldmap.nii" "${FMAPDIR}/${stem}_fieldmap.nii.gz" \
              "${FMAPDIR}/${stem}_magnitude.nii" "${FMAPDIR}/${stem}_magnitude.nii.gz"

        mv "$FMAP_SRC" "$OUT_FMAP"
        mv "$MAGB_SRC" "$OUT_MAG"

        # ---- Sidecar: units, what the EPIs are, how it was made ----
        INTENDED=""
        if [[ -n "${TASKS:-}" ]]; then
            IFS=',' read -ra _tasks <<< "$TASKS"
            for t in "${_tasks[@]}"; do
                t="${t// /}"
                [[ -z "$t" ]] && continue
                for bold in "${SESDIR}/func/${SUB}_${SES}"*"task-${t}"*_bold.nii "${SESDIR}/func/${SUB}_${SES}"*"task-${t}"*_bold.nii.gz; do
                    [[ -f "$bold" ]] || continue
                    INTENDED="${INTENDED}${INTENDED:+, }\"${SES}/func/$(basename "$bold")\""
                done
            done
        fi

        cat > "$OUT_JSON" <<JSONEOF
{
  "Units": "rad/s",
  "IntendedFor": [${INTENDED}],
  "EchoTime1": $(python3 -c "
import json
try: print(json.load(open('${PHASE_JSON}')).get('EchoTime1','null'))
except Exception: print('null')" 2>/dev/null || echo null),
  "EchoTime2": $(python3 -c "
import json
try: print(json.load(open('${PHASE_JSON}')).get('EchoTime2','null'))
except Exception: print('null')" 2>/dev/null || echo null),
  "DeltaEchoTimeMS": ${DELTA_TE},
  "GeneratedBy": {
    "Name": "make_fieldmaps.sh",
    "Description": "fsl_prepare_fieldmap on the phasediff, using a magnitude derived from the ${MAG_KIND}",
    "MagnitudeSource": "${MAG_KIND}",
    "Scanner": "${SCANNER}",
    "BetFractionalIntensity": ${BET_F},
    "Date": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  }
}
JSONEOF

        ok "wrote $(basename "$OUT_FMAP"), $(basename "$OUT_MAG"), $(basename "$OUT_JSON")"
        N_MADE=$((N_MADE+1))

        if [[ "$KEEP_WORK" == "true" ]]; then
            info "intermediates kept in ${WORK}"
        else
            rm -rf "$WORK"
        fi
    done
    printf "\n"
done

printf "${BLUE}=============================================${NC}\n"
printf "  Created: ${GREEN}%d${NC}   Skipped: ${YELLOW}%d${NC}   Failed: ${RED}%d${NC}\n" \
    "$N_MADE" "$N_SKIPPED" "$N_FAILED"
printf "${BLUE}=============================================${NC}\n"

if [[ "$N_MADE" -gt 0 && "$DRY_RUN_EFF" != "true" ]]; then
    printf "\nNext: set ${BLUE}PREPROC_MODE=precalc_fieldmap${NC} in pipeline_config.cfg,\n"
    printf "then run the preprocessing (menu option 2, or bash submit_preproc.sh).\n"
fi

[[ "$N_FAILED" -gt 0 ]] && exit 1
exit 0
