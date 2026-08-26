#!/usr/bin/env bash
# =============================================================================
#  SPM Preprocessing — standalone stage submitter
# =============================================================================
#  Runs the SPM preprocessing (or just some of its stages) without going
#  through the interactive menu in run_pipeline.sh. SLURM resources, paths and
#  defaults all come from pipeline_config.cfg.
#
#  Common uses:
#    # everything, whole subject list (same as menu option 2)
#    bash submit_preproc.sh
#
#    # you estimated the fieldmaps yourself — do the rest on the cluster
#    bash submit_preproc.sh --stages vdm,realign,coreg,smooth
#
#    # only re-smooth, e.g. after changing SMOOTH_FWHM
#    bash submit_preproc.sh --stages smooth
#
#    # one subject, in this shell, to see what happens
#    bash submit_preproc.sh --sub sub-01s13 --ses ses-01 --stages smooth --local
#
#    # only the first three lines of the subject list
#    bash submit_preproc.sh --array 1-3
#
#  Options:
#    --stages LIST   Stages to run (default: PREPROC_STAGES from the config).
#                    One or more of: topup, vdm, realign, coreg, smooth
#                    Shorthands: all, fieldmap, post_fieldmap
#    --sub SUB       Process a single subject...
#    --ses SES       ...and this session (both required together).
#    --local         Run here in the foreground instead of submitting to SLURM.
#                    Requires --sub/--ses.
#    --array RANGE   SLURM array range (default: 1-<lines in the subject list>).
#    --list FILE     Subject-session list (default: <script dir>/subses_list.txt).
#    --config FILE   Config file (default: <script dir>/pipeline_config.cfg).
#    --dependency D  SLURM dependency spec, e.g. afterok:12345.
#    --dry-run       Print the command that would run, then stop.
#    -h, --help      Show this help.
# =============================================================================

set -euo pipefail

# Where this script and its siblings live. Kept separate from the config's
# SCRIPT_DIR, which sourcing the config below overwrites.
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG_FILE="${SELF_DIR}/pipeline_config.cfg"
STAGES=""
SUB=""
SES=""
LIST=""
ARRAY=""
DEPENDENCY=""
RUN_LOCAL=false
DRY_RUN_SUBMIT=false

show_help() {
    awk 'NR>1 { if ($0 !~ /^#/) exit; print }' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --stages)     STAGES="$2";     shift 2 ;;
        --sub)        SUB="$2";        shift 2 ;;
        --ses)        SES="$2";        shift 2 ;;
        --list)       LIST="$2";       shift 2 ;;
        --config)     CONFIG_FILE="$2"; shift 2 ;;
        --array)      ARRAY="$2";      shift 2 ;;
        --dependency) DEPENDENCY="$2"; shift 2 ;;
        --local)      RUN_LOCAL=true;  shift ;;
        --dry-run)    DRY_RUN_SUBMIT=true; shift ;;
        -h|--help)    show_help; exit 0 ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            echo "       Run with --help to see the accepted options." >&2
            exit 1 ;;
    esac
done

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE" >&2
    exit 1
fi
source "$CONFIG_FILE"

: "${BIDS_ROOT:?ERROR: BIDS_ROOT not set in pipeline_config.cfg}"

SBATCH_SCRIPT="${SELF_DIR}/spm_preproc_array.sbatch"
if [[ ! -f "$SBATCH_SCRIPT" ]]; then
    echo "ERROR: spm_preproc_array.sbatch not found in ${SELF_DIR}" >&2
    exit 1
fi

# Sourcing the config set SCRIPT_DIR to the pipeline's script directory on the
# cluster — that is where subses_list.txt lives. Fall back to this script's own
# directory when the config does not set it.
LIST_DIR="${SCRIPT_DIR:-$SELF_DIR}"

[[ -z "$STAGES" ]] && STAGES="${PREPROC_STAGES:-all}"
STAGES="${STAGES// /}"

# ---- Validate the stage names before anything is submitted ----
VALID_STAGES="topup vdm realign coreg smooth all fieldmap fieldmaps fmap post_fieldmap postfieldmap post-fieldmap unwarp realign_unwarp"
IFS=',' read -ra _stage_parts <<< "$STAGES"
for _s in "${_stage_parts[@]}"; do
    [[ -z "$_s" ]] && continue
    if [[ ! " ${VALID_STAGES} " == *" ${_s} "* ]]; then
        echo "ERROR: Unknown stage: '${_s}'" >&2
        echo "       Valid stages: topup, vdm, realign, coreg, smooth" >&2
        echo "       Shorthands:   all, fieldmap, post_fieldmap" >&2
        exit 1
    fi
done

if [[ -n "$SUB" || -n "$SES" ]]; then
    if [[ -z "$SUB" || -z "$SES" ]]; then
        echo "ERROR: --sub and --ses must be given together." >&2
        exit 1
    fi
fi

if [[ "$RUN_LOCAL" == "true" && -z "$SUB" ]]; then
    echo "ERROR: --local needs --sub and --ses (it processes one session at a time)." >&2
    exit 1
fi

# =============================================================================
#  Local run — one session, in this shell
# =============================================================================
if [[ "$RUN_LOCAL" == "true" ]]; then
    cmd=(bash "$SBATCH_SCRIPT" --config "$CONFIG_FILE" --stages "$STAGES" --sub "$SUB" --ses "$SES")
    echo "Running locally: ${cmd[*]}"
    if [[ "$DRY_RUN_SUBMIT" == "true" ]]; then
        echo "(dry run — not executing)"
        exit 0
    fi
    PIPELINE_SCRIPT_DIR="$SELF_DIR" "${cmd[@]}"
    exit $?
fi

# =============================================================================
#  SLURM submission
# =============================================================================
if [[ "$DRY_RUN_SUBMIT" != "true" ]] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch not found. Use --local --sub SUB --ses SES to run here instead." >&2
    exit 1
fi

# Pass the stage list through to the job script
job_args=(--config "$CONFIG_FILE" --stages "$STAGES")

if [[ -n "$SUB" ]]; then
    # Single session: a one-task array so the SLURM output naming stays the same
    ARRAY="${ARRAY:-1-1}"
    job_args+=(--sub "$SUB" --ses "$SES")
    n_jobs=1
else
    [[ -z "$LIST" ]] && LIST="${LIST_DIR}/subses_list.txt"
    if [[ ! -f "$LIST" ]]; then
        echo "ERROR: Subject list not found: $LIST" >&2
        echo "       Run step0_make_subses_list.sh first, or pass --sub/--ses." >&2
        exit 1
    fi
    n_jobs=$(grep -cve '^[[:space:]]*$' "$LIST")
    if [[ "$n_jobs" -eq 0 ]]; then
        echo "ERROR: Subject list is empty: $LIST" >&2
        exit 1
    fi
    ARRAY="${ARRAY:-1-${n_jobs}}"
    job_args+=(--list "$LIST")
fi

# Job name carries the stages so partial runs are recognisable in squeue
job_name="spm_preproc"
[[ "$STAGES" != "all" ]] && job_name="spm_${STAGES//,/_}"

sbatch_args=(
    --export=ALL,PIPELINE_SCRIPT_DIR="${SELF_DIR}"
    --job-name="${job_name}"
    --partition="${SLURM_PARTITION}"
    --time="${PREPROC_TIME}"
    --cpus-per-task="${PREPROC_CPUS}"
    --mem-per-cpu="${PREPROC_MEM_PER_CPU}"
    --array="${ARRAY}"
)
[[ -n "$DEPENDENCY" ]] && sbatch_args+=(--dependency="${DEPENDENCY}")
if [[ -n "${NOTIFY_EMAIL:-}" ]]; then
    sbatch_args+=(--mail-user="${NOTIFY_EMAIL}" --mail-type=END,FAIL)
fi

echo "============================================"
echo "  Submitting SPM preprocessing"
echo "  Stages:   ${STAGES}"
echo "  Mode:     ${PREPROC_MODE:-<not set>}"
echo "  Array:    ${ARRAY}   (${n_jobs} session(s) in the list)"
echo "  Job name: ${job_name}"
echo "  Config:   ${CONFIG_FILE}"
echo "============================================"

if [[ "$DRY_RUN_SUBMIT" == "true" ]]; then
    echo "sbatch ${sbatch_args[*]} ${SBATCH_SCRIPT} ${job_args[*]}"
    echo "(dry run — not submitting)"
    exit 0
fi

sbatch "${sbatch_args[@]}" "$SBATCH_SCRIPT" "${job_args[@]}"
