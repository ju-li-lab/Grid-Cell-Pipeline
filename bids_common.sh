#!/usr/bin/env bash
# =============================================================================
#  bids_common.sh — shared helpers for the BIDS derivatives layout
# =============================================================================
#  Sourced by the shell scripts of the pipeline. Nothing here runs on its own.
#
#      source "${SCRIPT_DIR}/bids_common.sh"
#
#  It expects pipeline_config.cfg to have been sourced first, and provides:
#
#    Derivative paths
#      deriv_root                       <OUTPUT_ROOT>/derivatives, or DERIV_ROOT
#      deriv_dataset <name>             one derivative dataset inside it
#      deriv_preproc / deriv_rois / deriv_fieldmaps
#      deriv_session <dataset> <sub> <ses> [<modality>]
#      write_dataset_description <dir> <name> <description>
#
#    Resetting the preprocessing derivatives
#      stages_reset_from_scratch <stages> <mode>   is this a from-scratch run?
#      reset_preproc_deriv [<sub> <ses>]           wipe all of it, or one session
#
#    Reading BIDS filenames and sidecars
#      bids_entity <filename> <entity>  e.g. bids_entity "$f" run  ->  2
#      bids_suffix <filename>           e.g. _T2w.nii.gz          ->  T2w
#      bids_strip_ext <filename>
#      json_field <json> <field>
#      acq_seconds <nifti>              acquisition time from the JSON sidecar
#
#    Run selection (run_selection.tsv)
#      run_selection_get <tsv> <sub> <ses> <type>
#
#  Everything writes to stdout and returns non-zero when it has nothing to say,
#  so the callers can use  x="$(f ...)" || x=""  without set -e killing them.
# =============================================================================

# Guard against double-sourcing
[[ -n "${_BIDS_COMMON_SOURCED:-}" ]] && return 0
_BIDS_COMMON_SOURCED=1

# -----------------------------------------------------------------------------
#  Derivative paths
# -----------------------------------------------------------------------------

# The root that holds every derivative dataset this pipeline writes.
#
# Returns 1 when neither DERIV_ROOT nor OUTPUT_ROOT is set. Callers that assign
# the result under "set -e" would otherwise exit with no output at all, so say
# what is missing on the way out.
deriv_root() {
    if [[ -n "${DERIV_ROOT:-}" ]]; then
        printf '%s' "${DERIV_ROOT%/}"
    elif [[ -n "${OUTPUT_ROOT:-}" ]]; then
        printf '%s' "${OUTPUT_ROOT%/}/derivatives"
    else
        echo "ERROR: neither DERIV_ROOT nor OUTPUT_ROOT is set in pipeline_config.cfg" >&2
        echo "       The pipeline has nowhere to write its derivatives." >&2
        return 1
    fi
}

# One derivative dataset (a BIDS-style tree of its own) inside the root.
deriv_dataset() {
    local name="$1" root
    root="$(deriv_root)" || return 1
    printf '%s/%s' "$root" "${name#/}"
}

deriv_preproc()    { deriv_dataset "${DERIV_PREPROC:-spm-preproc}"; }
deriv_rois()       { deriv_dataset "${DERIV_ROIS:-rois}"; }
deriv_fieldmaps()  { deriv_dataset "${DERIV_FIELDMAPS:-fieldmaps}"; }

# deriv_session <dataset-dir> <sub> <ses> [modality]
#   -> <dataset-dir>/sub-XX/ses-YY[/anat|func|fmap]
deriv_session() {
    local dataset="$1" sub="$2" ses="$3" modality="${4:-}"
    local p="${dataset%/}/${sub}"
    [[ -n "$ses" ]] && p="${p}/${ses}"
    [[ -n "$modality" ]] && p="${p}/${modality}"
    printf '%s' "$p"
}

# A minimal BIDS derivatives dataset_description.json, written once per dataset.
write_dataset_description() {
    local dir="$1" name="$2" desc="${3:-}"
    local out="${dir%/}/dataset_description.json"
    mkdir -p "$dir"
    cat > "$out" <<EOF
{
  "Name": "${name}",
  "BIDSVersion": "1.8.0",
  "DatasetType": "derivative",
  "Description": "${desc}",
  "SourceDatasets": [
    { "URI": "${BIDS_ROOT:-}" }
  ],
  "GeneratedBy": [
    {
      "Name": "Grid-Cell-Pipeline",
      "Description": "${desc}",
      "CodeURL": "https://github.com/ju-li-lab/Grid-Cell-Pipeline"
    }
  ]
}
EOF
}

# -----------------------------------------------------------------------------
#  Resetting the preprocessing derivatives
# -----------------------------------------------------------------------------

# The first stage a given PREPROC_MODE actually performs. A stage list that
# includes it starts the preprocessing from scratch, so the previous output for
# that session is stale and can go. A list that does not (say realign,coreg,
# smooth) is a resume and must keep what is already on disk.
mode_entry_stage() {
    case "${1:-realign_only}" in
        topup)                          printf 'topup'   ;;
        realign_unwarp|precalc_fieldmap) printf 'vdm'    ;;
        *)                              printf 'realign' ;;
    esac
}

# stages_reset_from_scratch <stages> <preproc-mode>
#   Exit 0 when the given stage list re-runs the preprocessing from its start.
stages_reset_from_scratch() {
    local stages="${1:-all}" mode="${2:-}" entry part
    stages="${stages// /}"
    entry="$(mode_entry_stage "$mode")"

    local IFS=','
    for part in $stages; do
        case "$part" in
            all|fieldmap|fieldmaps|fmap) return 0 ;;
            "$entry")                    return 0 ;;
        esac
    done
    return 1
}

# Should the derivatives be wiped for this stage list? Honours DERIV_RESET.
deriv_reset_wanted() {
    local stages="${1:-all}" mode="${2:-}"
    case "${DERIV_RESET:-auto}" in
        never|false|no) return 1 ;;
        always|true|yes) return 0 ;;
        *) stages_reset_from_scratch "$stages" "$mode" ;;
    esac
}

# reset_preproc_deriv [<sub> <ses>]
#   No arguments: wipe the whole preprocessing derivative dataset.
#   With sub/ses: wipe only that session, which is what an array task does.
#
#   Refuses to delete anything that is not below the derivatives root, so a
#   mistyped DERIV_ROOT cannot take out a data directory.
reset_preproc_deriv() {
    local sub="${1:-}" ses="${2:-}"
    local base target root
    base="$(deriv_preproc)" || { echo "ERROR: cannot resolve the derivatives root (set OUTPUT_ROOT or DERIV_ROOT)" >&2; return 1; }
    root="$(deriv_root)"

    if [[ -n "$sub" ]]; then
        target="$(deriv_session "$base" "$sub" "$ses")"
    else
        target="$base"
    fi

    # Safety: the target must sit inside the derivatives root and be deeper
    # than "/", and the root must not be the BIDS root itself.
    case "$target" in
        "$root"/*|"$root") : ;;
        *) echo "ERROR: refusing to delete outside the derivatives root: $target" >&2; return 1 ;;
    esac
    if [[ -n "${BIDS_ROOT:-}" && "${target%/}" == "${BIDS_ROOT%/}" ]]; then
        echo "ERROR: refusing to delete BIDS_ROOT: $target" >&2
        return 1
    fi
    if [[ "${target}" != */*/* ]]; then
        echo "ERROR: refusing to delete a top-level path: $target" >&2
        return 1
    fi

    if [[ -d "$target" ]]; then
        rm -rf "$target"
        printf 'Removed previous preprocessing derivatives: %s\n' "$target"
    fi
    mkdir -p "$target"
    write_dataset_description "$base" "${DERIV_PREPROC:-spm-preproc}" \
        "SPM preprocessing (realign/unwarp, coregistration, smoothing)"
    return 0
}

# -----------------------------------------------------------------------------
#  Reading BIDS filenames
# -----------------------------------------------------------------------------

bids_strip_ext() {
    local n="${1##*/}"
    n="${n%.gz}"; n="${n%.nii}"; n="${n%.json}"; n="${n%.txt}"; n="${n%.tsv}"
    printf '%s' "$n"
}

# bids_entity <filename> <entity>   e.g. run -> "2" for ..._run-2_bold.nii
bids_entity() {
    local stem key val
    stem="$(bids_strip_ext "$1")"
    key="$2"
    if [[ "_${stem}_" =~ _${key}-([A-Za-z0-9]+)_ ]]; then
        val="${BASH_REMATCH[1]}"
        printf '%s' "$val"
        return 0
    fi
    return 1
}

# The BIDS suffix: the last _-separated chunk of the stem (T2w, bold, epi, ...).
bids_suffix() {
    local stem
    stem="$(bids_strip_ext "$1")"
    printf '%s' "${stem##*_}"
}

# The run label as it appears in run_selection.tsv: "run-2", or "no-run".
bids_run_label() {
    local r
    if r="$(bids_entity "$1" run)"; then
        printf 'run-%s' "$r"
    else
        printf 'no-run'
    fi
}

# -----------------------------------------------------------------------------
#  JSON sidecars
# -----------------------------------------------------------------------------

# json_field <json-file> <field>   — prints the value, or nothing if absent.
json_field() {
    local f="$1" k="$2"
    [[ -f "$f" ]] || return 1
    python3 - "$f" "$k" <<'PY' 2>/dev/null || return 1
import json, sys
try:
    with open(sys.argv[1]) as fh:
        d = json.load(fh)
except Exception:
    sys.exit(1)
v = d.get(sys.argv[2])
if v is None:
    sys.exit(1)
print(v)
PY
}

# The JSON sidecar that belongs to a NIfTI (same stem, .json).
json_sidecar() {
    local f="$1" base
    base="${f%.gz}"; base="${base%.nii}"
    if [[ -f "${base}.json" ]]; then
        printf '%s' "${base}.json"
        return 0
    fi
    return 1
}

# acq_seconds <nifti>
#   Seconds since midnight for this scan, read from AcquisitionDateTime or
#   AcquisitionTime in the JSON sidecar. Used to tell apart runs that were
#   acquired back to back from runs acquired after the subject left the
#   scanner. Prints nothing when the sidecar has no usable time.
acq_seconds() {
    local nii="$1" js v
    js="$(json_sidecar "$nii")" || return 1
    for key in AcquisitionDateTime AcquisitionTime; do
        v="$(json_field "$js" "$key")" || continue
        [[ -z "$v" ]] && continue
        python3 - "$v" <<'PY' 2>/dev/null && return 0
import re, sys
s = sys.argv[1]
m = re.search(r'(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)$', s.strip())
if not m:
    m = re.search(r'T?(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)', s.strip())
if not m:
    m = re.fullmatch(r'(\d{2})(\d{2})(\d{2}(?:\.\d+)?)', s.strip())
if not m:
    sys.exit(1)
h, mi, se = m.groups()
print('%.3f' % (int(h) * 3600 + int(mi) * 60 + float(se)))
PY
    done
    return 1
}

# series_number <nifti> — SeriesNumber from the sidecar, or nothing.
series_number() {
    local nii="$1" js
    js="$(json_sidecar "$nii")" || return 1
    json_field "$js" SeriesNumber
}

# A short human-readable acquisition stamp for the scan report: "10:42:13 #7".
acq_stamp() {
    local nii="$1" secs sn out=""
    if secs="$(acq_seconds "$nii")"; then
        out="$(python3 -c "s=float('$secs'); print('%02d:%02d:%02d' % (s//3600, (s%3600)//60, s%60))" 2>/dev/null)"
    fi
    if sn="$(series_number "$nii")"; then
        out="${out:+${out} }#${sn%.*}"
    fi
    [[ -z "$out" ]] && return 1
    printf '%s' "$out"
}

# -----------------------------------------------------------------------------
#  run_selection.tsv
# -----------------------------------------------------------------------------

# run_selection_get <tsv> <sub> <ses> <type>
#   Prints the selected run for one entry, e.g. "run-2". Nothing when the file
#   does not exist, the row is missing, or the choice was left blank.
#
#   Column positions are read from the header when there is one, so the old
#   five-column file (subject session type available_runs selected_run) and
#   the current one both work.
run_selection_get() {
    local tsv="$1" sub="$2" ses="$3" type="$4"
    [[ -f "$tsv" ]] || return 1

    python3 - "$tsv" "$sub" "$ses" "$type" <<'PY' 2>/dev/null || return 1
import sys

path, sub, ses, want = sys.argv[1:5]
cols = {'subject': 0, 'session': 1, 'type': 2, 'selected_run': 4}

with open(path) as fh:
    for lineno, raw in enumerate(fh):
        row = raw.rstrip('\n').split('\t')
        if not row or not row[0].strip():
            continue
        if lineno == 0 and row[0].strip().lower() in ('subject', 'sub'):
            cols = {name.strip().lower(): i for i, name in enumerate(row)}
            continue
        def get(name):
            i = cols.get(name)
            return row[i].strip() if i is not None and i < len(row) else ''
        if get('subject') != sub or get('session') != ses:
            continue
        if get('type') != want:
            continue
        sel = get('selected_run')
        if sel:
            print(sel)
            sys.exit(0)
sys.exit(1)
PY
}

# The path of the run selection file, honouring RUN_SELECTION_FILE.
run_selection_file() {
    if [[ -n "${RUN_SELECTION_FILE:-}" ]]; then
        printf '%s' "$RUN_SELECTION_FILE"
    else
        printf '%s' "${SCRIPT_DIR:-.}/run_selection.tsv"
    fi
}

# -----------------------------------------------------------------------------
#  Picking one file out of several runs
# -----------------------------------------------------------------------------

# pick_run <selected-run-label-or-empty> <file...>
#   Prints the file to use.
#     - one candidate                -> that one
#     - a selection like "run-2"     -> the candidate carrying that run entity
#     - "no-run"                     -> the candidate without a run entity
#     - nothing selected             -> the highest run number (sorted last)
#   Returns 1 when a selection was given but matches nothing, so the caller can
#   decide whether that is fatal.
pick_run() {
    local sel="$1"; shift
    local f label
    [[ $# -eq 0 ]] && return 1

    if [[ $# -eq 1 ]]; then
        printf '%s' "$1"
        return 0
    fi

    if [[ -n "$sel" ]]; then
        for f in "$@"; do
            label="$(bids_run_label "$f")"
            if [[ "$label" == "$sel" || "$label" == "run-${sel#run-}" ]]; then
                printf '%s' "$f"
                return 0
            fi
        done
        return 1
    fi

    # No selection: the last in sorted order, which is the highest run number.
    printf '%s' "$(printf '%s\n' "$@" | sort | tail -1)"
    return 0
}

# pick_nearest <anchor-seconds> <file...>
#   Prints the file acquired closest to the anchor time, and its gap in seconds
#   on the second line. Returns 1 when no candidate has an acquisition time, so
#   the caller can fall back to run numbers.
pick_nearest() {
    local anchor="$1"; shift
    local f secs best="" bestgap=""

    [[ -z "$anchor" || $# -eq 0 ]] && return 1

    for f in "$@"; do
        secs="$(acq_seconds "$f")" || continue
        local gap
        gap="$(python3 -c "print(abs(float('$secs') - float('$anchor')))" 2>/dev/null)" || continue
        if [[ -z "$bestgap" ]] || python3 -c "import sys; sys.exit(0 if float('$gap') < float('$bestgap') else 1)"; then
            best="$f"; bestgap="$gap"
        fi
    done

    [[ -z "$best" ]] && return 1
    printf '%s\n%s\n' "$best" "$bestgap"
    return 0
}
