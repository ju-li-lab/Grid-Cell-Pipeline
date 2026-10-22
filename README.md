# Grid Cell Pipeline

SPM preprocessing and GridCAT grid-cell analysis for BIDS fMRI data, built to
run as SLURM jobs on an HPC cluster.

Four stages: discover subjects, preprocess with SPM, flatten the output into the
layout GridCAT expects, run the GLMs.

---

## Quick start

Edit the two settings files, then run the pipeline:

```bash
bash run_all.sh
```

That runs subject discovery locally and queues the three cluster stages as one
SLURM dependency chain, so there is nothing to come back and launch.

To check the configuration without submitting anything:

```bash
bash run_all.sh --check
```

To run one stage on its own, call its script directly — each takes the same
options independently:

```bash
bash 00_list_subjects.sh                                   # rebuild the subject list
bash submit_preproc.sh --stages smooth                     # re-smooth only
bash submit_preproc.sh --sub sub-01s13 --ses ses-01 --local # one session, here, now
bash 02_prepare_gridcat.sh
bash 03_gridcat.sh
```

`00_list_subjects.sh` runs here and now — it only reads the BIDS directory. The
three cluster stages behave alike: **`bash` submits** a job and prints its id,
and `--local` runs it in the foreground instead, which is what you want when
debugging:

```bash
bash submit_preproc.sh --sub sub-01s13 --ses ses-01 --local
bash 02_prepare_gridcat.sh --local
bash 03_gridcat.sh --local
```

`sbatch 0N_...` also works for stages 2 and 3; SLURM runs the submitter, which
submits the real job.

To make a submitted stage wait for another job:

```bash
STAGE_DEPENDENCY=afterok:12345 bash 02_prepare_gridcat.sh
```

---

## Settings

Two files, each read by the language that uses it. Nothing else needs editing.

| File | Holds | Read by |
|---|---|---|
| `config.sh` | paths, SLURM resources, module names, subject filters, modes, filename patterns | sourced directly by the shell scripts; MATLAB reads the few keys it shares via `matlab/read_config.m` |
| `spm_options.m` | acquisition parameters and everything that goes into an SPM batch | MATLAB, as real values |

`spm_options.m` is a MATLAB function returning option structs in SPM's own
shape, which the batch builders hand to SPM wholesale:

```matlab
matlabbatch{1}.spm.spatial.realignunwarp.eoptions = CFG.spm.realignunwarp.eoptions;
```

So **any option SPM accepts can be set by adding a field to `spm_options.m`** —
no parser to extend, no mapping table, no validator to update. Field names are
SPM's own; look them up in the SPM batch editor or in `spm_cfg_realignunwarp.m`
and `spm_cfg_fieldmap.m`.

To see the current shell settings: `cat config.sh`.

Every run copies both files next to its output, so a result always carries the
settings that produced it.

---

## Input layout

```
BIDS_ROOT/
  sub-XX/
    ses-YY/
      anat/
        sub-XX_ses-YY_T2w.nii.gz                             anatomical reference
        sub-XX_ses-YY_hemi-left_label-ErC_mask.nii.gz        left ROI  (label-valued)
        sub-XX_ses-YY_hemi-right_label-ErC_mask.nii.gz       right ROI (label-valued)
      fmap/
        sub-XX_ses-YY_magnitude1.nii.gz        GRE fieldmap, for realign_unwarp
        sub-XX_ses-YY_phasediff.nii.gz
        sub-XX_ses-YY_phasediff.json
      func/
        sub-XX_ses-YY_task-run1_bold.nii.gz
        sub-XX_ses-YY_task-run1_bold.json
        sub-XX_ses-YY_task-reverse_bold.nii.gz reverse-PE EPI, for topup
```

The ROI filenames must contain `ROI_LABEL` from `config.sh` (`ErC` above). That
one setting is what the preprocessing, the flattening and the GLMs all filter on.

### Assumptions

- ROI masks are **label masks** (integer-valued), not probability maps.
- All tasks entering one GLM have comparable acquisition parameters —
  phase-encode direction, readout time, voxel size, TR. Runs that differ need
  separate preprocessing or exclusion.
- Each run has exactly one motion regressor file (`rp*.txt`).
- You know which space every file is in at each stage: native T2w, native EPI,
  unwarped EPI, resliced EPI.

---

## Where output goes

**Nothing is ever written into `BIDS_ROOT`.** SPM's Realign & Unwarp and
Coregister write the transforms they estimate into the headers of their *input*
images, and reading a `.nii.gz` leaves a decompressed twin behind. So every input
is copied into a derivatives tree first and the copy is what gets touched.

```
<DERIV_ROOT>/                       (default: OUTPUT_ROOT/derivatives)
  spm-preproc/sub-XX/ses-YY/
    func/   u*_bold.nii, rp_*.txt, vdm_task-*.nii, topup_*, meanu_session.nii, su*
    anat/   the T2w and the resliced ROI masks (r*_mask.nii)
    fmap/   the fieldmap inputs
    config.sh, spm_options.m, preproc_provenance.json
  rois/           ROI masks imported by tools/move_rois.sh    (kept)
  fieldmaps/      fieldmaps built by tools/make_fieldmaps.sh  (kept)

<OUTPUT_ROOT>/
  GLM_runauto/                      GridCAT input, flattened
  GLM_output[_<RUN_VARIANT>]/       GLM results, *_grid_metrics.txt
```

`spm-preproc` is regenerated from scratch on every full run — see `DERIV_RESET`
in `config.sh`. A resumed run keeps what is there, because that is its input.

`RUN_VARIANT` names the GLM output folder, so several analyses can coexist. Each
submission snapshots the settings, so you can change `RUN_VARIANT` and submit the
next variant while one is still running.

---

## Preprocessing stages

The SPM preprocessing is five stages. Each reads what it needs from disk, so you
can stop after any of them, inspect or replace the intermediate files, and resume.

| Stage | Does | Writes |
|---|---|---|
| `topup` | FSL topup distortion field | `func/topup_results_*` |
| `vdm` | one voxel displacement map per task | `func/vdm_task-<N>.nii` |
| `realign` | apply the correction + motion correction | `u*_bold.nii`, `rp_*.txt`, `meanu_session.nii` |
| `coreg` | T2w → session mean, reslice the ROIs | `anat/r*_mask.nii` |
| `smooth` | Gaussian smoothing | `func/su*_bold.nii` |

Set the default in `PREPROC_STAGES`, or override per submission:

```bash
bash submit_preproc.sh --stages vdm,realign,coreg,smooth   # resume from topup output
bash submit_preproc.sh --stages realign,coreg,smooth       # resume from ready-made VDMs
bash submit_preproc.sh --stages smooth                     # just re-smooth
```

Shorthands: `all`, `fieldmap` (= topup,vdm), `post_fieldmap` (= realign,coreg,smooth).

---

## Distortion correction modes

`PREPROC_MODE` in `config.sh`:

| Mode | Needs | How |
|---|---|---|
| `realign_unwarp` | `phasediff` + `magnitude1` in `fmap/` | SPM's FieldMap toolbox unwraps the phase and builds the VDM |
| `topup` | a reverse phase-encode EPI | FSL topup estimates the field; `TOPUP_APPLY_METHOD` then either converts it to a VDM for Realign & Unwarp, or corrects the EPIs with `applytopup` followed by Realign only |
| `precalc_fieldmap` | a ready-made B0 fieldmap | for sessions with a phasediff but **no** magnitude: `tools/make_fieldmaps.sh` builds the missing magnitude from the T1w with FSL, then SPM converts the fieldmap to a VDM |
| `realign_only` | nothing | motion correction only. Useful as a control when you are working out whether unwarping is helping |

Check the parameters against your JSON sidecars — wrong values here are silently
wrong results, not errors:

```bash
# readout time and echo times
grep -E 'TotalReadoutTime|EstimatedTotalReadoutTime' sub-XX/ses-YY/func/*_bold.json
grep -E 'EchoTime[12]' sub-XX/ses-YY/fmap/*_phasediff.json
# phase-encode direction
grep PhaseEncodingDirection sub-XX/ses-YY/func/*_bold.json
```

BIDS stores times in seconds. `TOTAL_READOUT_MS` and the echo times in
`spm_options.m` are in **milliseconds**; `C.topup.readout_sec` is in **seconds**.

`EPI_BASED_FIELDMAP` must be `0` for a GRE phasediff fieldmap and `1` only when
the fieldmap really is EPI-based.

---

## Multiple runs of the same scan

Run numbers are counted per modality and say nothing across modalities. If the
subject climbed out of the scanner between sequences, `run-2` of the T2w and
`run-2` of a task are from different visits, and pairing them coregisters the
wrong anatomy.

So the pipeline decides in this order:

1. your explicit choice in `RUN_SELECTION_FILE`
2. the run acquired closest in time, read from `AcquisitionTime` in the sidecars,
   within `RUN_MATCH_GAP_MIN` minutes
3. nothing else is safe — with `STRICT_RUN_MATCHING=true` it stops and asks

To write the selection file:

```bash
bash tools/scan_multirun.sh
```

It lists every run with its acquisition time and pre-fills a suggestion you can
correct. Whatever was chosen is recorded in `preproc_run_selection.tsv` next to
the output, which the GridCAT prep reads so the event files line up.

---

## ROI masks

```
ASHS segmentation
   └─ tools/ashs_extract_labels.sh   (configured by tools/ashs_extract_config.sh)
        └─ per-hemisphere masks
             └─ tools/move_rois.sh   imports them into <DERIV_ROOT>/rois/
                  └─ the coreg stage reslices them into EPI space
```

`tools/ashs_extract_config.sh` is a separate config on purpose: it is a separate
tool with its own vocabulary, run by hand, and folding it into `config.sh` would
recreate the problem that file used to have. The one thing the two must agree on
is the label name — `LABEL_NAMES` there, `ROI_LABEL` here.

`move_rois.sh` is dry-run by default; pass `--execute` to actually copy.

---

## Layout

```
config.sh               settings: paths, SLURM, filters          <- you edit this
spm_options.m           settings: acquisition and SPM parameters <- you edit this
run_all.sh              preflight checks + the whole SLURM chain

00_list_subjects.sh     discover subject-sessions, apply filters
01_preprocess.sbatch    SPM preprocessing, one array task per session
02_prepare_gridcat.sh   flatten the output into GridCAT's layout
03_gridcat.sh           GLM1 + GLM2

submit_preproc.sh       preprocessing submitter: stages, single sessions, --local
bids_common.sh          shared shell helpers: derivative paths, BIDS names

matlab/
  run_spm_preproc.m     the five stages, in order
  check_bids_params.m   compares the settings against the JSON sidecars
  inputs/               picking one run of each input, and recording which
  steps/                the SPM batches: realign, unwarp, coreg, smooth
  fieldmap/             building voxel displacement maps
  fsl/                  topup and applytopup wrappers
  bids_*.m              BIDS filenames, run selection, sidecars
  prepare_gridcat_directory.m
  run_gridcat_analysis.m

tools/                  run by hand, not part of the chain
R/                      reading the results back for plotting
docs/                   pipeline map (generated, untracked)
```

One function per file, named after the function. To change how smoothing works,
open `matlab/steps/smooth_bold_files.m`; to change the topup→VDM conversion,
`matlab/fieldmap/convert_topup_to_vdm.m`.

---

## Checking a preprocessing trial

Before committing a preprocessing run to GridCAT, measure how much of each ROI
mask will actually carry a value in the GLM:

```bash
bash tools/check_roi_coverage.sh --trial sep4_fwhm5
```

It reproduces SPM's implicit mask — a voxel is kept only if it clears
`maskingThreshold x global` in **every** scan, with the global computed per
scan — and reports, per subject/session/hemisphere, how many mask voxels
survive across a range of thresholds.

Two columns decide what to do:

- **`no_data`** — voxels with essentially no signal. These are dropout or a
  misplaced mask. No threshold recovers them; the fix is the acquisition or
  the coregistration.
- **`n_<thr>`** — voxels that would carry a value at that threshold. Set
  `C.GLM.maskingThreshold` in `spm_options.m` to the lowest threshold still on
  the flat part of the curve.

Every run appends to `roi_coverage.csv` with the `--trial` label as a column,
so trials are directly comparable:

```bash
bash tools/check_roi_coverage.sh --trial sep15_fwhm1
bash tools/check_roi_coverage.sh --trial sep4_fwhm5
```

```r
cov <- readr::read_csv("roi_coverage.csv")
cov |> dplyr::group_by(trial) |> dplyr::summarise(median(pct_0.80), median(pct_0.50))
```

Use `--sub sub-05s17` to check one subject while iterating, and `--deriv DIR`
to point at a derivatives tree other than the one `config.sh` names.

---

## Troubleshooting

**Check what a run actually did.** Every output session carries
`preproc_provenance.json` (what ran, on which inputs, which software versions)
plus verbatim copies of the `config.sh` and `spm_options.m` that produced it.
`preproc_stages.log` accumulates one line per run, so a resumed run does not
erase the history.

**Unwarping looks wrong.** The job log prints the peak VDM shift. A plausible EPI
distortion is a few voxels; tens of voxels means the readout time or the field
units are wrong. If the corrected images look *more* distorted than the
originals, that is a sign-convention error, not a magnitude error.

**The settings you changed had no effect.** Check they belong to the mode you are
running. In `topup` mode the SPM FieldMap toolbox is never called, so everything
under `C.fieldmap` in `spm_options.m` is inert; in `realign_unwarp` mode it is
`C.topup` that does nothing.

**Jobs.** `squeue -u $USER` to watch, `scancel <id>` to stop. SLURM logs land in
`logs/`.

**A stage printed its banner and then stopped.** `module` is usually a shell
function from the login profile, not a real command, so it does not exist
inside a non-interactive `bash script.sh`. Every stage now says
`NOTE: no 'module' command available` and carries on with whatever is on PATH,
then fails with a clear message if `matlab` is genuinely missing. If you see
`module: command not found` you are on an older copy of these scripts.
