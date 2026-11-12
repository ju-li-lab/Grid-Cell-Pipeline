# Grid Cell Pipeline

SPM preprocessing and GridCAT grid-cell analysis for BIDS fMRI data, built to
run as SLURM jobs on an HPC cluster.

Four stages: discover subjects, preprocess with SPM, flatten the output into
the layout GridCAT expects, run the GLMs.

---

## Quick start

Edit the two settings files, then run the pipeline:

```bash
bash run_all.sh
```

That runs subject discovery locally and queues the three cluster stages as one
SLURM dependency chain. `bash run_all.sh --check` checks the configuration
without submitting anything.

To run one stage on its own, call its script directly. **`bash` submits** a
job and prints its id; `--local` runs it in the foreground instead:

```bash
bash 00_list_subjects.sh                                    # rebuild the subject list
bash submit_preproc.sh --stages smooth                      # re-smooth only
bash submit_preproc.sh --sub sub-01s13 --ses ses-01 --local # one session, here, now
bash 02_prepare_gridcat.sh
bash 03_gridcat.sh
```

To make a submitted stage wait for another job:

```bash
STAGE_DEPENDENCY=afterok:12345 bash 02_prepare_gridcat.sh
```

---

## Settings

Two files, each read by the language that uses it:

| File | Holds | Read by |
|---|---|---|
| `config.sh` | paths, SLURM resources, module names, subject filters, modes, filename patterns | the shell scripts; MATLAB reads shared keys via `matlab/read_config.m` |
| `spm_options.m` | acquisition parameters and everything that goes into an SPM batch | MATLAB, as real values |

`spm_options.m` returns option structs in SPM's own shape, which the batch
builders hand to SPM wholesale — so **any option SPM accepts can be set by
adding a field there**, no parser or mapping table to extend.

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

The ROI filenames must contain `ROI_LABEL` from `config.sh`. Assumptions: ROI
masks are label masks (integer-valued); all tasks entering one GLM share
comparable acquisition parameters; each run has exactly one `rp*.txt`.

---

## Where output goes

**Nothing is ever written into `BIDS_ROOT`.** SPM writes the transforms it
estimates into the headers of its *input* images, so every input is copied
into a derivatives tree first and the copy is what gets touched.

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
    ROI/                            masks as resliced by the preprocessing
    ROI_tsnr/                       tSNR-filtered masks (tools/roi_tsnr)
  GLM_output[_<RUN_VARIANT>]/       GLM results
    run_info.tsv                    the settings that produced this run
    *_grid_metrics.txt              GridCAT's grid metrics
    *_voxel_ori.csv                 per-voxel grid orientation per run
```

`spm-preproc` is regenerated from scratch on every full run (`DERIV_RESET`).
`RUN_VARIANT` names the GLM output folder, so several analyses can coexist;
each submission snapshots the settings, so you can change `RUN_VARIANT` and
submit the next variant while one is still running. The variant name can stay
short — `run_info.tsv` records the settings, so nothing needs encoding in the
folder name.

---

## Preprocessing stages

| Stage | Does | Writes |
|---|---|---|
| `topup` | FSL topup distortion field | `func/topup_results_*` |
| `vdm` | one voxel displacement map per task | `func/vdm_task-<N>.nii` |
| `slicetime` | shift every slice to the same point in the TR | `a*_bold.nii` or `au*_bold.nii` |
| `realign` | apply the correction + motion correction | `u*_bold.nii`, `rp_*.txt`, `meanu_session.nii` |
| `coreg` | T2w → session mean, reslice the ROIs | `anat/r*_mask.nii` |
| `smooth` | Gaussian smoothing | `func/su*_bold.nii` |

Each stage reads what it needs from disk, so you can stop after any of them,
inspect or replace the intermediate files, and resume:

```bash
bash submit_preproc.sh --stages vdm,realign,coreg,smooth   # resume from topup output
bash submit_preproc.sh --stages smooth                     # just re-smooth
```

Shorthands: `all`, `fieldmap` (= topup,vdm), `post_fieldmap` (= everything
after the fieldmap stages), `stc` (= slicetime). Stages always run in pipeline
order whatever order you list them in.

### Slice timing

Off unless `DO_SLICE_TIMING=1`. Reads **`SliceTiming` from each run's own JSON
sidecar** and fails loudly if it is absent — correcting with the wrong order is
worse than not correcting. Everything is shifted to
`C.slicetime.ref_fraction` of the TR (default 0.5), which must agree with the
GLM's `microtimeOnset = ref_fraction * microtimeResolution`; the stage checks
the pair.

`SLICETIME_ORDER` decides whether the correction runs before realignment
(realignment sees corrected data, prefix `ua`) or after it (the motion
estimate is made on the untouched EPI, prefix `au`). Both are correct SPM
usage; for an ROI at the slab edge the `after` argument has teeth, so
comparing the two is worth the disk space.

The prefix the later stages read is **derived, never stored** — from
`DO_SLICE_TIMING`, `SLICETIME_ORDER` and `SMOOTH_FWHM`:

| | `before` | `after` |
|---|---|---|
| `DO_SLICE_TIMING=0` | `u` | `u` |
| + `SMOOTH_FWHM>0` | `su` | `su` |
| `DO_SLICE_TIMING=1` | `ua` | `au` |
| both | `sua` | `sau` |

(`func_prefix()` in `bids_common.sh`; the same rule in `run_spm_preproc.m` and
`run_gridcat_analysis.m`.) After changing `DO_SLICE_TIMING`, re-run from
`slicetime` onwards; after changing `SLICETIME_ORDER`, re-run
`--stages post_fieldmap` and keep the two orders in separate trees.

---

## Distortion correction modes

`PREPROC_MODE` in `config.sh`:

| Mode | Needs | How |
|---|---|---|
| `realign_unwarp` | `phasediff` + `magnitude1` in `fmap/` | SPM's FieldMap toolbox builds the VDM |
| `topup` | a reverse phase-encode EPI | FSL topup; `TOPUP_APPLY_METHOD` either converts the field to a VDM, or runs `applytopup` + Realign only |
| `precalc_fieldmap` | a ready-made B0 fieldmap | for sessions with a phasediff but **no** magnitude: `tools/make_fieldmaps.sh` builds it from the T1w |
| `realign_only` | nothing | motion correction only |

Check the acquisition parameters against your JSON sidecars — wrong values are
silently wrong results, not errors:

```bash
grep -E 'TotalReadoutTime|EstimatedTotalReadoutTime' sub-XX/ses-YY/func/*_bold.json
grep -E 'EchoTime[12]' sub-XX/ses-YY/fmap/*_phasediff.json
grep PhaseEncodingDirection sub-XX/ses-YY/func/*_bold.json
```

BIDS stores seconds; `TOTAL_READOUT_MS` and the echo times in `spm_options.m`
are **milliseconds**, `C.topup.readout_sec` is **seconds**.

---

## Multiple runs of the same scan

Run numbers are counted per modality, so `run-2` of the T2w need not belong
with `run-2` of a task. The pipeline decides in this order:

1. your explicit choice in `run_selection.tsv` (always in the pipeline
   directory, next to `config.sh`) — write it with `bash tools/scan_multirun.sh`
2. the run acquired closest in time, from `AcquisitionTime` in the sidecars
3. nothing else is safe — with `STRICT_RUN_MATCHING=true` it stops and asks

Whatever was chosen is recorded in `preproc_run_selection.tsv` next to the
output, which the GridCAT prep reads so the event files line up.

---

## ROI masks

```
ASHS segmentation
   └─ tools/ashs_extract_labels.sh   (configured by tools/ashs_extract_config.sh)
        └─ per-hemisphere masks
             └─ tools/move_rois.sh   imports them into <DERIV_ROOT>/rois/
                  └─ the coreg stage reslices them into EPI space
```

### When the coregistration will not hold

A T2w slab and a partial-volume task EPI overlap in only a small block of
tissue, and NMI over a small overlap has a degenerate optimum: SPM can drive
the images apart and write a mask that misses the brain without complaining.
Two settings deal with it:

```bash
COREG_VIA_TASK=rest        # route T2w -> whole-brain EPI -> session mean
COREG_MIN_ROI_VOXELS=5     # fail the job instead of writing an empty mask
```

The masks are still resampled **once**: both coregistration steps are
estimate-only and the single reslice applies the composed transform. After
changing either setting, `--stages coreg` is enough.

### tSNR filtering

`tools/roi_tsnr` drops low-tSNR outlier voxels from each mask (voxels more
than 1.5 IQR below Q1 of the subject's ROI tSNR distribution, all sessions)
and writes the result to `GLM_runauto/ROI_tsnr`. Set `ROI_TSNR_FILTER=1` in
`config.sh` to make the GLM use those masks; the flag is recorded in
`run_info.tsv`. See `tools/roi_tsnr/README.md`.

---

## GridCAT output and R analysis

Every GLM output folder carries `run_info.tsv` — one row with the settings
that produced it (smoothing, mask patterns, event usage, regressor method,
masking threshold, symmetry, averaging, tSNR filter, ...). The R importer
joins it onto every result row, so GridCAT runs are compared by their
settings, not by their folder names.

```
R/
  import_gridcat.R        import one or MORE GLM output dirs ->
                          run_info, grid_runs, grid_stability, voxel_ori
  plot_gridcat_run.R      one run: contrast / Rayleigh z / NaN share per run
  compare_gridcat_runs.R  several runs: which settings differ, paired metric
                          comparison, phi agreement between runs
  plot_voxel_stability.R  voxel grid orientation across runs as a spaghetti
                          plot per subject x session
  subject_info.csv        subject-level info (group, age, ...), joined on
                          subject_id; extend it with any columns you need
  look_at_tsnr.R          tSNR tables from tools/tsnr_report.sh
  motion.R                within- and between-task motion QC
  paths.R                 what route each run travelled (event files)
```

Set `GRIDCAT_DIRS` in `import_gridcat.R` to the output folder(s), source it,
then source whichever plot script you need. Session-level scores (task
performance etc.) can later be joined onto any table by `subject_id` +
`session`.

---

## Layout

```
config.sh               settings: paths, SLURM, filters          <- you edit this
spm_options.m           settings: acquisition and SPM parameters <- you edit this
run_all.sh              preflight checks + the whole SLURM chain

00_list_subjects.sh     discover subject-sessions, apply filters
01_preprocess.sbatch    SPM preprocessing, one array task per session
02_prepare_gridcat.sh   flatten the output into GridCAT's layout
03_gridcat.sh           GLM1 + GLM2, grid metrics, voxel orientations, run_info

submit_preproc.sh       preprocessing submitter: stages, single sessions, --local
bids_common.sh          shared shell helpers: derivative paths, BIDS names

matlab/
  run_spm_preproc.m     the preprocessing stages, in order
  check_bids_params.m   compares the settings against the JSON sidecars
  inputs/               picking one run of each input, and recording which
  steps/                the SPM batches: realign, unwarp, coreg, smooth
  fieldmap/             building voxel displacement maps
  fsl/                  topup and applytopup wrappers
  bids_*.m              BIDS filenames, run selection, sidecars
  prepare_gridcat_directory.m
  run_gridcat_analysis.m

tools/                  run by hand, not part of the chain
  roi_tsnr/             tSNR outlier filter for the ROI masks
  grid_maps/            voxel-wise maps from GridCAT output
  scan_multirun.sh      write run_selection.tsv
  check_roi_coverage.sh how much of each ROI survives SPM's implicit mask
  ...
R/                      importing and plotting the results
```

One function per file, named after the function.

---

## Troubleshooting

**Check what a run actually did.** Every output session carries
`preproc_provenance.json` plus verbatim copies of the `config.sh` and
`spm_options.m` that produced it; every GLM output folder carries
`run_info.tsv`.

**Unwarping looks wrong.** The job log prints the peak VDM shift. A plausible
EPI distortion is a few voxels; tens of voxels means the readout time or the
field units are wrong. Corrected images looking *more* distorted is a
sign-convention error.

**A setting had no effect.** Check it belongs to the mode you are running: in
`topup` mode everything under `C.fieldmap` is inert, in `realign_unwarp` mode
it is `C.topup` that does nothing.

**Jobs.** `squeue -u $USER` to watch, `scancel <id>` to stop. SLURM logs land
in `logs/`.
