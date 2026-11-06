# roi_tsnr — keep only the high-tSNR voxels of each ROI

Writes filtered copies of the left/right masks, rebuilds `hemi-bilat` as the
union of the two **filtered** sides, and writes the tSNR image it used.
Filenames are unchanged, so `config.sh` needs no edit — point the ROI
directory at the new folder. The folder name carries the rule, and
`filter_info.txt` carries the rest.

```bash
# 1. look before you cut
bash tools/roi_tsnr/run_roi_tsnr_filter.sh --dry-run --local

# 2. write the masks; left and right can differ
bash tools/roi_tsnr/run_roi_tsnr_filter.sh --keep 60,40
```

One SLURM array task **per subject × session**, then a merge job. `--local`
runs them serially, which is what you want inside an interactive session:

```bash
srun --pty -c 2 --mem 8G bash
```

`--help` lists every option. Deleting this folder breaks nothing.

---

## Measured on your session-1 data

Right ErC tSNR is lower than left in **6 of 6 subjects** — median 12.1 vs 16.1,
so the right side runs at 75% of the left. Your asymmetry is real:

| subject | left | right | R/L | brain |
|---|---|---|---|---|
| sub-01s17 | 14.2 | 11.3 | 0.80 | 10.9 |
| sub-02s17 | 11.7 | 10.1 | 0.87 | 10.5 |
| sub-04s17 | 17.3 | 12.1 | 0.70 | 11.1 |
| sub-05s17 | 17.2 | 12.1 | 0.70 | 11.3 |
| sub-06s17 | 15.5 | 14.7 | 0.94 | 10.9 |
| sub-08s17 | 16.7 | 15.1 | 0.90 | 9.9 |

`princeErC` is 28–46 voxels on the left and 40–53 on the right, so there is
something to select from. (pmErC is 8–10 per side — far too few. If you filter
that ROI, fix `mthresh` and use the whole ErC first.)

### Smoothing erases the asymmetry you are trying to see

Same subject, same masks, only the input series differs:

| series | brain | left | right |
|---|---|---|---|
| `usub-*` (unsmoothed) | 10.9 | **14.2** | **11.3** |
| `susub-*` (5 mm) | 54.1 | 59.5 | 59.3 |

Smoothing multiplies tSNR by **~5.0** (predicted 5.4 for an ideal kernel on
independent noise) — and it flattens left and right to within 0.3%. The 26%
asymmetry is simply gone. A 5 mm kernel averages ~30 voxels' worth of variance,
which is wider than a 40-voxel ROI, so smoothed tSNR measures how much tissue
surrounds a voxel, not the voxel.

**So the ~74–82 ErC tSNR figures in your notes are a smoothed number.**
Unsmoothed the ErC sits at 11–17 against a whole-brain 10.9 — still at or above
whole-brain, so still not dropout, but nowhere near 75. The launcher refuses an
`s*` prefix for this reason.

---

## Choosing the rule

**`--keep L,R` (rank) is the default and usually right.** n stays constant
across subjects and sessions, which is what a group test wants.

**`--min-tsnr L,R` (absolute)** equalises *quality* instead and lets n vary,
which weights subjects unequally. Use it only with a reason.

Those are the two horns of the asymmetry. `--keep 50,50` on sub-01 keeps 20
left voxels at tSNR ≥ 14.3 and 27 right at ≥ 11.3 — equal counts, unequal
quality. An absolute cut flips that. Either is defensible; pick one, write it
down, and apply it to everybody.

`--combine min` instead of `mean` is worth considering: SPM's implicit mask is
an **intersection** across runs, so the worst run already decides whether a
voxel survives at all.

### Pick the threshold on reproducibility, never on the contrast

Run `--dry-run`, read the sweep table, and choose. Then judge the choice by
whether φ becomes more repeatable — run-to-run spread, or median `R-runs` from
`tools/grid_maps` — **not** by whether the contrast grew. The contrast is
`cos(6Δφ)` with a sign flip at 15°, so sweeping thresholds and keeping the best
one finds significance in noise reliably. If reliability is flat or worsens as
n shrinks, there is nothing to select, and you have learned that without
spending a forking path.

### The pre-post trap

Default is `--sessions ses-01`. Re-selecting voxels per session makes the ROI a
different set each time, so a tSNR shift from shim, coil loading or motion
reads as an intervention effect. Fix the set on session 1 and carry it forward;
if you do re-select, compare the sets and know the overlap before it becomes a
result.

---

## What the numbers mean

`roi_tsnr_sweep.csv` — one row per subject × session × hemisphere × threshold:
`keep_pct`, `n_kept`, `cut_tsnr`. This is the decision table.

`roi_tsnr.csv` — what was actually applied: `n_roi`, `n_kept`, `pct_kept`,
`cut_tsnr`, and `med_kept` / `med_dropped` / `med_all` so you can see how far
apart the two groups are. **If `med_kept` and `med_dropped` are close, the ROI
has no bad tail and filtering is just shrinking n.**

`<sub>_<ses>_tsnr.nii` — the tSNR image. Open it with the ROI on top: a solid
dark block on one side is dropout, scattered dark voxels are noise.

### How tSNR is defined here

- **numerator**: mean over time of the raw series.
- **denominator**: SD of the residual after removing a DCT basis at the GLM's
  own cutoff — the noise the GLM actually sees, not drift it removes anyway.
  `tools/tsnr_report.sh` makes the same choice with an FSL Gaussian filter; a
  DCT set is what SPM itself uses, so this matches the GLM more closely.

tSNR is **direction-blind**, which is the point: selecting on it cannot bias
the GLM2 contrast. Two things it is not blind to, both small and pointing
opposite ways — the SD includes task-driven variance, so a strongly responding
voxel reads slightly *lower*; whereas `ResMS` has the grid regressors removed,
so a responding voxel would read slightly *higher*. Never select on
`voxelwiseAmplitude`: that one is direction-dependent and is what GLM2 tests.

---

## A note for whoever edits the sbatch

It inherits the launcher's exports and **then** sources `config.sh`, so any
exported name that `config.sh` also defines is silently overwritten by the
config value. That is why the task list is exported as `RTF_TASKS`:
`config.sh` defines `TASKS=1,2,3`, which arrives as `1,2,3` and is not MATLAB
syntax. Check for a collision before adding a variable.
