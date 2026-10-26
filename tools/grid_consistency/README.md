# Grid orientation consistency — preprocessing QC

**This directory is optional and self-contained.** Nothing in the pipeline
calls it. Delete `tools/grid_consistency/` and everything else keeps working.

## What it measures, and why it is separate

One question: *does a voxel report the same grid orientation in run 1 as in
run 2?*

That is a **preprocessing** question. If registration or unwarping is off, a
voxel samples different tissue in each run, so its orientation stops
reproducing. Whether aligned beats misaligned is an **analysis** question,
and it belongs to GLM2. Keeping them apart means a change to one does not
muddy your reading of the other.

```bash
bash tools/grid_consistency/run_check_grid_consistency.sh --trial sep4_fwhm5
```

## It does not depend on your GLM2 settings

It reads only the per-run images GLM1 writes:

```
voxelwiseGridOri_<event>_run<N>_deg.nii
```

and never `voxelwiseGridOri_<event>_allRunsAvg_deg.nii`. That matters twice:

- **allRunsAvg contains the run it would be compared against**, so agreement
  with it is circular. In one dataset the circular pairs read 78–88% where
  the honest run-vs-run pairs read 49–56%.
- It means nothing here depends on `GLM2_averageMeanGridOriAcrossRuns`. You
  can set it to 0 or 1 and these numbers stay comparable across trials.

## Reading the output

```
subject     session  hemi   runs      n_mask n_valid pct_stable  r_diff  z_diff
sub-05s17   ses-01   left   1_vs_2        14      14      57.1%   0.412    2.38
```

| | |
|---|---|
| `pct_stable` | GridCAT's own measure, directly comparable with your `grid_metrics` files. **Chance is printed in the header** — for 6-fold with the default 15° tolerance it is exactly 50%. |
| `r_diff` | resultant length of the orientation *difference*, 0–1. Threshold-free, so more sensitive than counting hits. |
| `z_diff` | `n_valid × r_diff²`. Rayleigh's statistic: ~1 is chance, >3 is roughly p<0.05. **This is the one to compare trials on** — it wastes no information and is not biased by small n the way `r_diff` is. |
| `r_run<N>` | between-voxel coherence within a run. This is the n-independent form of the Rayleigh z in `grid_metrics`, where `z = n·r²`. |

A constant offset between runs is perfectly *consistent* but not *aligned*:
it gives `r_diff = 1` with `pct_stable = 0`. That is why both are reported.

## Two rules for interpreting it

**Compare the same subject across preprocessing trials, never between
subjects.** Consistency only registers if there is grid signal to reproduce.
A subject with no grid modulation sits at chance however good the
preprocessing.

**Low n makes it noisy.** With 15–30 voxels per ROI a stability percentage
carries a standard error around 9–12 points. Read the median across sessions,
not a single row.

## Defaults worth knowing

- **Session 1 only.** Sessions after an intervention can differ for real
  reasons, which would confound a preprocessing comparison. `--all-sessions`
  when you want the rest.
- **Hemispheres are never pooled.** The two ErC masks differ in size and in
  how much signal survives; a pooled number would hide that asymmetry.
- **The shift test is on.** It translates the mask along each axis and
  prints the z profile:

  ```
  sub-E  ses-01  left   shift x:    1.21    3.44  111.08  192.00   88.31   peak +1
  ```

  A peak off centre means the mask and the functional data disagree about
  where the ROI is — a coregistration shift, with the sign telling you which
  way. **It prints the profile rather than a verdict on purpose:** scanning
  five shifts on three axes will find a peak in pure noise (the maximum of 15
  draws from the null sits around z = 3), so only a large, clearly peaked
  profile means anything. A flat profile is the answer you want. It also
  scans each axis on its own, not all combinations, so treat it as a pointer
  rather than a search.

## Comparing trials on a fixed voxel set

A mask mean changes when the voxel set changes, so a difference between
trials can be sampling rather than signal. `n_valid` is reported so you can
see when sets differ. To force one fixed set:

```matlab
% once per trial
check_grid_consistency(glmDir, 'Trial', 'A', 'WriteValidMask', true)
```

then intersect the resulting `validvox_*.nii` across trials
(`fslmaths a -mas b out`) and re-run every trial with
`'CommonMaskDir', <dir holding the intersection>`.

Every run appends to `grid_consistency.csv` with the `--trial` label as a
column, so trials stack up for comparison in R.

## Not included

Grid *amplitude* — `sqrt(β_sin² + β_cos²)` per voxel. GridCAT computes it in
`betas2ori` but only writes the orientation to disk, so it would have to be
recomputed from the GLM1 beta images. It is also not cross-validated, so
noise alone produces positive amplitude: useful for asking *where* signal is
relatively stronger, not for establishing that grid signal exists.
