# roi_tsnr — drop low-tSNR outlier voxels from each ROI

**Rule:** within each subject × session, pool the tSNR values of all ROI
voxels (left + right) and drop every voxel more than 1.5 interquartile ranges
below the first quartile (`tSNR < Q1 − 1.5·IQR`). A standard outlier fence,
computed within each subject — no threshold is carried between subjects.
All sessions are filtered by default.

Writes filtered copies of the left/right masks (filenames unchanged), rebuilds
`hemi-bilat` as the union of the two **filtered** sides, and writes the tSNR
image it used. Output defaults to `GLM_runauto/ROI_tsnr`, which `03_gridcat.sh`
reads when `ROI_TSNR_FILTER=1` in `config.sh`.

```bash
# 1. look before you cut
bash tools/roi_tsnr/run_roi_tsnr_filter.sh --dry-run --local

# 2. write the masks
bash tools/roi_tsnr/run_roi_tsnr_filter.sh

# 3. use them
#    config.sh: ROI_TSNR_FILTER=1, then re-run 03_gridcat.sh
```

One SLURM array task per subject × session, then a merge job. `--local` runs
them serially (for an interactive session: `srun --pty -c 2 --mem 8G bash`).
`--help` lists every option. Deleting this folder breaks nothing.

## How tSNR is defined here

- **numerator**: mean over time of the raw series.
- **denominator**: SD of the residual after removing a DCT basis at the GLM's
  own cutoff — the noise the GLM actually sees, not drift it removes anyway.

tSNR is **direction-blind**, so selecting on it cannot bias the GLM2 contrast.
Never select on `voxelwiseAmplitude`: that is what GLM2 tests.

**Always filter on unsmoothed data** (prefix `u`/`ua`). A 5 mm kernel averages
~30 voxels' worth of variance — wider than the ROI — so smoothed tSNR measures
how much tissue surrounds a voxel, not the voxel. The launcher refuses an `s*`
prefix for this reason.

## Reading the output

`roi_tsnr.csv` — one row per subject × session × hemisphere: `n_roi`,
`n_kept`, `pct_kept`, `cut_tsnr`, and `med_kept` / `med_dropped` / `med_all`.
If `med_kept` and `med_dropped` are close, the ROI had no bad tail and the
filter did (correctly) little.

`<sub>_<ses>_tsnr.nii` — the tSNR image. Open it with the ROI on top: a solid
dark block on one side is dropout, scattered dark voxels are noise.

Judge the filter by whether φ becomes more repeatable (run-to-run spread, or
the voxel-stability spaghetti plot in `R/plot_voxel_stability.R`), **not** by
whether the contrast grew — the contrast is `cos(6Δφ)` with a sign flip at
15°, so optimising it finds significance in noise reliably.
