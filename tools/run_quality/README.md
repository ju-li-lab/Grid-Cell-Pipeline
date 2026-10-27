# Per-run data quality

**This directory is optional and self-contained.** Nothing in the pipeline
calls it. Delete `tools/run_quality/` and everything else keeps working.

```bash
bash tools/run_quality/run_check_run_quality.sh --trial sep4_fwhm5
```

## Read the reliability line before acting on anything

A per-run number is only useful if the same run scores the same twice. This
tool measures that for you, every time, with a split-half check — and the
metric it was originally built around **failed**:

```
---- can you act on task_r2? ----
  split-half correlation : r = 0.099   95% CI [-0.39 0.54]   (18 runs)
  Spearman-Brown         : r = 0.180
  -> NO EVIDENCE it works. ... Use tsnr_roi and mean_fd.
```

The check needs **at least 10 runs** and refuses to pronounce below that — on
a 6-run subset of the same data it read `r = 0.891`, the opposite conclusion.
A correlation from a handful of points is itself noise, so the verdict uses
the upper end of the confidence interval, not the point estimate.

The first and second half of the same run disagreed almost completely, so a
run's `task_r2` is mostly noise. Dropping runs on it would be dropping them at
random. It is still computed, because on data with a stronger task response it
may well be reliable — but the tool checks rather than assumes.

## What it reports

| | |
|---|---|
| `tsnr_brain` | median temporal SNR over the brain. Model-free, normally the most reliable number here. |
| `tsnr_roi` | the same inside the ROI mask. **This is the one that matters for an ROI analysis** — the brain can look fine while the ErC does not. |
| `mean_fd` | framewise displacement (Power), mm: `|Δtrans| + 50·|Δrot|`. |
| `pct_fd_hi` | share of volumes above `--fd-threshold` (default 0.5 mm). |
| `dvars` | median standardised temporal derivative, as % of mean signal. Catches spikes and dropout that motion parameters miss. |
| `task_r2` | partial R² of the task regressors in the 500 most responsive voxels, over and above motion and drift. See above. |

One caveat on `tsnr_roi`: tSNR is mean ÷ temporal SD, so a region with little
signal *and* little temporal variance can score well. Read it alongside the
voxel counts from `tools/check_roi_coverage.sh` rather than on its own.

## Two things that make task_r2 honest

Both were measured on the dataset this was built against, and both matter:

- **Voxels are chosen leave-one-run-out**, from a contrast computed on the
  *other* runs of that session. Choosing them from all runs inflated `task_r2`
  by **+131%** on average, and by **23×** for individual runs — the circular
  version read 0.224 for a run whose honest score was 0.0096.
- **Motion regressors are in the nuisance model.** Leaving them out inflated it
  by a further **+62%**, because a run whose head motion happens to follow the
  task otherwise scores well for the wrong reason.

## Comparing trials

Every run appends to `run_quality.csv` with the `--trial` label as a column,
alongside `task_r2_half1` and `task_r2_half2` so the reliability can be
recomputed in R.

`--sessions` defaults to **all** sessions here, unlike the grid-consistency
tool: data quality is not affected by an intervention the way grid signal
might be.
