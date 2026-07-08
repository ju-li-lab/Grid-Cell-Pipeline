# Looking inside a run (temporal resolution for GridCAT)

GridCAT reports **one contrast value per run**. This folder adds tools to see
how that grid signal behaves *over the course of a run*, so you can chase the
observation that early runs show a lower mean contrast.

## Why one value per run is baked into GridCAT

The pipeline runs two GLMs (see `run_gridcat_analysis.m` and the GridCAT
toolbox):

1. **GLM1 (estimate).** Each grid event (`translate`, ~77/run) is modelled
   with two parametric regressors, `sin(6θ)` and `cos(6θ)`, where θ is the
   movement direction. Per voxel, the two betas give a preferred orientation
   `φ = atan2(β_sin, β_cos)/6` and an amplitude
   (`betas2ori.m`). Averaging orientation across ROI voxels gives the **mean
   grid orientation φ** (`calcMeanGridOri.m`). With `AVG_ORI_ACROSS_RUNS=1`
   this φ is pooled across all runs.

2. **GLM2 (test).** For the held-out test events (config uses odd events for
   GLM1, even for GLM2), each event gets an alignment value
   `cos(6(θ − φ))`. With `GLM2_REGRESSOR_METHOD=aligned_misaligned` the test
   events are split into **aligned** (`cos ≥ 0`) and **misaligned** (`cos < 0`)
   regressors (`generateMultiCondFile.m`), and the contrast **aligned −
   misaligned** is estimated (`createDefaultContrastSet.m`). The "mean contrast
   value" you read is that contrast image averaged over ROI voxels
   (`gridMetric_gridCodeResponseMagnitude.m`).

**The key point:** that contrast is a *single GLM beta* estimated from ~19
aligned + ~19 misaligned events pooled together. The alignment value itself is
known a priori — it is the design, not data. What actually varies and gets
estimated is the **BOLD response**. To get within-run resolution you must
estimate the response at a finer temporal grain. There is no within-run
information to "extract" from the existing per-run contrast — you have to build
a different estimate.

## The ladder of approaches (simple → rigorous)

| Approach | What it estimates | Power on small data | Effort |
|---|---|---|---|
| **1. Peristimulus ROI read-out** (this tool) | ROI-mean BOLD in a window after each event | Uses all events; noisy per-event but trends survive averaging | Low — no new GLM |
| **2. Temporal-bin GLM** | aligned/misaligned contrast per *time bin* per run | Good; each bin still pools events | Medium — edit `generateMultiCondFile` + `createDefaultContrastSet` |
| **3. Single-trial (LSS) GLM** | one deconvolved beta per event | Cleanest per-event, but single-trial betas are noisy | High — new SPM model |

Start at **1** to *look*. It is the most transparent mechanism: read the ROI
signal around each event, tag it with within-run time and alignment, and see
whether aligned−misaligned drifts across the run. Once you know the shape of
the effect, decide whether **2** (keeps GridCAT's statistics, just at finer
grain — the natural "rethink how contrasts are built") or **3** (maximum
resolution) is worth the power cost.

Because we work with small datasets, prefer **few time bins** (2 = early/late)
before going finer — every extra bin halves the events per estimate.

## The tools

### `withinrun_extract_events.m`
Produces the fundamental substrate: a long-format CSV with **one row per grid
event** — its onset, within-run fraction, orientation θ, the mean grid
orientation φ (reproduced with GridCAT's own `calcMeanGridOri`, so it matches
the pipeline), the alignment `cos(6(θ−φ))`, aligned/misaligned label, GLM1/GLM2
membership, and the ROI response.

```matlab
% path setup identical to a normal pipeline run
addpath('/Users/juli/Documents/MATLAB/GridCAT');
addpath('/Users/juli/Documents/MATLAB/GridCAT/CircStat2012a');
addpath('<your SPM dir>'); spm('defaults','fmri');

T = withinrun_extract_events( ...
      '<OUTPUT_ROOT>/sub-01s13_ses-01_GLM1', ...           % GLM1 result dir
      '<OUTPUT_ROOT>/ROI/rsub-01s13_ses-01_ErC-bilat.nii',...% same mask as GLM2
      'within_run/sub-01s13_ses-01_bilat_events.csv');
```

### `withinrun_summarise.m`
Bins the events by within-run time and reports/plots aligned, misaligned, and
their difference per bin (the within-run analogue of the per-run contrast),
with a between-run SEM.

```matlab
S = withinrun_summarise('within_run/sub-01s13_ses-01_bilat_events.csv', 'nBins', 2);
```

For a group view, pass a cell array of CSVs from several sessions.

## Caveats for the peristimulus tool

* Events are ~8 s apart, so HRFs overlap; a single event's number is
  contaminated by its neighbours. Fine for run-level *trends* (contamination
  averages out over many events), not for clean single-trial values.
* Motion and non-grid events are not regressed out (GridCAT's GLM removes
  them). The aligned−misaligned *difference* cancels most shared nuisance, so
  trust the difference more than the absolute levels.
* Validate on one session first and eyeball the ROI timeseries before trusting
  the batch.
