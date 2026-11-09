# onset_check — does event time 0 land on the first saved volume?

**This directory is optional and self-contained.** Nothing in the pipeline
calls it. Delete `tools/onset_check/` and everything else keeps working.

## The question

Your logs open with `Time: 0;First Scanner Trigger received.`, so event time 0
is the first trigger the stimulus program received. What no log can tell you is
whether that trigger was the first *saved* volume or one of the volumes the
scanner discarded. If it was a discarded one, every onset in every GLM is late
by a whole number of TRs.

That matters more here than in an ordinary fMRI study. Successive translate
events are about 7.9 s apart and each carries a different heading, so a 2 s
offset does not merely weaken the grid estimate — it starts pairing the BOLD
response of one heading with the sine and cosine weights of the next, and the
expected result is a clean null.

## Why it does not use the grid analysis

Because the grid signal is the thing under investigation. A timing test whose
read-out is grid reproducibility cannot tell "the clock is wrong" apart from
"there is no grid", and would return a flat curve in both cases.

This uses the task's visual structure instead. The rotate / stand / translate
events tile the movement phase continuously with durations from 0.3 s to 6.4 s,
so a design built from them carries strong, irregularly spaced structure in the
band the haemodynamic response passes — and a passive walk through a virtual
environment drives visual cortex hard whatever the entorhinal cortex is doing.
Build that design, convolve it with the canonical HRF, slide it against the
data one microtime bin at a time, and record the F statistic of the task block
at every shift. The curve peaks at the true offset.

Voxels are chosen by the standard deviation of their nuisance-residualised time
course. That criterion contains no timing information, so it cannot pull the
peak toward any particular shift.

## Reading the answer

`delta_hat` is where the curve peaks, in seconds. The correction is always

```
corrected_onset = onset - delta_hat
```

| `delta_hat` | meaning |
|---|---|
| ~0 | event time 0 is the first saved volume. Nothing to do. |
| ~ +2.06 | one dummy was discarded after the trigger fired. Subtract 2.06 s from every onset. |
| ~ +4.12 | two were. |
| ~ −2.06 | the onsets sit one volume late. Add 2.06 s. |

## How much to trust it — read this before the number

Simulated from a real s17 event file with AR(1) noise and the same high-pass
the GLM uses, using the tool's own estimator (normalise each run's curve,
average the curves, interpolate the peak), 54 runs per repeat:

| model R² | median F peak | per-run sd | pooled over 54 | pooled sd |
|---|---|---|---|---|
| 0.015 | 4.9 | 6.1 s | +0.22 | **1.99 s** |
| 0.02 | 4.7 | 5.9 s | +0.58 | **1.95 s** |
| 0.05 | 4.7 | 5.0 s | −0.06 | **0.095 s** |
| 0.12 | 8.2 | 2.9 s | −0.02 | **0.053 s** |

There is a cliff between R² = 0.02 and R² = 0.05. Above it the pooled estimate
is good to a tenth of a second and settles the question outright. Below it the
pooled estimate is worthless — a 2 s standard deviation on a 2.06 s question.

**`F_peak` cannot tell you which side you are on.** At R² = 0.02 the median
F peak is 4.7; at R² = 0.05 it is also 4.7, because at low signal the peak is
just the largest of 161 noise draws. The statistic that discriminates is
`peak_ratio` — the peak divided by the median of the curve. A curve with a real
bump has `peak_ratio` well above 1.3; a flat curve sits near 1.0 wherever its
maximum happens to land. The tool prints it per run, prints the pooled value,
and **refuses to give a verdict when the pooled curve is flat** rather than
reporting a meaningless number.

## If the curve comes back flat

The first version of this tool averaged the top 2% of voxels by residual
standard deviation. On a partial slab that selects motion-edge and vessel
voxels, not task-responsive ones, and it produced flat curves with peaks
scattered across the whole search range. It now takes the mean F across every
in-brain voxel instead — selection-free, so it cannot be biased, and the
dilution from unresponsive voxels costs far less than a wrong selection did.

If it is still flat, in order of effort:

1. **Check the field of view.** The whole method assumes the task drives visual
   cortex, and this is an entorhinal-targeted slab. If occipital cortex is not
   covered, there is no strong visual response to find and no amount of
   statistics will invent one.
2. **Use the optic-flow regressors.** The raw Unity logs carry `Position` and
   `Rotation` on every frame, so angular and linear speed can be reconstructed
   at roughly 60 Hz. A continuous speed regressor is far stronger than three
   boxcars, because rotate-versus-translate is only a modest differential
   contrast — both are optic flow. This needs the logs on the cluster.
3. **Fall back to metadata.** The DICOM instance count against the protocol's
   planned measurements answers the same question without any model.

## Running it

```bash
bash tools/onset_check/run_onset_check.sh
```

Defaults: the preprocessed tree and series prefix from `config.sh`, all subjects,
all sessions, tasks 1–3, microtime onset 9 of 16, a sweep of ±5 TR in steps of
one microtime bin. Writes `onset_check.csv` (one row per run) and
`onset_check_curve.csv` (the pooled curve) to `OUTPUT_ROOT/onset_check`.

Start with one subject to check the paths resolve:

```bash
bash tools/onset_check/run_onset_check.sh --sub sub-01s17 --sessions ses-01
```

The event tables are found the way `run_gridcat_analysis.m` finds them: it
tries `GLM_runauto/EventFiles/Event_tables`, then `GLM_runauto/EventData`, then
anything holding `*EventData.txt` anywhere below `GLM_runauto`. The banner
prints the directory it settled on and how many tables are in it — check that
count matches the number of runs you expect. `--events DIR` overrides it.

`--help` lists every option. There is no voxel-selection option any more;
`--brain` only sets the in-brain intensity threshold.

## The independent confirmation

The default mode reads the slice-time-corrected series and samples the design
at `MicrotimeOnset` of `MicrotimeRes`. With `--t0 9`, 16 bins and
`slicetime.ref_fraction = 0.5` those two reference points coincide exactly, so
`delta_hat` is the event-clock offset alone — but it does assume the slice
timing correction did what it says.

`--per-slice` assumes none of that. It reads a series with **no** slice timing
correction and fits one curve per slice. Each slice's own acquisition time then
shows up in its own `delta_hat`, so

```
delta_hat(slice) = offset + SliceTiming(slice) - (MicrotimeOnset-1) * TR/MicrotimeRes
```

A regression of `delta_hat` on `SliceTiming` should return **slope 1**, which
confirms the `SliceTiming` vector and the acquisition order straight from the
data, and its intercept gives the event-clock offset without assuming the
reference fraction. One fit, three confirmations.

```bash
bash tools/onset_check/run_onset_check.sh --per-slice \
  --json /path/to/sub-01s17_ses-01_task-1_bold.json \
  --sub sub-01s17 --sessions ses-01
```

With `SLICETIME_ORDER=before` the only series without slice timing correction
is the staged raw, so the default prefix for this mode is empty and the data
are unrealigned (the motion parameters are still regressed out). If you ever
run with `SLICETIME_ORDER=after`, pass `--prefix u` instead: that series is
realigned but not slice-time corrected, which is strictly better here.
