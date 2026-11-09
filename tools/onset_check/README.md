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

## How much to trust it

Measured on simulated data built from a real s17 event file, with AR(1) noise
and the same high-pass the GLM uses. The estimator is the one the tool uses:
normalise each run's F curve, average the curves, then interpolate the peak.

| true offset | model R² | recovered | sd over 30 repeats |
|---|---|---|---|
| 0 | 0.05 | +0.035 | 0.185 |
| 0 | 0.12 | −0.022 | 0.060 |
| +2.06 (1 volume) | 0.05 | +2.050 | 0.135 |
| +2.06 | 0.12 | +2.065 | 0.058 |
| +4.12 (2 volumes) | 0.12 | +4.139 | 0.084 |
| −2.06 | 0.12 | −2.080 | 0.054 |

Unbiased, and it separates zero from one volume by more than ten standard
deviations. A single run is *not* decisive — per-run spread is about 1.5 s at
R² = 0.12 — which is exactly why the pooled curve is the headline number and
the per-run column is only there to show consistency.

**What it cannot do.** The precision above is statistical. The *accuracy* is
limited to roughly ±1 s by the canonical HRF assumption, because the real
response in the selected voxels may peak earlier or later than the canonical
shape and the test cannot separate that from a clock offset. So this settles
whole-volume offsets beyond argument and says nothing useful about
`microtimeOnset = 8` versus `9` — a 0.129 s question that does not matter.
Read any `delta_hat` within ±1 s of zero as zero.

If `F_peak` comes back below about 5, the visual response is too weak for the
run to contribute; check that `--prefix` really points at a preprocessed series
and that `--keep` is not so small that only a handful of voxels survive.

## Running it

```bash
bash tools/onset_check/run_onset_check.sh
```

Defaults: the preprocessed tree and series prefix from `config.sh`, event
tables from `OUTPUT_ROOT/GLM_runauto/EventFiles/Event_tables`, all subjects,
all sessions, tasks 1–3, microtime onset 9 of 16, a sweep of ±5 TR in steps of
one microtime bin. Writes `onset_check.csv` (one row per run) and
`onset_check_curve.csv` (the pooled curve) to `OUTPUT_ROOT/onset_check`.

Start with one subject to check the paths resolve:

```bash
bash tools/onset_check/run_onset_check.sh --sub sub-01s17 --sessions ses-01
```

`--help` lists every option.

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
