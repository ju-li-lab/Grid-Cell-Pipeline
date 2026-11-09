# grid_maps — where in the brain to look

Turns what GridCAT already wrote into voxel-wise maps. **Fits no model**, reads
no functionals, writes nothing into the GLM directories. Deleting this folder
breaks nothing.

GridCAT's ROI enters GLM2 only as the scalar φ, so every image it writes is
whole-brain. These tools read those images and nothing else — by **image
name**, never by beta number, so `keepUnusedGridEvents` and any other change
to the design's column order is irrelevant here.

```bash
bash tools/grid_maps/run_grid_maps.sh --trial uw1_sm0 --glm /path/to/GLM_output_x
```

One SLURM array task per subject, then a merge job that waits for the array and
concatenates the per-subject CSVs. `--local` runs them here, one after another.
Defaults to **ses-01 only** — sessions 2 and 3 are post-intervention, so a
cross-session number there mixes measurement reliability with a real effect.
`--help` lists the options.

---

## Read them in this order

| # | image | the question it answers |
|---|---|---|
| 1 | `roi-kept` | is there any ROI left to analyse? |
| 2 | `sens` | where could a signal be seen at all? |
| 3 | `nbr-dev` | do neighbouring voxels agree on orientation? |
| 4 | `R-runs` | does a voxel keep its orientation across runs? |
| 5 | `clust-size`, `clust-ori` | *where* do the agreeing voxels sit? |
| 6 | `coh-local`, `z-local` | the same as 3 and 5, in Rayleigh's units |
| 7 | `F-grid`, `pct-grid` | does the voxel's wobble actually fit a grid? |

Stop at the first one that is empty. If there is no ROI left (1), nothing below
it means anything.

---

## 1. `roi-kept_<hemi>_<label>.nii`

**What it is.** 0 outside the ROI, **1 = ROI voxel the GLM threw away**,
**2 = ROI voxel with an estimate**.

**Look for:** whether the 1s form one solid block or are sprinkled through the
ROI. A solid block on the ventral or anterior side is signal dropout — a
physics problem, fixed by acquisition or by moving the ROI. Scattered 1s are a
thresholding problem, fixed by `mthresh`.

**The number:** `pct_kept` in `maps.csv`. Below ~50% there is not enough ROI
left for anything downstream to be stable.

---

## 2. `sens.nii`

**What it is.** `1/sqrt(ResMS)` — the GLM's own sensitivity. High = low residual
noise = a voxel where an effect could be seen.

**Look for:** whether the ROI is dark compared with nearby cortex. If it is,
every orientation estimate inside it is noise-dominated and no amount of
analysis fixes that.

**Why this one matters more than it looks.** It is the **only map here that is
blind to direction**, so it is the only one you may select voxels on without
biasing GLM2. Do **not** select on amplitude: amplitude is the thing you are
testing.

**But it is NOT a tSNR, and on its own it is a poor selection criterion.**
`1/sqrt(ResMS)` has no mean normalisation, and noise scales with signal level,
so it largely tracks *inverse* intensity. Measured inside one GLM, in a
princeErC ROI whose mean level spans 14.6x (`thr00`, so the dim tail is still
there): `Spearman(1/sqrt(ResMS), mean level) = -0.763`, against **+0.783** for
`beta_constant/sqrt(ResMS)`. The two criteria correlate **-0.245** with each
other and share only 45% of their top half. **For selection, divide by the run
constant** (`beta_constant/sqrt(ResMS)`) or use the real tSNR from
`tools/roi_tsnr`. Use this map for looking, not for cutting.

It gets worse once `mthresh = -Inf`. In the `thr08` GLM the level range is only
1.6x and the two agree much better (+0.616, 75% overlap) — because the 0.8 cut
had already removed the dim voxels. Removing that cut widens the range and
makes the unnormalised map actively misleading.

---

## 3. `nbr-dev.nii` — the one to open first

**What it is.** For each voxel: how far its grid orientation sits from its six
face neighbours', **in orientation degrees**, 0–30. Nothing to threshold.

**How to read it, in one line:**

| value | meaning |
|---|---|
| **0°** | identical to its neighbours |
| **15°** | exactly chance — no relationship |
| **30°** | maximally opposed |

15° is also the sign-flip line: two voxels 15° apart would give GLM2 contrasts
of opposite sign.

**Look for:** patches well below 15. Check whether any of them touch the ROI.

**The trap — and it is the big one.** Smoothing makes neighbouring voxels share
signal *by construction*, so on smoothed data this map goes dark everywhere and
means nothing. The six session-1 subjects in
`GLM_output_6f_pmod_keep0_avg0_thr8` come out at **6.6–7.3° against a chance
level of 15°**, and that GLM was built on `susub-*` images — smoothed. That
number is measuring the kernel, not the brain. **Run this on an unsmoothed
trial.** Running it on both and asking whether the smoothed map is merely a
blurred copy of the unsmoothed one is itself the cleanest test of how much
smoothing is doing for you.

---

## 4. `R-runs.nii`, `R-runs-w.nii`, `n-runs.nii`

**What it is.** Per voxel, how well the runs agree on orientation. 1 = every run
points the same way, 0 = they cancel. `R-runs` uses unit vectors so one loud
run cannot carry it; `R-runs-w` is amplitude-weighted and equals
`nRuns × A_allRunsAvg / ΣA_run`. **Read `R-runs` first.**

**Look for:** voxels high in `R-runs` that are also bright in `sens`. Agreement
in a noisy voxel is luck.

**The null is not zero.** For random orientations, R averages **0.63 at n = 2**
and **0.52 at n = 3** (the tool prints this). So R = 0.5 across three runs is
*nothing*, not half an effect. Always check `n-runs` first: R from two runs and
R from three are not comparable.

**The deeper caveat.** All three runs replay the same route. A voxel can score
high because something route-locked repeats, not because it codes space. Only a
permutation null — shuffle the direction labels within run, refit GLM1,
recompute — separates those, and that is not in this tool.

---

## 5. `clust-label.nii`, `clust-size.nii`, `clust-ori.nii`

**What it is.** Join two face-adjacent voxels when their orientations agree to
within `--angle` (default 5°), then take the connected components. Unlike a
Rayleigh sphere there is **no radius to choose**. 1° and 59° count as 2° apart,
and the tool never does modular arithmetic to get that right — each voxel is a
unit vector and agreement is a dot product, so the wrap is built in.

- `clust-label` — cluster id, largest first, 0 outside.
- `clust-size` — every voxel carries its own cluster's size, so you can
  threshold at the null's 95th percentile in any viewer and see only what beats
  chance.
- `clust-ori` — every voxel carries its cluster's mean orientation, so clusters
  read as flat patches of colour. **Neighbouring patches with the same colour
  are one grid; patches with different colours are not.**

**Look for:** a compact patch, above `max_null_p95` in size, overlapping the
ROI, with `R_within` near 1.

**Read `frac_largest` and `max_clust_R` before anything else.** A real patch is
small and internally tight (`R_within` ≈ 1). A **percolating chain** is enormous
and internally loose (`R_within` ≈ 0) — it wanders through orientation space one
tolerable 5° step at a time and is a lattice artefact, not a finding. On the
smoothed session-1 data the largest cluster holds **91–94% of all valid voxels
with `R_within` = 0.07–0.45**: pure percolation. The tool prints a warning when
this happens.

**Why `--angle` defaults to 5 and not 15.** Two random orientations agree within
*t* with probability *t*/30. Each voxel has 6 neighbours, so the shuffled graph
has mean degree 6*t*/30. Bond percolation on a cubic lattice turns on at
p = 0.2488, i.e. *t* ≈ 7.5°. At the tempting *t* = 15 (the sign-flip line),
p = 0.5 and **pure noise produces one cluster the size of the brain** — extent
stops being a statistic. Do not raise `--angle` past 7.

---

## 6. `coh-local.nii`, `z-local.nii`, `n-local.nii`

**What it is.** The same question as `nbr-dev`, in the units your
`grid_metrics.txt` files already use. `coh-local` is the resultant length inside
a `--radius` sphere (default 5 mm); `z-local` is `n × coh²`, Rayleigh's own
statistic, so **≈1 is chance and >3 is roughly p < 0.05 for one neighbourhood**.

**Look for:** nothing new — use these to connect the maps above to the Rayleigh
numbers you already have. `nbr-dev` says the same thing in degrees and needs no
radius, which is why it comes first.

**Always check `n-local`.** An edge voxel with 6 valid neighbours and a deep one
with 60 are not comparable, and z grows with n.

**Same trap as 3:** meaningless on smoothed data.

---

## The two CSVs from the cluster tool

`clusters_summary.csv` — one row per subject × session:

| column | read it as |
|---|---|
| `dev_deg` | mean neighbour disagreement. **Chance is 15.** |
| `dev_null_deg`, `dev_p` | the shuffled reference and a one-sided p |
| `max_clust`, `max_null_p95` | largest cluster, and the size chance produces |
| `frac_largest` | share of valid voxels in it. **>0.25 = percolating** |
| `max_clust_R` | its internal tightness. **Near 0 = a chain, not a patch** |
| `roi_dev_deg` | `dev_deg` restricted to the ROI — the one to compare across trials |

`clusters.csv` — one row per cluster: size, mm³, mean orientation, `R_within`,
centroid in voxel **and world** coordinates (paste those into a viewer), and
`n_in_roi`.

**`dev_deg` and `dev_p` survive percolation; the size columns do not.** If
`frac_largest` is high, read the deviation columns and ignore the rest.

---

## What a pass does not mean

All of these say the orientation field is **spatially structured**, or that a
voxel is **self-consistent**. Neither is grid coding — any source of spatial
autocorrelation produces the first, and anything route-locked produces the
second. Even unsmoothed EPI carries autocorrelation from reconstruction and from
realignment interpolation, so a small effect is expected.

Whether aligned beats misaligned is still GLM2's question. These maps say
**where to look**. Having looked, test in an independently defined ROI or
against a permutation null — otherwise the forking paths have just moved into
image space.

The tools are run whole-brain on purpose. If the field is structured in visual
cortex but not in the ErC, the method works and the ErC is noisy. If it is
structured nowhere, something upstream is wrong. Both are worth knowing and
neither is visible from inside an ROI.

---

## Appendix: selecting ROI voxels on tSNR

### Unsmoothed. The gain settles it.

Smoothing raises tSNR by averaging away independent noise. At FWHM 5 mm on
2 x 2 x 3.45 mm voxels the gain is **5.4x**, equivalent to averaging **~30
voxels**:

| FWHM | tSNR gain | effective voxels averaged |
|---|---|---|
| 4 mm | x3.7 | 14 |
| **5 mm** | **x5.4** | **30** |
| 6 mm | x7.3 | 53 |
| 8 mm | x11.3 | 127 |

pmErC is about **10 voxels per hemisphere** in this grid. A kernel whose
variance-equivalent support is 30 voxels is wider than the structure, so a
smoothed "ErC voxel" is mostly not ErC. Its tSNR is then a statement about how
much tissue surrounds it — near the collateral sulcus, mostly distance to the
nearest air — and not about the voxel.

Two more reasons, either of which is enough on its own:

- Smoothing makes neighbouring tSNR values nearly identical, so inside a
  10-90 voxel ROI the ranking collapses. You end up selecting a contiguous blob
  at the ROI's centre. That may be a reasonable thing to do, but it is a
  *distance-to-edge* criterion, not a tSNR one, and should be called that.
- Keep selection and analysis in the same space. Selecting on unsmoothed tSNR
  and then analysing smoothed data means the number you selected on does not
  describe the voxel you test.

So: **select unsmoothed, analyse unsmoothed**, computed on the same
preprocessing branch immediately before smoothing, with the GLM's high-pass
(`HPF_SEC=128`, which `tsnr_report.sh` already applies - raw tSNR counts
scanner drift as noise that the GLM removes anyway).

### What to select on

| criterion | safe? | the catch |
|---|---|---|
| `beta_constant/sqrt(ResMS)` | yes | **the one to use.** A GLM-native tSNR. The grid pmods are in GLM1's model, so a strong response very slightly lowers ResMS — measured at 0.16%, and it is GLM1's half, not the half GLM2 tests. |
| `1/sqrt(ResMS)` (the raw `sens` map) | **no** | not mean-normalised, so it tracks *inverse* intensity (Spearman -0.763 with mean level at `thr00`). Use it to look, not to cut. |
| tSNR from the timeseries | yes | SD includes **task-driven** variance, so a strongly responding voxel reads *lower*. Also small, but it points the wrong way. |
| `mean(beta_constant)/sqrt(ResMS)` | yes | best of both - normalises by the voxel's own mean, so coil sensitivity drops out. Read the constants by name from `SPM.xX.name`; `keep1` shifts their numbers. |
| **amplitude** (`voxelwiseAmplitude_*`) | **NO** | direction-dependent. It is the thing GLM2 tests. Selecting on it is the circular move. |

### Order of operations, and why it is this order

`mthresh = 0.8` **is already a voxel filter, and a bad one** - it ranks voxels
by absolute intensity, which is dominated by coil sensitivity, not by noise. It
currently removes 50-87% of the ROI. So this is not about adding a filter; it
is about replacing a bad one with a principled one.

1. `mthresh = -Inf` - stop discarding voxels on brightness.
2. Whole ErC as the ROI, pmErC as secondary.
3. *Then* tSNR selection.

Measured on `sub-01s17 ses-01`: pmErC is 10 + 8 = 18 voxels, of which **6**
survive `mthresh = 0.8`. Taking the top half by tSNR of that leaves **3
voxels**. Selection needs the mask fixed first or there is nothing to select
from.

### Pick the threshold on reproducibility, never on the contrast

Sweep the percentile - keep the top 100, 75, 50, 25% - and at each step read a
**direction-blind** reliability number: run-to-run spread of phi, or median
`R-runs`. Does the orientation estimate become more repeatable as you tighten?

- **A knee** means there is a real bad tail, and it says where to cut.
- **Flat, or worse as n shrinks**, means there is nothing to select, and you
  have learned that *without looking at the contrast* - no forking path spent.

Only once the cut is fixed do you run GLM2, once. The contrast is
`cos(6*dphi)` with a sign flip at 15 degrees, so sweeping thresholds and
picking by contrast finds "significance" in noise reliably.

### The pre-post trap

Re-selecting voxels per session makes the ROI a different set of voxels in each
session. Any session-to-session change in tSNR - shim, coil loading, motion -
then masquerades as an intervention effect. Direction-blindness does not
protect against this.

- Use a **rank rule** (top k%), not an absolute threshold, so n is constant
  across subjects and sessions. An absolute cut lets n vary and weights
  subjects unequally.
- Fix the set from **session 1** and carry it forward. The session affines
  differ within a subject (verified), so this needs reslicing.
- If you do re-select per session, **report the overlap** of the selected sets.
  A low overlap means the session comparison is compromised, and that is worth
  knowing before it becomes a result.

### What is missing in the tooling

`tsnr_report.sh` computes a per-voxel tSNR image at line 81 and then deletes it
at the end of the function - it only keeps `roi_tsnr`, the ROI **mean**.
Selection needs the image. Keeping it is a two-line change (an output path plus
sparing `$tmp/tsnr.*` from the `rm`). Or skip it and use the `sens` map, which
is already written per subject x session and needs nothing new.

---

## 7. `F-grid.nii`, `pct-grid.nii`, `tsnr-glm.nii` — does the model fit?

### First, what ResMS is

A voxel's data is one number per scan: 1041 intensities across the three runs.
The GLM says that series is a weighted sum of known shapes plus a leftover:

```
y(t) = b1*rotate(t) + b2*stand(t) + b3*unused-translate(t)
     + b4*translate(t) + b5*pmodSIN(t) + b6*pmodCOS(t)
     + b7..b12*motion(t) + b13*constant + drift  +  e(t)
```

SPM picks the `b`s that make the leftover `e(t)` as small as possible.
**ResMS is the variance of that leftover** — `sum(e^2)` over the error degrees
of freedom, one number per voxel. `sqrt(ResMS)` is its SD, in intensity units:
the typical size of the wobble the model could not account for.

So ResMS is "the part of the voxel's fluctuation not explained by **the whole
model**". That includes motion, rotate/stand and the translate events — and
also the drift, **and the grid modulation itself**, because `pmodSIN` and
`pmodCOS` are two of the model's columns. The expected grid modulation has
already been subtracted before ResMS is computed, which is exactly why ResMS
cannot tell you where the grid signal is.

Measured in princeErC (run 1, 93 voxels, unsmoothed), as a share of the
variance left after the high-pass:

| component | share |
|---|---|
| **residual (ResMS)** | **91.06%** |
| motion (6 regressors) | 5.70% |
| rotate / stand | 1.04% |
| translate, unused half | 0.50% |
| translate, main effect | 0.37% |
| **grid pmods (sin + cos)** | **0.31%** |

### `F-grid` — how well the grid model fits

The grid hypothesis says the response to translating is larger at some
headings than others, with a 60 degree period. In the GLM that is the fitted
component `b_sin*pmodSIN(t) + b_cos*pmodCOS(t)`. The question "how well does
this voxel match that?" is answered by comparing the size of that component
with the voxel's own noise:

```
F = ( Sb_sin^2/sum(Vss) + Sb_cos^2/sum(Vcc) ) / ResMS / 2
```

where `Sb` are the betas summed over runs and `Vss`, `Vcc` are those columns'
parameter variances from `SPM.xX.Bcov`. **This is an F statistic on 2 and
erdf degrees of freedom**, the same thing an F-contrast over the two pmod
columns gives you in SPM.

Two things make it the right map for this question. It is **orientation-free**:
it asks whether *any* 6-fold modulation is present, not whether the grid points
a particular way, so you do not have to trust phi to read it. And each
component is divided by **its own** variance, which matters here because the
sin and cos columns are not equally estimable (variance ratio ~1.35 in runs 1
and 3).

**How to read it.** For F(2, 969): the null **median is 0.694**, the null mean
is 1.00, and p < 0.05 uncorrected is **F > 3.01**. Measured in princeErC,
sub-01s17 ses-01:

| | value |
|---|---|
| median F | **0.42** (null median 0.69) |
| 90th percentile | 1.96 |
| max | 4.47 |
| ErC voxels at p < 0.05 | **3 of 93** (4.7 expected by chance) |
| whole brain at p < 0.05 | **3.9%** (5.0% expected by chance) |

So the fit is at chance or just under it, in the ErC and everywhere else. That
is the honest headline, and it is what you would want to know before tuning
anything further.

**What it is not.** Because the orientation is fitted to the same data, a high
F in GLM1 is not evidence of grid coding — it is "a 6-fold modulation at the
orientation that best fits this data". That is precisely why the design needs
the GLM1/GLM2 split. Read `F-grid` as a troubleshooting map, never as a result,
and never select voxels on it.

### `pct-grid` — the same thing in physiological units

```
pct-grid = 100 * A_allRunsAvg / mean(beta_constant)
```

The size of the fitted grid modulation as a **percentage of the voxel's own
mean signal**. This is the map that "takes the amplitude of the voxels into
account" in the sense that matters: dividing by the run constant removes the
voxel's brightness, so a dim ErC voxel and a bright one are on the same scale.
(Within a single voxel there is nothing to weight — GridCAT's amplitude
weighting is for *combining* voxels, not for making a map.)

Measured in princeErC: **median 2.72%**, 90th percentile 5.6%, max 11.3%. Note
how implausible that is for a real BOLD effect — a genuine entorhinal task
response is a few tenths of a percent. A 2.7% "modulation" is the signature of
a quantity made of noise, and it agrees with the independent finding that
amplitude is 97-99% noise.

The two maps answer different questions and you want both: **`pct-grid` is how
big, `F-grid` is how sure.**

### `tsnr-glm` — the mean-normalised companion to `sens`

```
tsnr-glm = mean(beta_constant) / sqrt(ResMS)
```

A GLM-native tSNR. Use this, not `sens`, whenever you want a sensitivity
number: `sens` is `1/sqrt(ResMS)` with no mean normalisation and therefore
tracks inverse intensity (see section 2).
