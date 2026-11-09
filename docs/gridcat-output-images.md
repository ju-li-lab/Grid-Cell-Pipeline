# What is in a GridCAT output directory

Reference for `GLM_output_*/sub-XXsYY_ses-ZZ_GLM1/` and `..._GLM2_ErC-*/`.
Verified against `sub-01s17_ses-01` in `GLM_output_6f_pmod_keep0_avg0_thr8`
and against the GridCAT source in `~/Documents/MATLAB/GridCAT/`.

All images are **native EPI grid**: 96 x 96 x 30, 2 x 2 x 3.45 mm, float32
(except `mask.nii` uint8 and `ResMS`/`RPV` float64). One directory per
subject x session; **the affine differs between sessions** (verified: sub-01
ses-01 vs ses-02 do not match), so anything cross-session needs reslicing first.

---

## GLM1 — where orientation is estimated

36 design columns = 3 runs x (`rotate`, `stand`, `GridEvent-translate`,
`pmodSIN`, `pmodCOS`, 6 motion) + 3 run constants. With `keep0` the
`UnusedGridEvent-translate` column is dropped, so the half of the translation
events not used by GLM1 is **unmodelled** and sits in the baseline.

| file | intensity |
|---|---|
| `beta_0001 .. 0036.nii` | GLM parameter estimate for that design column, in grand-mean-scaled BOLD units. Order is the `SPM.xX.name` order above. The grid pair is **run 1: 0004 (sin) / 0005 (cos); run 2: 0015 / 0016; run 3: 0026 / 0027.** |
| `voxelwiseGridOri_translate_runN_deg.nii` | Grid orientation, **degrees in [0, 60)**. `atan2(beta_sin, beta_cos)`, wrapped to [0, 2pi), divided by 6. Circular: 59 deg and 1 deg are 2 deg apart. Never smooth, average or threshold this image with ordinary tools. |
| `voxelwiseAmplitude_translate_runN.nii` | `sqrt(beta_sin^2 + beta_cos^2)`, the depth of the 6-fold modulation, in beta units. Non-negative, **not** divided by error, so not a statistic. This is the weight used when `useVoxelWeighting = 1`. **It is a noise magnitude - see the section at the end.** |
| `voxelwiseGridOri/Amplitude_*_allRunsAvg*` | The run betas are averaged **first** (`nanmean` over runs, sin and cos separately), *then* converted. So this is a vector mean, not a mean of angles — it already is the per-voxel circular mean across runs, and its amplitude is the **resultant length**, which shrinks when runs disagree. See "stability map" below. |
| `mask.nii` | 0/1 implicit mask: voxel kept only where `Y(v,t) > mthresh * global(t)` in **every** scan of **all three** runs. Here 104,674 of 276,480. Every beta / con / ResMS voxel outside it is NaN. |
| `ResMS.nii` | Residual mean square (sigma-hat squared) per voxel, squared beta units. The GLM's own noise map and the denominator of every t. |
| `RPV.nii` | Resels per voxel — SPM's smoothness estimate for RFT correction. Eroded at the mask edge (93,073 finite here). Not useful for orientation work. |
| `SPM.mat` | Design, filter, whitening, `Vbeta` descriptors. The only reliable regressor-to-beta map. |
| `GridCAT_GLM1.mat` | The whole `gridEventType` struct: `oriPerVoxel_rad` and `amplitudePerVoxel` as full volumes for every run and for `allRunsAvg`. Everything below can be done from this one file in MATLAB. |
| `multipleConditions_runN.mat` | The `names` / `onsets` / `durations` / `pmod` handed to SPM. Where to verify the event split and the orientation values actually used. |

## GLM2 — where alignment is tested

33 design columns = 3 runs x (`rotate`, `stand`, `GridEvent-translate`,
`PMOD-alignmentWithROI1meanOri`, 6 motion) + 3 constants. The single pmod is
`cos(6 * (alpha - phi_ROI))`.

| file | intensity |
|---|---|
| `beta_0004 / 0014 / 0024.nii` | The alignment betas for runs 1 / 2 / 3. |
| `con_0001 .. 0003.nii` | Those same betas via a one-entry T contrast — numerically identical to the betas. |
| `con_0004.nii` | c = `[1 1 1]`: the **sum** across runs, not the mean. Verified: 0.184 + (-0.163) + 0.664 = 0.685 vs the reported 0.686. |
| `spmT_000N.nii` | `con / sqrt(sigma-hat^2 * c' (X'X)^-1 c)`. Note spmT is **0** outside the mask while con is **NaN** — watch this when doing image arithmetic. |
| `mask.nii` | Byte-identical to GLM1's (verified), so GLM1 orientation images and GLM2 con images can be combined voxelwise with no re-masking. |

**The enabling fact:** the ROI enters GLM2 only through `phi`, one scalar per run
per ROI. Every other part of GLM2 is estimated at **every voxel in the mask**.
So `con_0001` is already a whole-brain map of "alignment with the entorhinal
phi", and the GLM1 orientation and amplitude images are whole-brain too. The
`grid_metrics.txt` files are ROI summaries of whole-brain images that are
sitting on disk.

---

## Derived images you can build from what is already there

Work in **components**, never in angles: `u = (beta_cos, beta_sin)`. Filtering,
averaging and reslicing are all legitimate on `u` and all wrong on the degree
image.

### 1. Within-voxel stability across runs (R map)

Per voxel, with run vectors `v_r = (beta_cos,r, beta_sin,r)`:

```
R = |sum_r v_r| / sum_r |v_r|            in [0, 1]
  = 3 * Amplitude_allRunsAvg / (A_run1 + A_run2 + A_run3)
```

1 = all three runs point the same way in 6-alpha space; 0 = they cancel. Every
term is already on disk, so this costs one image-arithmetic pass and no new GLM.
This is the direct, continuous version of GridCAT's `% STABLE VOXELS`, which
only asks whether a binary angular window was met.

Two cautions. **The null is not 0**: for three vectors with uniform random
phases E[R] is about 0.5, and R is biased up when one run's amplitude dominates
— so compare against a permutation null (shuffle the direction labels within
run, refit GLM1, recompute R), not against zero. And because all three runs
replay the **same route**, a voxel can score high from route-locked artefact
rather than grid coding; only the permutation null separates those.

Also note the `allRunsAvg` vs `runN` comparisons in `grid_metrics.txt` are
circular — `allRunsAvg` contains the run it is compared to. That is where the
79-80% figures come from, and they mean nothing.

### 2. Local spatial coherence (where Rayleigh z clusters)

Rayleigh z is a property of a *set* of voxels, so to make it an image compute it
in a moving neighbourhood:

1. Unit-vector field `u_hat = (beta_cos, beta_sin) / A` — two real images.
2. Box- or Gaussian-filter each component over a small neighbourhood (3x3x3, or
   a ~6 mm sphere). Also filter the validity mask to get `n(v)`, the count of
   non-NaN voxels in each neighbourhood.
3. `R_local = sqrt(u_bar^2 + v_bar^2)`, and `z_local = n * R_local^2` — the same
   statistic `circ_rtest` computes.

Filter `(beta_cos, beta_sin)` instead of the unit vectors to get the
amplitude-weighted version, matching `useVoxelWeighting = 1`. Use the per-voxel
valid count for `n` or mask-edge voxels get inflated z.

**Do this on unsmoothed data.** On 5 mm-smoothed data neighbouring voxels share
signal by construction, so the map lights up wherever there is any signal at
all. Running it at sm0 and sm5 and asking whether the sm5 map is merely a
blurred copy of the sm0 map is itself the cleanest test of the smoothing
confound.

### 3. Cross-session stability

The same R recipe, but reslice first — the session affines differ. Reslice each
session's `(beta_cos, beta_sin)` pair into one space (the subject's T1, or
ses-01's EPI grid) with **trilinear on the components**. Then R across sessions.

A grid orientation is supposed to be stable within an environment, and the
environment is identical across your sessions, so a session-to-session R map is
the strongest within-subject test available here — and it does not touch the
group design.

### 4. Direction-blind maps for deciding where to look

- `beta_constant / sqrt(ResMS)` — a GLM-native tSNR, and the one to select on.
  **Direction-blind**, so selecting voxels on it does not bias GLM2. This is
  Stangl's tSNR suggestion computed inside the GLM you already ran.
  `Amplitude / sqrt(ResMS)` is more sensitive but is *not* direction-blind — it
  is the thing you are testing, so selecting on it biases the result.
- `1 / sqrt(ResMS)` **alone is not a tSNR.** With no mean normalisation it
  tracks *inverse* intensity, because noise scales with signal level. Measured
  inside one GLM over a princeErC ROI spanning 14.6x in mean level, it
  correlates **-0.763** with that level while `beta_constant/sqrt(ResMS)`
  correlates **+0.783**; the two share only 45% of their top half. Always
  divide by the run constant. (In the `thr08` GLM the level range is just 1.6x
  and the two agree far better, +0.616 — the 0.8 cut had already removed the
  dim voxels, so removing that cut makes the unnormalised map worse, not
  better.)
- `beta_0003 / 0013 / 0023` (`GridEvent-translate`) — the translation main
  effect. A map of these says whether the ErC responds to the task at all. The
  step-0 positive control, nearly free.
- `mask.nii` overlaid on `rsub-*_label-*pmErC_mask.nii` — shows literally which
  ErC voxels the 0.8 threshold deleted, and whether they form a contiguous
  ventral/anterior slab (susceptibility) or are scattered (noise). One figure,
  settles a lot. The resliced ROI masks are in
  `derivatives/spm-preproc*/sub-*/ses-*/anat/`.

All four are descriptive and all four are cherry-pickable. Use them to decide
where to look, then test in an independently defined ROI or against a
permutation null.

---

## Two things to fix

- **`con_0004` is a sum, not a mean.** `R/import_gridcat.R` takes
  `MEAN CON-VALUE WITHIN ROI` verbatim, so the `run == "avg"` rows in
  `grid_runs.csv` are on a 3x scale relative to runs 1-3. Correlations and
  within-`avg` tests are unaffected; magnitudes and any pooling of `avg` with
  per-run rows are not.
- **`mthresh = -Inf` admits zero-variance voxels.** With all-zero data SPM
  returns beta = 0, and `atan2(0, 0)` is 0 in MATLAB — so those voxels get a
  valid-looking orientation of exactly 0 deg and amplitude 0, which spikes an
  unweighted Rayleigh test at 0 deg. Either pass an explicit brain mask
  (`spmModel_spec.m:51` hardcodes `mask = {''}` and would need patching) or drop
  `A == 0` voxels in the metrics.

---

## Does implicit masking during smoothing help?

**What the option actually does.** SPM's Smooth module has an "Implicit masking"
flag (`spm.spatial.smooth.im`, 0/1). Its documented job is narrow: after
smoothing, it puts the *input* image's missing-data pattern back onto the
output — voxels that were 0 (integer images) or NaN (float images) before
smoothing are made missing again afterwards.

**What it does not do, and this is the part that matters.** It does **not**
exclude those zeros from the averaging. `spm_smooth` convolves the raw data,
air included. So an ErC voxel next to the ear canal still gets averaged with
~0, with or without `im = 1`. All the flag stops is signal leaking *outward*
into voxels that had no data.

So the answer to the question as asked is no: this flag does not protect the
ErC from its dark neighbours. Worth one check against your own SPM, since it is
cheap:

```bash
grep -n "job.im" $SPM_DIR/config/spm_run_smooth.m
```

**The thing you were actually reaching for** is mask-aware smoothing, also
called normalised convolution: smooth the data, smooth a binary mask with the
same kernel, divide one by the other. Each voxel then averages only real
tissue, reweighted. AFNI's `3dBlurInMask -mask` does exactly this; by hand it is
two `spm_smooth` calls and a division.

**But think about what dilution does before reaching for it.** Averaging an ErC
voxel with air multiplies `beta_sin` and `beta_cos` by the *same* factor, and
orientation is `atan2(beta_sin, beta_cos)` — a ratio. **Uniform dilution cancels
out of φ entirely.** It lowers the amplitude, which is only a weight, and it
raises noise relative to signal. It does not move the orientation.

What *does* move φ is averaging with neighbouring **tissue** that has a
different orientation, and mask-aware smoothing does nothing about that. That
is the real cost of a 5 mm kernel on a structure a few voxels across, and the
reason unsmoothed is the honest default here.

**One version would genuinely help:** smoothing *inside the ErC mask only*
(`3dBlurInMask -mask <ErC>`). Then the averaging is confined to voxels that are
supposed to share an orientation, which is a defensible thing to assume. Note
what it costs: it couples ROI voxels to each other, so the **between-voxel
Rayleigh test becomes circular** on such data. Use it for the GLM2 magnitude if
you use it at all, never for the coherence metrics.

**And it is not a fix for the masking problem.** Any of these change which
voxels clear `0.8 x global`, but that is reachable directly at `mthresh`, with
no kernel and no new assumption. At `mthresh = -Inf` smoothing should add zero
voxels; if it still adds some, the extra ones were never about signal.

---

## Amplitude is a noise magnitude, not an effect size

Under the null, `beta_sin` and `beta_cos` are both `N(0, ResMS * v_jj)` with
`v_jj` from `SPM.xX.Bcov`. The design matrix is identical for every voxel, so
`v_jj` is a constant and **sigma_beta is proportional to sqrt(ResMS)**. Then
`A = hypot(beta_sin, beta_cos)` is Rayleigh with that scale and
`E[A] = sigma_beta * sqrt(pi/2)`. Amplitude tracks the residual by
construction, with no signal involved.

Measured on `sub-01s17 ses-01`:

| quantity | runs 1 / 2 / 3 | exact null |
|---|---|---|
| `corr(A, sqrt(ResMS))` brain-wide | 0.78 / 0.46 / 0.50 | - |
| `E[A~]` | 1.215 / 1.241 / 1.238 | 1.2533 |
| `CV[A~]` | 0.61 / 0.56 / 0.60 | 0.5227 |
| **`E[A~^2]`** | **2.02 / 2.02 / 2.09** | **2.000** |

where `A~ = sqrt(beta_sin^2/v_ss + beta_cos^2/v_cc) / sqrt(ResMS)`, which is
exactly Rayleigh(1) under the null. The excess of `E[A~^2]` over 2 is the mean
squared signal-to-noise, so **the 6-fold modulation accounts for 1-4% of the
amplitude** and the rest is noise. A second cause acts brain-wide:
`global = 'None'`, so beta and sigma both scale with raw voxel intensity.

**Inside the ROI the correlation mostly vanishes.** Across princeErC,
sqrt(ResMS) spans only 1.3-1.7x while A spans 3.6-53x, so within-ROI amplitude
variation is almost entirely Rayleigh scatter (CV 0.523 by construction) rather
than residual variation - `corr(A, sqrt(ResMS))` is 0.24 / -0.23 / 0.04. Picking
low-ResMS voxels inside the ErC therefore does *not* meaningfully pick
low-amplitude ones.

### What this means for useVoxelWeighting

With weighting on, `calcMeanGridOri` computes
`atan2(sum beta_sin, sum beta_cos)` - the weight on a voxel is literally its
amplitude, hence proportional to its sigma, where inverse-variance weighting
would use `1/sigma^2`. That is the wrong way round. Amplitude weighting is a bet
that A is signal-dominated; it holds at high SNR and fails here.

It makes no measurable difference, though. Run-to-run phi spread, weighted vs
unweighted - both columns are already in every `grid_metrics.txt`, so this needs
no rerun:

| ROI | weighted | unweighted | unweighted better in |
|---|---|---|---|
| bilat | 12.99 deg | 13.00 deg | 3/6 |
| left | 10.38 deg | 11.55 deg | 2/6 |
| right | 11.45 deg | 10.42 deg | 3/6 |

A coin flip, with both at the ~11 deg null. And the top three voxels carry only
10-19% of the total weight across 25-61 voxels, so no small set of noisy voxels
is running away with phi. Set `useVoxelWeighting = 0` because it is the
defensible choice, not because it will produce a result.

## The two modulators are not equally estimable

From `SPM.xX.Bcov`, `var(beta_cos)/var(beta_sin)` is **1.30-1.38 in run 1,
1.00-1.06 in run 2, 1.34-1.46 in run 3** - nearly identical in all six
subjects, because they replay the same three routes. Run 2's route balances the
two modulators; runs 1 and 3 do not, leaving one estimated up to 21% more
precisely than the other.

So the beta noise is an ellipse, not a circle, and **the null distribution of
phi is not uniform**: simulated with the measured ratios, phi density varies
+/-7% with a 30 deg period, peaking near 0 and 30 deg, in runs 1 and 3.

This does *not* inflate the between-voxel Rayleigh test - the distortion is
purely 2-fold (axial), so the 1-fold resultant `circ_rtest` measures stays at
chance (simulated z = 0.2-1.3). But it does rule out assuming an analytic
uniform null for phi: build it by permuting direction labels and refitting.

---

## `keepUnusedGridEvents`: what it does to the design

`generateMultiCondFile.m:40-60` splits the `translate` events by
`useEventForGLM(GLMnr).logicalIds` into two **disjoint** sets:

- the events this GLM uses become a **new** column, `GridEvent-translate`, and
  that column is the only one carrying pmods (sin/cos in GLM1, the alignment
  modulator in GLM2);
- the complement stays in the original column, renamed
  `UnusedGridEvent-translate`, **with no pmod** — a plain boxcar nuisance.

`keep = 1` puts both in the design; `keep = 0` adds the unused column to
`deleteRegIds` and drops it. It is **never** one combined regressor over all
translation events. Verified on sub-01s17 ses-01:

| | columns | per run | erdf |
|---|---|---|---|
| `keep0` (`GLM_output_6f_pmod_keep0_avg0_thr8`) | 36 | rotate, stand, GridEvent-translate, pmodSIN, pmodCOS, 6 motion | 972 |
| `keep1` (`GLM_output_pmod_thr00_...`) | 39 | + `UnusedGridEvent-translate` | 969 |

The column is also dropped automatically when the GLM uses *every* translate
event, since the complement is then empty.

### What that does to ResMS

`keep = 1` projects the unused events' evoked response out of the residual, so
ResMS **falls** at translation-responsive voxels, at a cost of 1 DOF per run
(972 to 969, 0.3% — negligible).

`keep = 0` leaves that response **in** the residual. Two consequences, and the
second is the one that matters:

1. ResMS rises at translation-responsive voxels.
2. The residual now carries a task-locked, block-structured component. SPM
   estimates the AR hyperparameters (`xVi.form = 'AR(0.2)'`, measured
   `h = [0.410 0.696 | 0.321 0.719 | 0.278 0.793]`, two per session) by
   **pooling over voxels**, then applies the resulting whitening matrix to
   *every voxel in that session*. So unmodelled task structure degrades the
   whitening globally, not just where the task response lives. With a
   first-half/second-half split the unmodelled block is contiguous and ~300 s
   long, which is the worst case for that estimate.

Also, with `keep = 0` the run constant no longer estimates rest — it estimates
rest plus the unmodelled translation epochs. The pmods are mean-centred within
their parent regressor, so a constant baseline shift does not bias phi, but
every variance goes up.

`keep = 1` is correct and is what `config.sh` uses. The `thr8` folder is the
odd one out. The size of the ResMS change is **not** measurable from what is on
disk: the two folders were built from different inputs (`susub-*`, smoothed and
unwarped with no STC, vs `uasub-*`, STC and unwarped unsmoothed). One extra GLM1
on the same inputs with the flag flipped would give the number.

---

## ResMS, precisely

`ResMS = (r'r) / trRV`, where `r` is the residual after (a) session grand-mean
scaling, (b) high-pass filtering, (c) whitening, and (d) projecting out the
**full** design.

- **The pmods are in that design.** So ResMS is noise *after* the grid model,
  not before — 2 grid DOF per run removed in GLM1, 1 in GLM2. It stays
  **direction-blind**, because how much variance sin and cos remove depends on
  the amplitude, not the phase; so selecting voxels on ResMS cannot bias phi or
  the sign of the GLM2 contrast.
- **The 128 s is a frequency band, not a time window.** The HPF is implemented
  as removal of a DCT basis — 11 columns per run, 33 in all. Verified:
  1041 scans - 36 design - 33 DCT = 972 = `erdf`. So ResMS is variance between
  1/128 Hz (0.0078) and Nyquist 1/(2 x 2.06) = 0.243 Hz. Everything slower than
  128 s is already gone before the residual is formed.
- **Units are comparable across subjects and sessions.** `sGMsca = 'session
  specific'`, `GM = 100`, and `gSF` has exactly 3 unique values
  (0.3878 / 0.3914 / 0.3922) — one per run. Each session is scaled so its mean
  voxel value is 100, so `sqrt(ResMS)` is in *% of session grand mean*. This is
  what makes a ResMS map usable as a tSNR surrogate. (`global = 'None'` switches
  off *per-scan* proportional scaling only; session grand-mean scaling still
  happens.)
- **The mask threshold is per-scan and in scaled units.** `xM.I = 0` (SPM's own
  implicit-zero masking is off), `xM.VM = []` (no explicit mask), and `xM.TH`
  holds 1041 values (~79.1-79.3 here) = `mthresh x` that scan's scaled global. A
  voxel is kept only if it clears **its own scan's** threshold in **every** scan
  of **every** run.
- **One ResMS per GLM, pooled over all 1041 scans.** There is no per-run noise
  map here; per-run tSNR has to come from the timeseries.

### All-events (`specifier 8`) vs first-half + `keep1` (`specifier 2`)

Per run, the two designs are:

```
A  specifier 8:  rotate, stand, T_all, T_all(x)sin_all, T_all(x)cos_all, R1-R6   (11 cols)
B  specifier 2 + keep1:  rotate, stand, T_2nd, T_1st, T_1st(x)sin_1st, T_1st(x)cos_1st, R1-R6   (12 cols)
```

With specifier 8 the complement is empty, so `all(eventsForGridAnalysis_logIds)`
is true and `UnusedGridEvent-translate` is dropped **regardless of the keep
flag** — in design A, `keep` is inert and `keep0` and `keep1` are identical.

**ResMS should be nearly identical.** Both model every translate epoch as
task-related, so neither leaves the ~300 s contiguous block of unmodelled task
response that makes `specifier 2 + keep0` pathological. Neither design nests the
other (`T_all = T_1st + T_2nd` is in B's span, but `T_1st` alone is not in A's,
and `T_2nd(x)sin_2nd` is in neither), so they differ by what each uniquely
explains:

- A additionally fits the **6-fold modulation of second-half events**. Measured,
  that modulation is 1-4% of amplitude (`E[A~^2] = 2.02` against an exact null
  of 2.000), so it captures essentially nothing real.
- B additionally fits a **separate mean amplitude for each half** — it absorbs
  any drift in the translate response across the run, which A forces to one
  height.

ResMS is an unbiased sigma-hat-squared (`r'r / trRV`), so adding regressors that
explain nothing real leaves its expectation unchanged — `r'r` and `trRV` shrink
together. Expect a difference well under a couple of percent. If median
`sqrt(ResMS)` in the ROI moves more than that, the likely cause is a real
task-response drift across the run, which would be worth knowing about on its
own.

**They are not interchangeable for the analysis, though.** Specifier 8 puts
every event GLM2 tests into the estimate of phi, so the GLM2 contrast becomes
**circular** and will come out positive by construction. Do not read it as a
result.

Run it as a **phi instrument** instead, where no held-out contrast is needed.
phi noise scales as ~1/sqrt(n events), so doubling the events behind phi
(~38 to ~76 per run) should shrink the run-to-run spread by sqrt(2):

- partly-real-but-swamped phi → spread drops clearly below ~9 deg (from the
  current 9.4-11.4 deg against a uniform null of 11.0 deg);
- pure-noise phi → stays pinned at the null however many events it is fed.

Same prediction for session-to-session `|dphi|`: 13.5-14.3 deg now against a
15.0 deg chance line, so ~9.5-10 deg if real. That one comparison separates "not
enough data" from "no signal" — the question the seven-trial preprocessing grid
could not answer.
