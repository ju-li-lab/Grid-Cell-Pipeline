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
| `voxelwiseAmplitude_translate_runN.nii` | `sqrt(beta_sin^2 + beta_cos^2)`, the depth of the 6-fold modulation, in beta units. Non-negative, **not** divided by error, so not a statistic. This is the weight used when `useVoxelWeighting = 1`. |
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

- `1 / sqrt(ResMS)` — a GLM-native sensitivity map. **Direction-blind**, so
  selecting voxels on it does not bias GLM2. This is Stangl's tSNR suggestion
  computed inside the GLM you already ran. `Amplitude / sqrt(ResMS)` is more
  sensitive but is *not* direction-blind — it is the thing you are testing, so
  selecting on it biases the result.
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
