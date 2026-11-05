# Why there is no grid signal yet, and what to do about it

Written 2026-11-04. Compares the SOLIS pipeline against Stangl et al. 2018
(*Curr Biol* 28:1108–1115), using the SOLIS100_Stahn_64ch protocol, the BIDS
JSON sidecars, `grid_runs.csv`, `tsnr_viarest.csv` and the 57 event files in
`/Volumes/train_juli/SOLIS/GLM_runauto/EventData`.

Everything marked **[measured]** comes from those files. Everything marked
**[estimate]** is arithmetic on published numbers and is labelled as such.

---

## 0. Where the analysis actually stands

From `R/grid_runs.csv` (`princeErC`, 6 subjects × 3 sessions × 3 runs,
`is_average == FALSE`): **[measured]**

| ROI | ses 1 | ses 2 | ses 3 | all |
|---|---|---|---|---|
| bilat | +0.013 (t=0.28) | −0.011 (t=−0.23) | +0.039 (t=0.73) | +0.013 (t=0.35) |
| left | −0.012 (t=−0.40) | −0.018 (t=−0.26) | −0.004 (t=−0.04) | −0.011 (t=−0.21) |
| right | +0.041 (t=0.55) | −0.087 (t=−1.50) | +0.048 (t=0.78) | +0.001 (t=0.01) |

One-sample *t* over 6 subjects, each subject's 3 runs averaged. Run-level SD
is **0.230**; the between-subject SD of the session means is about **0.09**.

Two things follow immediately.

**The measurement is noisier than the thing it measures.** A single run moves
by ±0.23 while real between-subject differences are ~0.09 wide. Averaging 3
runs gives a per-session SEM of 0.23/√3 = 0.133 — still larger than the
between-subject spread. Nothing downstream (session effects, group
comparisons) can work until this ratio is inverted.

**The point estimate is ~10× too small, and that is the real problem — not the
p-value.** Stangl's young adults: magnitude 0.14, SEM 0.038, t₁₉ = 3.71, i.e.
Cohen's *d* = 0.83. With *d* = 0.83 and n = 6 you have **38% power**
**[estimate]**. So a null at n = 6 was never surprising. But a *point estimate*
of 0.013 is not a power problem — with a working pipeline the mean should sit
near 0.14 regardless of whether it clears significance. It does not.

Power at *d* = 0.83, for later reference: n=6 → 0.38, n=10 → 0.65, n=13 →
0.78, n=20 → 0.94. If the per-run measurement noise stays at 0.23, *d* drops
to 0.65 and n=20 only reaches 0.79 **[estimate]**.

---

## 1. SOLIS vs Stangl, side by side

### Acquisition

| | Stangl 2018 | SOLIS **[measured]** |
|---|---|---|
| Scanner / coil | 3T Prisma, 64ch, GRAPPA 2 | 3T MAGNETOM XT, 64ch, GRAPPA 2 |
| Task voxel | 2 × 2 × 2 mm = **8 mm³** | 2 × 2 × 3.0 mm + 15% gap → **3.45 mm** spacing = **13.8 mm³** |
| Slices | 24 | 30 |
| Slab orientation | **parallel to the hippocampal long axis** | near-axial, tilted ~13.5° from the DICOM direction cosines |
| TR / TE / flip | 1500 / 30 / 80° | 2060 / 30 / 90° |
| PE / readout | not stated | `j-` (A≫P), TotalReadoutTime 33.25 ms, eff. echo spacing 0.35 ms |
| Run structure | 4 × 16 min = 64 min | 3 × 11.9 min (347 vol) = 35.7 min |
| Whole-brain EPI for coreg | yes, 84 slices, same 2 mm grid | not in the protocol (rest is 3 mm iso, 44 slices) |
| T2w | TSE 0.5 × 0.5 × 2.0, 64 sl, ⟂ hippocampal axis | TSE 0.4 × 0.4 × 2.0, 40 sl, C>T−25° |
| Fieldmap | none used | GRE 3 mm iso, TE 4.92/7.38 ms, once per session |

The through-plane voxel is **1.7× larger** and there is a 0.45 mm slice gap
that is never sampled. Entorhinal cortex is a thin folded ribbon; in the
through-slice direction it spans only a few voxels here.

### Task

| | Stangl 2018 | SOLIS **[measured]** |
|---|---|---|
| Movement | **active**, joystick, subject steers | **passive replay**, 3 fixed routes in a Latin square |
| Translate/rotate | mutually exclusive — straight legs only | same: `rotate` / `stand` / `translate`, no overlap |
| Objects | same 3 objects, **pre-trained to criterion** before scanning | objects change every run |
| Translation events | not stated | **77 per run**, mean 2.57 s, SOA 7.92 s |
| Translation time | **1810 ± 123 s** (young) | **198 s per run → 595 s per session** |
| Heading sampling | continuous | 36 headings at 10° steps, well counterbalanced (resultant 0.026–0.047) |

The heading counterbalancing is genuinely good and rules out sampling-bias
explanations. The **translation budget is the problem**: one SOLIS session
carries **33% of Stangl's translation data**. All three sessions together
roughly equal Stangl's single session.

### Preprocessing and GLM

| | Stangl 2018 | SOLIS |
|---|---|---|
| Realign | SPM12 | SPM25 Realign & Unwarp |
| Unwarping | not used | GRE fieldmap VDM |
| Slice timing | not mentioned (TR 1500) | yes, ref 0.5, microtimeOnset 8/16 matched |
| Normalisation | none, native space | none, native space ✔ |
| Smoothing | 5 mm FWHM | 5 mm FWHM ✔ |
| Nuisance | 6 motion params, unused grid events, feedback phase | 6 motion params, unused grid events (`keepUnusedGridEvents = 1`) ✔ |
| Short events | **translations < 1 s not modelled** | **not excluded** — 20/57 runs end in a leg < 0.5 s, 2 at exactly 0.000 s |
| ROI | **whole ErC, bilateral**, manual ITK-SNAP (Berron) | `princeErC` from ASHS Princeton (earlier runs used pmErC) |
| φ estimated from | **all 4 runs pooled**, first half of each | **per run** (`avgOriAcrossRuns = 0`) |
| Masking threshold | SPM default (0.8) on a large ROI | **0.8 on a 87-voxel ROI** |

The last three rows are where the real damage is.

---

## 2. The mechanism: three budgets, all short, all multiplying

The grid contrast is `A_test · cos(6·(φ_test − φ_train))`. The 6-fold term
means the full range of possible orientation error (0–30°) maps onto the
cosine's **entire +1 → −1 range**, with the sign flipping at **15°**. There is
no "slightly wrong" regime: a φ that is off by 16° contributes a *negative*
number to your group mean, as strongly as a φ that is off by 29°.

So noise in φ does not attenuate the effect. **It randomises its sign.** A
group mean of ~0 with a run-level SD of 0.23 is exactly what an unbiased but
noisy φ produces. That is the signature in `grid_runs.csv`.

How noisy is φ? Its precision scales roughly as `√(events × voxels)`.

**Budget 1 — usable voxels.** `grid_runs.csv` **[measured]**:

| ROI | in mask | lost to implicit mask | usable |
|---|---|---|---|
| bilat | 86.9 | 47.5 (54.1%) | **39.4** |
| left | 38.6 | 20.3 (51.7%) | 18.2 |
| right | 48.4 | 27.2 (56.2%) | 21.2 |

Stangl's whole bilateral ErC at 8 mm³ per voxel is roughly 400–500 voxels
**[estimate]**, assuming ~1700–2000 mm³ per hemisphere from the Berron
protocol. So ≈ **10× fewer usable voxels**.

**Budget 2 — the data behind φ.** `eventUsageSpecifier_GLM1 = 2` (first half
of each run's events) and `avgOriAcrossRuns = 0` means φ is estimated from
**one run's first half: ~38 events, ~99 s of translation**. Stangl pooled four
runs: ~905 s. That is **≈ 9× less data** behind the single number everything
else is conditioned on.

Combined, the orientation estimate is roughly **√(9 × 10) ≈ 9–10× noisier**
than Stangl's **[estimate]**. At that level φ is effectively uniform on
[0°, 60°), `cos(6Δφ)` is effectively uniform on [−1, +1], and the group mean
goes to zero. Which is what you see.

**Budget 3 — subjects.** n = 6 against *d* = 0.83 is 38% power.

### Why the implicit mask is the specific culprit, not dropout

From `tsnr_viarest.csv` **[measured]**:

| mask | n_vox | ROI tSNR | brain tSNR | ROI mean / (0.8 × global) |
|---|---|---|---|---|
| lpmErC | 10.4 | **81.9** | 75.4 | 1.106 |
| rpmErC | 13.4 | **74.5** | 75.4 | 1.165 |

**Entorhinal tSNR is as good as whole-brain tSNR.** There is no signal-quality
problem in this acquisition. Stangl's worry about tSNR, and his advice to
pick high-tSNR voxels, is aimed at a problem you do not have.

The last column is the one that matters. SPM keeps a voxel only if, in
**every scan of every run in the GLM**,

```
value(voxel, t) > mthresh × global(t)
```

With `mthresh = 0.8`, the ROI *mean* sits at **1.11× the cut** (range
0.92–1.43). Within-ROI intensity varies by far more than ±11% — edge voxels
partial-volume with white matter, the collateral sulcus and the ambient
cistern. So the cut passes straight **through the middle of the ROI**, which
is why 54% disappears.

Three consequences, all bad:

1. **The survivors are not a random half.** They are the brightest tip, so the
   ROI centroid shifts by several mm and you are no longer measuring the
   entorhinal cortex you segmented.
2. **It is an intersection across all three runs** — one run with slightly
   lower signal removes the voxel from all of them.
3. **Smoothing changes which voxels survive.** Smoothing pulls dim voxels up
   toward brighter neighbours, so a bigger kernel silently admits more voxels.
   If you increase smoothing and the result improves, you will not know
   whether the kernel helped or the mask did.

This is a default SPM setting interacting badly with a small ROI near the
skull base. It is not a biological finding and it is one line in
`spm_options.m`.

---

## 3. One thing to rule out before anything else

**113.8 s of each run has no events in it. [measured]**

- Scan: 347 volumes × 2.06 s = **714.8 s**
- Event file: first event at 1.01 s, last ends at **~601 s**, with inter-event
  gaps of exactly 0.00 s throughout (rotate/stand/translate tile continuously)

So the task occupies a contiguous 601 s block and 114 s is unaccounted for. If
that block sits at the **end** of the scan (scanner set to 12:01 min for a
~10 min task), everything is fine. If it sits at the **start** — a lead-in,
instructions, a re-encoding phase — then every onset in every GLM is wrong by
up to 114 s, the design matrix is fitted to unrelated volumes, and the result
is a perfect null with the exact statistics you are seeing.

A constant offset of this size reproduces *every* symptom: zero mean contrast,
run-level SD swamping the subject effect, split-half agreement explained
entirely by shared events, and orientation that never reproduces across runs.

**It is cheap to settle and nothing else is worth doing until it is.** See
step 0 below.

---

## 4. The preprocessing questions you asked, answered

### Should I use a larger smoothing kernel?

**No — and this is the one place where I think your instinct points the wrong
way.** Test it, but expect it to hurt.

The argument *for* a bigger kernel is per-voxel SNR. The argument against is
geometry. Your ROI is 87 voxels on a **2 × 2 × 3.45 mm** grid — about
1200 mm³ bilaterally, so roughly 10 × 8 × 7 mm per hemisphere. A 5 mm FWHM
Gaussian already has a kernel width comparable to the ROI's smallest
dimension. Three things follow:

1. **The ROI average stops being an average.** Averaging betas over voxels
   only buys you anything if the voxels carry partly independent noise. At
   5 mm on this grid, neighbouring voxels are already correlated by
   construction; the effective N is far below 39. Going to 8 mm makes the
   effective N approximately 1.
2. **Beyond ROI size, smoothing substitutes other tissue for yours.** A
   smoothed ROI mean is a weighted average of the raw voxels *including the
   neighbours* — hippocampus, white matter, CSF, temporal pole. You would be
   measuring the grid code of the structures next to the entorhinal cortex.
3. **It corrupts your stability metrics.** Spatial stability (Rayleigh *z*)
   measures how coherent voxel-wise orientations are across the ROI. Smoothing
   manufactures exactly that coherence. Stangl recomputed spatial stability on
   **unsmoothed** data for this reason. Your current `rayleigh_z` of 3.3–5.9,
   flat across every session and hemisphere, is consistent with smooth noise
   rather than a grid code.

Stangl used 5 mm on 2 mm isotropic data — a kernel 2.5 voxels wide. On your
grid, 5 mm is 2.5 voxels in-plane and 1.45 through-plane. You are already at
the same place in voxel units. **Keep 5 mm, and include 0 mm and 3 mm as arms
in the comparison grid** (step 4) — the unsmoothed arm is also what you need
for an honest spatial-stability number.

If your goal is per-voxel SNR, the right instrument is **growing the ROI and
recovering the masked-out voxels**, not a wider kernel. Those get you real
independent samples; smoothing gets you copies.

### Should I drop slice-timing correction?

**Keep it.** Stangl's silence is not evidence — at TR 1500 ms STC matters much
less. Your `SliceTiming` **[measured]** spans **0 → 1.9725 s of a 2.06 s TR**
(30 slices, descending interleaved). Two slices can be almost a whole TR
apart. With 2.57 s mean events at 7.92 s SOA, an uncorrected ±1 s slice offset
is a substantial fraction of the event.

The pipeline is doing this correctly: it reads per-run `SliceTiming` from the
JSON rather than assuming an order (which matters — 30 slices, *descending*
interleaved, is exactly the case where guessing goes wrong), and
`microtimeOnset = 8` of `microtimeResolution = 16` matches `ref_fraction =
0.5`. Leave it.

One stale comment to fix: `config.sh` says slices are acquired "sequentially".
They are interleaved descending. The behaviour is right; only the comment is
wrong.

### Should I drop unwarping?

**The question is not "unwarp or not" — it is "is the ROI in register with the
EPI as acquired?"** There are two internally consistent answers and you must
pick one of them, not a mixture.

- **Stangl's way:** no unwarping, and coregister T2w (+ ROIs) *into* the
  distorted EPI. The ROI inherits the distortion, so it still lands on the
  right tissue — provided the local distortion over the MTL is close to a
  rigid shift that the coregistration can absorb.
- **Your way:** unwarp the EPI, then coregister into undistorted space.
  Correct in principle, and it also fixes the non-rigid part a rigid
  coregistration cannot.

Your acquisition numbers now check out, which removes the main reason to be
nervous **[measured]**:

- `TotalReadoutTime = 0.0332499 s`, and the protocol gives
  (96 − 1) × (0.70 ms / 2) = **33.25 ms**. These agree exactly, so
  `TOTAL_READOUT_MS = 33.2499` in `spm_options.m` is right and
  `C.topup.readout_sec = 0.05409` is the wrong one (it is inert in
  `realign_unwarp` mode, but delete or fix it so it cannot mislead later).
- `PhaseEncodingDirection = j-` → `BLIP_DIRECTION = -1` is right.

That closes two of the three open unwarping questions. The remaining risk is
not the parameters but the **fieldmap's shelf life**: one GRE map at the start
of a session, applied to three 12-minute runs spanning ~40 minutes. B0 drifts
with head position and gradient heating, so run 3 is corrected with a map that
may no longer describe it.

Magnitude check: at 33.25 ms readout, a 50–100 Hz off-resonance — ordinary for
the MTL — displaces tissue by 1.7–3.3 voxels, i.e. **3.3–6.6 mm along A–P**
**[estimate]**. Entorhinal cortex is about 10 mm deep in that direction. So
this is not a detail; it is comparable to half your ROI either way.

**Do not decide this from the literature. Run both arms** (step 4) and let the
ROI-coverage count, tSNR and the positive control decide.

### Stangl's advice: pick high-tSNR voxels

Legitimate and non-circular — tSNR is blind to direction, so selecting on it
cannot manufacture a hexadirectional effect. But on your data it is **aimed at
the wrong problem**: your ErC tSNR (75–82) already equals whole-brain tSNR.

The version of this advice that *does* apply to you is its mirror image: you
are already losing 54% of the ROI to a crude intensity criterion (the 0.8
implicit mask), applied blindly. Replace it with an explicit, documented,
direction-blind criterion:

1. `maskingThreshold = -Inf` so SPM stops deleting voxels, then
2. rank ROI voxels by tSNR computed on the **rest** run or the session mean —
   never on the task runs used for GLM1/GLM2 — and keep the top *k*.

That gives you the same protection with a threshold you chose, measured on
independent data, and reported.

### Stangl's advice: focus on posterior ErC

Biologically well-founded — posteromedial EC is where hexadirectional signal
is strongest (Navarro Schröder et al. 2015) — but **note that it is not what
the paper did.** Stangl used the whole bilateral ErC, and that is the method
that produced t₁₉ = 3.71.

On your voxel grid the refinement is expensive: pmErC gives **10–13 voxels per
hemisphere before masking** and 4–7 after **[measured]**, against 39 usable
bilateral for whole `princeErC`. You cannot estimate a stable orientation from
5 voxels.

**Run whole bilateral ErC as primary** (it is both the better-powered choice
and the faithful replication) and **pmErC as a pre-registered secondary.** Once
the primary works and you know what a real effect looks like in these data,
the posterior subdivision becomes a meaningful question rather than a
power-destroying one.

---

## 5. What is not a preprocessing problem

Be clear-eyed that two things on this list cannot be fixed downstream.

**The movement is passive.** Hexadirectional BOLD has been shown repeatedly
for active virtual navigation, imagined navigation and visual exploration. A
passive optic-flow replay with no motor efference copy is a weaker
manipulation, and the published evidence base for it is thinner. Your design
has a real protection — the environment is identical across runs, so an
environment-anchored orientation *must* repeat if it exists — but you should
state the passive-movement caveat in the paper rather than discover it in
review.

**The objects change every run, so subjects are encoding, not retrieving.**
Stangl pre-trained to criterion (every object within 20 vm, eight green
smileys in a row) *before* scanning, so his subjects navigated a map they
already had. Grid-like signal is read out most reliably when a stable
allocentric reference frame is in place. Encoding a new object set each run is
a different cognitive state, and it is a plausible reason for a genuinely
smaller effect even with perfect preprocessing.

**Design facts already settled, so do not re-litigate them** (see
`R/paths.R`): subjects do not steer; all 57 runs replay 3 routes in a Latin
square; headings are counterbalanced to resultant 0.026–0.047; route and run
number are confounded within a session but balanced across sessions. There is
no behavioural or sampling account of the orientation variance left to find.

**n = 6, single-arm.** `R/subject_info.csv` has 5 rows, all cohort s17, and
`group` is empty. There is no control arm in this cohort yet, and no row at
all for `sub-06s17`. Session effects in a single arm with an unreliable
baseline are regression to the mean until proven otherwise.

---

## 6. The plan, in order

Each step has a decision rule. Do not run step *n*+1 until step *n* passes.

### Step 0 — Positive control. Does the GLM see the task at all?

Two checks, both cheap, both independent of grid cells.

**0a. Settle the 114 s.** Find out from the task/presentation logs whether the
601 s event block starts at scanner volume 1. Ask: what is the first event's
onset relative to the first saved volume, and were dummies discarded by the
scanner (`NumberOfVolumesDiscardedByScanner` is absent from the JSON) or by
the task?

**0b. Fit a throwaway GLM with `translate` vs `stand`** on one subject,
session 1, all 3 runs, in native space. Translation is continuous optic flow;
standing is not. This should light up **V1 and MT+/V5** at a trivially low
threshold.

> **Decision rule.** Strong, correctly located visual/MT+ response → timing and
> preprocessing are sound, move to step 1. Weak, absent, or displaced response
> → stop. Nothing else matters until it is fixed, and a timing offset is the
> first suspect.

This is the single highest-value hour in the whole plan. Your current result is
indistinguishable from "the design matrix is misaligned with the data", and
you cannot tell those apart from grid metrics alone.

### Step 1 — Stop deleting half the ROI

Run `matlab/check_roi_coverage.m` (it already sweeps thresholds and reproduces
SPM's per-scan global exactly) across all subject-sessions, for `princeErC`
and `pmErC`.

Then set `C.GLM.maskingThreshold = -Inf` in `spm_options.m` and re-run GridCAT.
Report `n_voxels` and `n_nan_voxels` before and after.

> **Expected:** usable bilateral voxels go from ~39 to ~87. **Decision rule:**
> if `n_nan_voxels` is not near zero afterwards, the remaining losses are true
> dropout or mask-off-brain and need the ROI reslice fixed instead.

Note the pipeline cannot currently pass an explicit mask — `spmModel_spec.m:51`
hardcodes `fmri_spec.mask = {''}`. With `-Inf` you do not need one for an ROI
analysis, but be aware SPM will then estimate air voxels too (slower, harmless).

### Step 2 — Estimate φ from all three runs

Set `C.GLM.avgOriAcrossRuns = 1`. This switches GridCAT from
`GLM1_results...run(runIdx)` to `...allRunsAvg`
(`generateMultiCondFile.m:94`), giving one orientation per session per ROI
instead of one per run.

This is Stangl's method and it triples the data behind φ: ~38 events → ~115,
99 s → 297 s. Combined with step 1 the orientation estimate should improve by
roughly **√(3 × 2.2) ≈ 2.6×** **[estimate]**.

It also removes the route/run confound: with one φ per session, the Latin
square no longer matters.

> **Decision rule:** the diagnostic is not the contrast but the *reproducibility
> of φ*. Use `tools/grid_consistency/`. φ should start to repeat across
> sessions for the same subject. If it still does not, no amount of extra data
> will help and the problem is upstream.

### Step 3 — Adopt Stangl's event rules

- **Drop translation events shorter than 1 s.** Stangl did; you do not. 20 of
  57 runs end in a leg under 0.5 s and two end at exactly 0.000 s, because the
  scanner stops mid-leg. A zero-duration event is a degenerate delta in an SPM
  design and a high-leverage point in a parametric modulation. Cost: 1 event in
  77.
- **Check that the feedback/cue phase is modelled.** Stangl included it as a
  regressor of no interest. Your event files contain only
  `rotate`/`stand`/`translate` — confirm that `stand` genuinely covers the
  cue/feedback periods and that nothing is falling into the implicit baseline.

### Step 4 — A small preprocessing grid, judged on ROI coverage and the control

Now, and only now, vary preprocessing. Keep the GLM settings from steps 1–3
fixed. Four arms:

| arm | unwarp | smoothing |
|---|---|---|
| A | yes | 5 mm (current) |
| B | **no** (Stangl) | 5 mm |
| C | yes | **0 mm** |
| D | yes | 3 mm |

Judge each on **direction-blind** criteria only — never on the grid contrast,
or you are selecting preprocessing on your effect of interest:

1. ROI voxels surviving (`check_roi_coverage.m`)
2. ROI tSNR (`tools/tsnr_report.sh`)
3. Step 0b's translate-vs-stand control: peak *t* in V1/MT+
4. Visual check that the resliced ROI sits on entorhinal cortex in the EPI

Arm C also gives you the unsmoothed data you need for an honest spatial
stability number.

While you are here, **check the ROI reslice itself**: `princeErC` resliced is
86.9 voxels × 13.8 mm³ = **1200 mm³ bilaterally**, about 600 mm³ per
hemisphere. Published whole-ErC volumes are nearer 1700–2000 mm³ per
hemisphere **[estimate]**. Compare the ASHS mask volume in native T2 space
against the resliced volume. If a third of the anatomy is vanishing in the
reslice from 0.4 × 0.4 × 2.0 mm onto a 2 × 2 × 3.45 mm grid, that is another
factor-of-two in voxel count sitting on the table.

### Step 5 — Build the right null before reading any stability number

Two of your metrics are currently compared against nulls that do not apply.

- **`pct_stable` reads ~65% against a 50% analytic chance line.** That line
  assumes two independent uniform orientation draws per voxel. They are not
  independent — the same voxel in two runs shares anatomy, vasculature and the
  smoothing kernel's local geometry, so any voxel-specific response to a
  non-grid directional feature reproduces and inflates agreement with no grid
  code present. Build an **empirical null by permuting direction labels within
  each run** and recomputing `pct_stable`.
- **GridCAT's "vs allRunsAvg" stability (79%/76%) is circular** — the average
  contains the run being compared against it. Run-vs-run is at chance
  (55.2%/46.7%). Use run-vs-run.

Also build a **split null**: 1000 random event partitions via single-trial
betas (not by re-running SPM), giving the split-induced SD and a per-run
permutation *p*. The 10° heading quantisation makes this clean — there are
only 6 distinct wedge positions, so the permutation space is small and exactly
enumerable.

### Step 6 — The reliability gate

**Do not proceed to sessions or groups until this passes.**

Compute the split-half reliability of the per-subject grid magnitude, the way
`tools/run_quality/` already does for its own metrics (it refuses to pronounce
below 10 runs and reports a CI, which is the right behaviour — reuse it).

> **Gate:** the per-session estimate must be reliable enough that
> between-subject variance exceeds within-subject measurement variance. Today
> it is 0.09 vs 0.133 and fails. Equivalently: `sd(run contrast)/√3` must drop
> below ~0.09.

The arithmetic for why this gate matters **[estimate]**: at Stangl's effect
(0.14, between-subject SD 0.17), perfect measurement gives *d* = 0.82 and 94%
power at n = 20. Your current measurement noise drops that to *d* = 0.65 and
79%. Measurement noise is not a nuisance here — it is the difference between a
study that works at n = 20 and one that does not.

### Step 7 — Then, and only then, the designed comparisons

1. **Session 1, whole bilateral ErC, one-sample t vs 0.** This is the
   replication. It must work before anything else is interpretable.
2. **Sessions 2 and 3**, same pipeline, pre-specified.
3. **Control group.** `R/subject_info.csv` has no `group` values and no row for
   `sub-06s17`. Fill it in before running any group model.
4. **Intervention effect.** With a single arm and a baseline of zero
   reliability, a session drop is regression to the mean. The control arm is
   not optional for this claim.

The best test your design offers, and one Stangl did not have, is the
**same-route repeat**: estimate φ on one session's run of route P2 and test the
aligned contrast on another session's run of P2. The design matrix is
identical, no event is shared, and the environment is the same. If a real
environment-anchored grid orientation exists, that is where it shows.

---

## 7. Summary: the ranked list

| # | Finding | Status | Effort | Expected gain |
|---|---|---|---|---|
| 1 | 114 s of each run is unaccounted for; a start-offset would produce exactly this null | **unverified** | 1 h | decisive either way |
| 2 | `maskingThreshold = 0.8` deletes 54% of the ROI; ROI mean sits 11% above the cut | **measured** | 1 line + rerun | ~2.2× voxels |
| 3 | φ estimated per run from 99 s / 38 events; Stangl pooled 4 runs | **measured** | 1 line + rerun | ~3× data for φ |
| 4 | Translations < 1 s not excluded; 2 runs have a 0.000 s event | **measured** | small | removes leverage points |
| 5 | ROI reslice may be losing a third of the anatomy | **unverified** | 1 h | possibly ~2× voxels |
| 6 | pmErC is too small at this voxel size; whole ErC should be primary | **measured** | config | ~3× voxels vs pmErC |
| 7 | Smoothing: keep 5 mm, do not increase; add 0 mm arm | **reasoned** | rerun | protects stability metrics |
| 8 | Unwarping: parameters verified correct; test both arms empirically | **measured** | rerun | unknown, bounded |
| 9 | Slice timing: keep, it is correctly implemented | **measured** | none | — |
| 10 | tSNR voxel selection: ErC tSNR already equals whole-brain; low priority | **measured** | — | small |
| 11 | Stability nulls are wrong; neither metric currently supports an inference | **measured** | analysis | correct inference |
| 12 | n = 6 gives 38% power at Stangl's effect size | **estimate** | recruitment | — |

Items 2, 3, 4 and 6 are four configuration changes. Together they should move
the orientation estimate by roughly **4–5×** **[estimate]**, which is the bulk
of the 9–10× gap to Stangl. That is the fastest path to knowing whether there
is a grid signal in these data at all.

---

## 8. What to go back to Stangl with

His two suggestions were reasonable defaults that happen not to fit this
dataset, and it is worth telling him why:

- **High-tSNR voxel selection** — entorhinal tSNR here is 75–82 against a
  whole-brain 75. There is no dropout problem to select around. The real
  voxel-loss mechanism is SPM's implicit mask, which is a different fix.
- **Posterior ErC** — at 13.8 mm³ voxels, pmErC is 10–13 voxels per hemisphere
  before masking and 4–7 after. It is the right biology at the wrong
  resolution; whole ErC has to be the primary ROI here.

The questions worth asking him instead:

1. Has he seen hexadirectional modulation under **passive** movement, and if
   so at what effect size?
2. Does he regard **pre-training to a stable object map** as necessary, or
   would encoding a new object set each run be expected to work?
3. For the orientation estimate, how much translation data does he regard as
   the practical floor? (Stangl 2018 pooled ~905 s; SOLIS has 297 s per session
   once runs are pooled, and 99 s per run as currently configured.)
