# What the 7 preprocessing trials show

Analysed 2026-11-05 from `/Users/juli/Desktop/SOLIS/`. Companion to
`troubleshooting-plan.md`.

## What was run

Seven of the eight cells of `unwarp {0,1} × smooth {0,5} × avg {0,1}`
(`unwarp0_smooth0_avg1` is missing), all with `thr08`, `prince`, `halves`,
`pmod`.

**Balance.** The two `unwarp0` arms are each missing 2 unilateral cells, and
different ones: `sub-05 ses-02 L+R` in `smooth0`, `sub-06 ses-02 L+R` in
`smooth5`. Bilateral survived in both. Everything below uses the **50 cells
present in all seven arms**. Worth knowing on its own: without unwarping,
some unilateral ROIs collapse below the usable-voxel floor, and which ones is
not stable.

**The variable with the largest known effect was held constant.** All seven
are `thr08`. The masking threshold costs ~54% of the ROI
(`troubleshooting-plan.md` §2) and it is not in this grid.

---

## 1. The headline

**The measurement does not survive the analysis choices.** Same subject, same
session, same run, same GLM — flip one preprocessing switch and you get a
nearly unrelated number.

Paired across the identical 54 bilateral runs:

| switch | mean \|Δcontrast\| | correlation across the two arms |
|---|---|---|
| smooth 0 → 5 (uw1, avg0) | 0.462 | **r = 0.25** |
| smooth 0 → 5 (uw0, avg0) | 0.426 | **r = 0.30** |
| unwarp 0 → 1 (sm0, avg0) | 0.600 | **r = 0.12** |
| unwarp 0 → 1 (sm5, avg0) | 0.150 | r = 0.62 |
| avg 0 → 1 (uw1, sm0) | 0.351 | r = 0.49 |
| avg 0 → 1 (uw1, sm5) | 0.154 | r = 0.55 |
| all three at once | 0.464 | **r = 0.09** |

For scale, in a single arm: SD of one run's contrast = **0.239**, SD of the
six subject means = **0.106**.

So flipping the smoothing kernel moves an individual run by **0.46** — about
**4× the entire between-subject spread** you are trying to measure, and
roughly **2× the within-arm SD**. The two corners of the grid agree at
**r = 0.09**.

Variance decomposition of the run-level contrast across all 7 arms:

| source | levels | % of variance |
|---|---|---|
| subject | 6 | 1.7% |
| session | 3 | 0.2% |
| run within session | 3 | 1.9% |
| preprocessing arm | 7 | 1.4% |
| **residual** | | **94.8%** |

Note what the 1.4% means. Preprocessing barely moves the *group mean* — it
re-randomises *which run gets which value*. That is the signature of a
quantity with no stable content: you are not biasing an estimate, you are
redrawing it.

---

## 2. Smoothing's gain is threshold evasion, and here is the proof

ROI coverage, bilateral, common cells:

| arm | mask voxels | lost to mthresh | usable | % lost |
|---|---|---|---|---|
| uw0 sm0 | 85.4 | 71.1 | **14.3** | 83.3 |
| uw1 sm0 | 86.9 | 62.8 | **24.2** | 72.2 |
| uw0 sm5 | 85.4 | 60.2 | **25.2** | 70.5 |
| uw1 sm5 | 86.9 | 47.5 | **39.4** | 54.6 |

Read naively: both switches help and they stack — smoothing adds ~11–15
voxels, unwarping adds ~10–14.

But smoothing cannot create data. It adds voxels because it **pulls dim ROI
voxels above the `0.8 × global` cut by mixing in brighter neighbours**. If
those 15 extra voxels carried real signal, something downstream would improve.
Nothing does:

| | uw1 sm0 | uw1 sm5 |
|---|---|---|
| usable voxels | 24.2 | **39.4** (+63%) |
| φ run-to-run spread (deg) | 11.35 | 10.92 |
| φ session-to-session \|Δφ\| (deg) | **13.46** | 14.33 (worse) |
| run-vs-run temporal stability | 49.4% | 51.8% |

**+63% voxels, no improvement in anything.** That is what threshold evasion
looks like: the added voxels are weighted averages of voxels already in the
set.

### The Rayleigh z result is the clearest finding in the whole grid

Between-voxel orientation coherence within the ROI:

| arm | Rayleigh z |
|---|---|
| uw0 sm0 | **1.30** |
| uw1 sm0 | **1.15** |
| uw0 sm5 | 6.13 |
| uw1 sm5 | 7.39 |

**Smoothing multiplies spatial coherence by 5–6×.** Nothing about the brain
changed; a Gaussian kernel made neighbouring voxels report the same number.

The unsmoothed value is the honest one, and z ≈ 1.2 (p ≈ 0.3) says **there is
no coherent grid orientation across entorhinal voxels at all**. Every
spatial-stability number computed on the 5 mm data — including the ones in
`grid_runs.csv` — is an artifact of the kernel. This is the measurement
Stangl deliberately recomputed unsmoothed, and now you can see why.

**Falsifiable prediction, worth one run:** at `mthresh = -Inf`, smoothing
should add **zero** voxels. If `smooth0` and `smooth5` then report the same
count, the entire apparent benefit of smoothing in this grid was threshold
evasion and nothing else.

---

## 3. Unwarping does something real, but only visibly at 0 mm

Unwarping is the one switch with a defensible mechanism here: it improves
ROI–EPI registration, so the mask lands on entorhinal cortex rather than on
its edge. Evidence that it is doing something real rather than cosmetic:

- At **sm0** it adds 10 voxels with **no** smoothing confound (83.3% → 72.2%
  lost). That gain cannot be threshold evasion.
- It gives the **best session-to-session φ reproducibility** in the grid
  (13.46° at sm0, against a 15.0° chance line) — the only number anywhere
  that moves in the right direction for the right reason.
- Without it, unilateral ROIs drop out of the analysis entirely, and
  unpredictably.

But note the interaction: at **sm0** the unwarp switch moves a run by 0.60
(r = 0.12), at **sm5** by 0.15 (r = 0.62). Smoothing is *hiding* the
instability, not removing it. Agreement bought by blurring is not agreement.

**Keep unwarping on.** It is the one variable in this grid that earns its
place.

---

## 4. φ never reproduces, in any arm

This is the criterion that matters, and it is the reason the grid is a dead
end.

| arm | run-to-run spread (null 11.0°) | session-to-session \|Δφ\| (chance 15.0°) |
|---|---|---|
| uw0 sm0 | 11.36 | 15.91 |
| uw1 sm0 | 11.35 | **13.46** |
| uw0 sm5 | **9.43** | 14.18 |
| uw1 sm5 | 10.92 | 14.33 |

**The entire 2 × 2 preprocessing grid spans 1.9° of φ reproducibility**, on a
measure whose chance line is 11.0°. The best arm (9.43°) beats chance by less
than the spread of the grid itself — and its advantage comes from smoothing,
which §2 shows is manufactured.

Temporal stability tells the same story. Honest run-vs-run: 46.6 / 49.4 / 51.8
/ 54.7% against a 50% chance line. The circular "vs allRunsAvg" number reads
**79–80% in every single arm** — it is insensitive to preprocessing because
the comparison contains the run being compared, so it measures nothing.

**You have hit the floor of what preprocessing can do.** The environment is
identical across runs, so an environment-anchored grid orientation *must*
repeat. It does not, in any configuration you have tried.

---

## 5. Do not pick an arm by its contrast

Across 21 arm × hemisphere combinations, **16 of 21 point estimates are
negative**, and the three largest |t| values are all negative:

- `uw1 sm0 avg1` bilat: **−0.1375, t = −3.35**
- `uw1 sm5 avg1` bilat: −0.0354, t = −2.84
- `uw0 sm0 avg0` left: −0.1321, t = −2.37

A negative grid contrast means aligned directions produced *less* BOLD than
misaligned — the opposite of a grid code. At the individual-run level the sign
split is near chance in every arm (42–69% negative, only one arm at p < 0.05
out of 21 tests), so this is not yet a demonstrated systematic bias. But the
direction is the wrong one everywhere, and it rules out the one thing you
might be tempted to do.

**If you had selected the arm with the strongest effect, you would have
selected a significant anti-grid result.** With 21 non-independent tests, a
measurement that reshuffles at r = 0.09 across arms, and no direction-blind
criterion separating them, arm selection on the contrast is pure
garden-of-forking-paths. Judge arms on ROI coverage, tSNR, φ reproducibility
and the positive control. Never on the thing you are trying to measure.

One genuine signal in the noise: **`avg1` makes the subject means more
negative and more significant** (`uw1 sm0`: −0.083, t = −0.81 → −0.138,
t = −3.35). That is consistent with `avg1` working as designed — pooling runs
reduces noise in φ, so whatever systematic offset exists stops being masked.
Worth following up as a *bias* question, not as a result: if there is a
reproducible Δφ between the first-half and second-half event sets, the most
likely cause is that the halves are not equally balanced for heading. The
route is identical in every run and subject, so any such imbalance would be
identical too — and therefore systematic. Check
`corr(sin6θ, cos6θ)` and the heading histogram **within each half** of each
route, not just over the whole run.

---

## 6. What to test instead

In rough order of expected return:

1. **`maskingThreshold = -Inf`.** The variable this grid held constant.
   Expected: usable bilateral voxels 39 → ~87 at sm5, and 24 → ~87 at sm0,
   which also neutralises the smoothing confound in every metric above.
2. **The step-0 timing check.** Still unverified, still capable of explaining
   everything on its own. 114 s of each run has no events in it.
3. **The split**: `halves` vs `oddeven` vs `Q1+Q3 / Q2+Q4`, and the
   leakage-free `estimate on run 1, test on runs 2–3`. These are *analysis*
   variables — cheap, and the `halves` split is the prime suspect for §5's
   negative bias.
4. **ROI size.** 86.9 voxels × 13.8 mm³ = 1200 mm³ bilaterally, against
   published whole-ErC volumes near 1700–2000 mm³ *per hemisphere*. Compare
   the ASHS mask in native T2 space against the resliced mask. Another factor
   of two may be recoverable.
5. **Short-event exclusion** (< 1 s, Stangl's rule).
6. **An empirical null** by permuting direction labels within run, via
   single-trial betas rather than re-running SPM.

---

## 7. How to run far fewer trials

Six rules, in order of how much they save.

**1. Put the variable with the biggest known effect in the grid.** Seven
trials went to two variables worth 1.9° of φ reproducibility between them,
while `mthresh` — worth a factor of two in usable voxels — sat at the SPM
default in all of them. Rank candidates by expected effect *before* running
anything, and spend the budget at the top.

**2. Never vary a GLM-only switch as a preprocessing trial.** `avg0/avg1`
changes nothing before GLM2; GLM1 is bit-identical. Four of your seven folders
are two pairs differing only in `avg`. One GridCAT run can emit both. Same for
`xFold` (your archived 5f/7f control models), `gridRegressorMethod` and the
event split. **That alone would have cut this grid from 7 runs to 4.**

**3. Screen on direction-blind metrics before running GLM2 at all.** A ladder,
not a factorial:

| stage | needs | cost | decides |
|---|---|---|---|
| A | preprocessing output only | minutes | ROI coverage, ROI tSNR, visual ROI–EPI check |
| B | + GLM1 | ~1 run | Rayleigh z, φ run-to-run reproducibility |
| C | + GLM2 | full | the contrast |

Stage A alone would have ranked these four preprocessing arms correctly and
cost no GridCAT runs. Only arms that survive B deserve C. `check_roi_coverage.m`
already does stage A and sweeps thresholds in one pass — for `mthresh` you do
not need separate trials at all, just one sweep.

**4. Hold smoothing at 0 while diagnosing anything.** Smoothing confounds
every metric in this document — it inflates voxel counts through the mask, it
multiplies Rayleigh z by 5–6×, and it masks instability in other switches
(unwarp r = 0.12 at sm0 vs 0.62 at sm5). Diagnose unsmoothed, then add
smoothing once at the end as the *last* variable, if at all.

**5. Vary one thing from a fixed reference; reserve factorials for suspected
interactions.** You found a real one (unwarp × smooth), so the 2 × 2 was not
wasted — but 2ᵏ grows fast and most pairs do not interact. Change one switch
from a frozen baseline, and only expand to a factorial where a mechanism
predicts the interaction.

**6. Always keep the cell set balanced, and run the full set or none.** Two
arms here are missing different unilateral cells, so any unrestricted
comparison across arms is partly a comparison of which subjects are in it.
Intersect to the common set before reading anything.

### What the next round should look like

Three runs, not seven:

| run | config | question |
|---|---|---|
| 1 | `uw1 sm0 thr-Inf halves` | does removing the mask recover the ROI, and does smoothing's gain vanish? |
| 2 | `uw1 sm0 thr-Inf oddeven` | is the halves split producing the negative bias? |
| 3 | `uw1 sm5 thr-Inf halves` | does smoothing still do anything once the mask is gone? |

`avg0` and `avg1` come out of all three for free. Stage A (ROI coverage, tSNR)
runs before any of them.
