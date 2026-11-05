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
biasing GLM2. This is Stangl's "pick high-SNR voxels" suggestion, computed
inside the GLM you already ran. Do **not** select on amplitude: amplitude is
the thing you are testing.

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
