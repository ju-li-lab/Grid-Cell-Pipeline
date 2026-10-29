# ============================================================================
#  Group-level summary bars, in the style of Stangl et al. (2018) Figure 2
# ----------------------------------------------------------------------------
#  Three metrics, each as a group mean +- SEM against its own chance line,
#  with the published values drawn on as a reference so you can see at a
#  glance whether your numbers are unusual or completely ordinary.
#
#  Stangl, Achtzehn, Huber, Dietrich, Tempelmann & Wolbers (2018).
#  Compromised Grid-Cell-like Representations in Old Age as a Key Mechanism
#  to Explain Age-Related Navigational Deficits.
#  Current Biology 28, 1108-1115.  doi:10.1016/j.cub.2018.02.038
#
#  That paper is the right yardstick: it uses GridCAT v1.04, written by the
#  same group, on entorhinal cortex, with the same two stability metrics.
#
#  Only subjects with a complete set of runs across all sessions are used,
#  so every bar describes the same people. See REQUIRE_COMPLETE below.
#
#  Needs grid_runs and grid_stability from import_gridcat.R - run that first.
#
#  The plots are wrapped in print() on purpose: RStudio's Source button uses
#  source(), which throws away a bare ggplot without drawing it.
# ============================================================================

library(dplyr)
library(ggplot2)

# complete_subjects() lives in its own file so that editing or restoring
# import_gridcat.R cannot take it with it. Sourced automatically.
if (!exists("complete_subjects")) {
  .cs <- c("complete_subjects.R", "R/complete_subjects.R",
           "../R/complete_subjects.R")
  .cs <- .cs[file.exists(.cs)]
  if (!length(.cs)) stop(
    "Cannot find complete_subjects.R. In RStudio use\n",
    "  Session > Set Working Directory > To Source File Location,\n",
    "then Source this script again.")
  source(.cs[1]); rm(.cs)
}


# ---- settings --------------------------------------------------------------

HEMI <- "bilat"        # "left", "right" or "bilat". Never pool hemispheres:
                      # each has its own grid orientation.

SESSIONS <- NULL      # NULL = all sessions, or e.g. c(1) for session 1 only.
                      # For the stability panels you had a reason to use
                      # session 1 alone: the intervention may act on exactly
                      # the thing these metrics measure.

GROUP_BY <- "session" # which column makes the bars. Switch to "group" once
                      # subject_info.csv has the group column filled in.

REQUIRE_COMPLETE <- TRUE  # Only keep subjects that have every run of every
                      # session for this hemisphere. A subject who ran well
                      # in session 1 but dropped out of session 3 would
                      # otherwise make session 1 look better than session 3
                      # purely by being in it. Set to FALSE to see what the
                      # incomplete subjects would have added.

SHOW_BENCHMARK <- TRUE


# ---- the benchmark ---------------------------------------------------------
#  t and p are quoted exactly from the paper. The means are read off Figure
#  2A and 2C, so treat them as approximate (about +-0.01 and +-0.5 pp).
#
#  Young adults, n = 20, healthy, 19-30 years. Four runs of 16 minutes each,
#  bilateral entorhinal cortex traced by hand on a 0.5 x 0.5 x 2 mm T2w,
#  2 x 2 x 2 mm EPI at TR 1.5 s, 5 mm FWHM smoothing, realignment only -
#  no slice timing correction, no fieldmap unwarping, no normalisation.
#
#  Old adults, n = 21, are included because they are the paper's own null:
#  a group that genuinely has no detectable grid code. That is the shape of
#  "nothing there", and it is worth seeing next to your own bars.

BENCH <- list(
  magnitude_young = 0.14,   # Fig 2A;  t(19) =  3.71, p = 0.002
  magnitude_old   = -0.04,  # Fig 2A;  t(20) = -0.79, p = 0.438   (null)
  stability_young = 53.5,   # Fig 2C;  t(19) =  4.93, p < 0.001  vs 50%
  stability_old   = 48.3    # Fig 2C;  t(20) = -1.55, p = 0.136   (null)
)

#  A caveat on the magnitude benchmark. It is a parameter estimate, so its
#  units depend on how the data were scaled. SPM applies grand mean scaling
#  during estimation, which makes betas roughly comparable between studies,
#  but only roughly. Read that line as an order of magnitude, not a target.
#  The stability benchmark is a percentage and carries no units at all, so
#  that one is a fair, direct comparison.


# ---- helpers ---------------------------------------------------------------

# Stangl's unit of analysis is the participant: one value per person, then a
# group mean, then a one-sample t-test against chance. Averaging over runs
# first is what keeps a subject with more runs from dominating.
bar_stats <- function(d, null_value) {
  d %>%
    group_by(grp) %>%
    summarise(
      n_sub = sum(!is.na(v)),
      mean  = mean(v, na.rm = TRUE),
      sem   = sd(v, na.rm = TRUE) / sqrt(n_sub),
      t     = tryCatch(unname(t.test(v, mu = null_value)$statistic),
                       error = function(e) NA_real_),
      p     = tryCatch(t.test(v, mu = null_value)$p.value,
                       error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(star = if_else(!is.na(p) & p < 0.05, "*", ""))
}

bar_plot <- function(st, null_value, ylab, title, subtitle = NULL) {
  # Star sits just outboard of the error bar, on whichever side the bar
  # points. A bar below the null gets its star below, not floating above
  # the chance line where it would look like it belonged to nothing.
  span <- diff(range(c(st$mean + st$sem, st$mean - st$sem, null_value),
                     na.rm = TRUE))
  if (!is.finite(span) || span == 0) span <- 1
  st <- st %>%
    mutate(star_y = if_else(mean >= null_value,
                            mean + sem + 0.08 * span,
                            mean - sem - 0.08 * span))
  ggplot(st, aes(x = grp, y = mean, fill = grp)) +
    geom_hline(yintercept = null_value, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15) +
    geom_text(aes(y = star_y, label = star), size = 7, vjust = 0.4) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_minimal() +
    theme(legend.position = "none")
}

# Draws the paper's young and old values as reference lines.
benchmark_lines <- function(young, old, label_x = 0.55) {
  if (!SHOW_BENCHMARK) return(NULL)
  list(
    geom_hline(yintercept = young, linetype = "dotted", colour = "#2166ac"),
    geom_hline(yintercept = old,   linetype = "dotted", colour = "#b2182b"),
    annotate("text", x = label_x, y = young, label = "Stangl young",
             hjust = 0, vjust = -0.6, size = 2.8, colour = "#2166ac"),
    annotate("text", x = label_x, y = old, label = "Stangl old (null)",
             hjust = 0, vjust = 1.4, size = 2.8, colour = "#b2182b")
  )
}

report <- function(st, null_value, what) {
  cat(sprintf("\n=== %s  (one-sample t-test vs %g) ===\n", what, null_value))
  st %>%
    mutate(mean = round(mean, 4), sem = round(sem, 4),
           t = round(t, 2), p = signif(p, 3)) %>%
    as.data.frame() %>% print(row.names = FALSE)
}


# ---- data ------------------------------------------------------------------

if (!exists("grid_runs")) stop("Run import_gridcat.R first - grid_runs is missing.")

pick <- function(d) {
  d <- d %>% filter(hemisphere == HEMI)
  if (!is.null(SESSIONS)) d <- d %>% filter(session %in% SESSIONS)
  if (!GROUP_BY %in% names(d)) stop("No column called '", GROUP_BY, "'.")
  d %>% mutate(grp = factor(.data[[GROUP_BY]]))
}

# Magnitude and spatial stability: one row per run, so drop GridCAT's own
# "avg" row - it is a summary of those same runs, not another measurement.
runs <- pick(grid_runs) %>% filter(!is_average)

# To use GridCAT's pooled estimate instead of the mean over runs:
# runs <- pick(grid_runs) %>% filter(is_average)

# Stability: one row per PAIR of runs. Pairs involving "avg" are not
# independent measurements, so they go.
stab <- NULL
if (exists("grid_stability") && !is.null(grid_stability) && nrow(grid_stability)) {
  stab <- pick(grid_stability) %>% filter(run_a != "avg", run_b != "avg")
}


# ---- completeness ----------------------------------------------------------
#  A subject contributes to every bar or to none of them. complete_subjects()
#  lives in import_gridcat.R so this script and session_effects.R apply the
#  same rule; it prints which subjects it dropped and why.

if (REQUIRE_COMPLETE) {
  keep <- complete_subjects(runs, stab)
  runs <- runs %>% filter(subject_id %in% keep)
  if (!is.null(stab)) stab <- stab %>% filter(subject_id %in% keep)
}


# ---- 1. grid code magnitude ------------------------------------------------
#  Stangl: "the average parametric modulation regressor's parameter estimates
#  within entorhinal cortex voxels". That is exactly your `contrast` column.
#  Chance is zero.

mag <- runs %>%
  group_by(subject_id, grp) %>%
  summarise(v = mean(contrast, na.rm = TRUE), .groups = "drop")

mag_st <- bar_stats(mag, 0)
report(mag_st, 0, "Grid code magnitude")

p1 <- bar_plot(mag_st, 0, "parameter estimate",
               paste0("Grid code magnitude (6-fold) - ", HEMI, " ErC"),
               "bars are group mean +- SEM; * is p < 0.05 vs zero") +
  benchmark_lines(BENCH$magnitude_young, BENCH$magnitude_old)
print(p1)


# ---- 2. temporal stability -------------------------------------------------
#  Stangl: a voxel is "stable" if its orientation changes by less than +-15
#  degrees between data halves; the metric is the proportion of such voxels,
#  averaged over runs. Grid orientations live on 0-60 degrees, so the largest
#  possible change is 30 and +-15 is exactly half of it - which is why chance
#  is 50% and not something you have to simulate.
#
#  One difference worth keeping in mind: Stangl compared the two halves of a
#  single run, whereas GridCAT's stability block compares whole runs against
#  each other. Minutes separate their two estimates; tens of minutes separate
#  yours. Yours is the harder test, so it should read lower.

if (!is.null(stab) && nrow(stab)) {

  tstab <- stab %>%
    group_by(subject_id, grp) %>%
    summarise(v = mean(pct_stable, na.rm = TRUE), .groups = "drop")

  tstab_st <- bar_stats(tstab, 50)
  report(tstab_st, 50, "Temporal stability (% stable voxels)")

  # Paper's axis is 40-60. Widen it rather than clip if anything falls
  # outside, so a group mean can never quietly disappear off the top.
  ylim2 <- range(c(40, 60, tstab_st$mean + tstab_st$sem,
                   tstab_st$mean - tstab_st$sem), na.rm = TRUE)

  p2 <- bar_plot(tstab_st, 50, "% stable voxels",
                 paste0("Temporal stability - ", HEMI, " ErC"),
                 "dashed line is the 50% chance level; * is p < 0.05 vs chance") +
    coord_cartesian(ylim = ylim2) +
    benchmark_lines(BENCH$stability_young, BENCH$stability_old)
  print(p2)

} else {
  message("No stability rows for this slice - skipping the temporal panel.")
}


# ---- 3. spatial stability --------------------------------------------------
#  Stangl's other metric: how much the voxels inside the ROI agree with each
#  other about orientation at one point in time, as Rayleigh's z for
#  non-uniformity. Higher z = more coherent.
#
#  Chance is z = 1, not zero: under the null 2z follows a chi-square with two
#  degrees of freedom, so z has an expectation of 1 whatever the voxel count.
#  z = 2.995 is the conventional p < 0.05 line.
#
#  The catch is that z = n * r^2, so it grows with the number of voxels for
#  the same true coherence. Panel 4 is there to show whether that is what you
#  are looking at.

if ("rayleigh_z" %in% names(runs)) {

  sstab <- runs %>%
    group_by(subject_id, grp) %>%
    summarise(v = mean(rayleigh_z, na.rm = TRUE), .groups = "drop")

  sstab_st <- bar_stats(sstab, 1)
  report(sstab_st, 1, "Spatial stability (Rayleigh z)")

  p3 <- bar_plot(sstab_st, 1, "Rayleigh z",
                 paste0("Spatial stability - ", HEMI, " ErC"),
                 "dashed line is chance (z = 1); dotted line is p < 0.05 (z = 2.995)") +
    geom_hline(yintercept = 2.995, linetype = "dotted", colour = "grey40")
  print(p3)

  # ---- 4. is the z just ROI size? ----
  # Same idea as the diagnostic in compare_trials.R: if the points slope, a
  # subject's apparent coherence is partly a count of their voxels.
  p4 <- runs %>%
    group_by(subject_id, grp) %>%
    summarise(z = mean(rayleigh_z, na.rm = TRUE),
              vox = mean(n_voxels - n_nan_voxels, na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(x = vox, y = z, colour = grp)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    labs(title = "Spatial stability against usable voxels",
         subtitle = "a slope means z is partly counting voxels, not coherence",
         x = "usable voxels in ROI", y = "Rayleigh z", colour = NULL) +
    theme_minimal()
  print(p4)

} else {
  message("No rayleigh_z column - skipping the spatial panels.")
}


# ---- how to read this ------------------------------------------------------

cat("
---- reading this ----

The single most useful number here is Stangl's young adult group: healthy
20-somethings, n = 20, 64 minutes of task each, whole bilateral entorhinal
cortex, running this same toolbox. They landed at 53.5% temporal stability
against a 50% chance line.

Three and a half percentage points. That is what a real, publishable,
t = 4.93 grid code looks like on this metric. It is significant because
n = 20 and the between-subject SEM is around 0.7 pp, not because the effect
is large. Their older group sat at 48.3%, which is what nothing looks like.

So a result in the high 40s to mid 50s is not a failed analysis. It is the
range the measure lives in. What separates signal from noise here is not the
distance from 50% in one ROI, it is whether the group mean sits above chance
consistently, with enough subjects to see 3 points of difference.

Your test is harder than theirs in three ways, all of which should push your
numbers down and your error bars out:
  - they compared halves of one run; GridCAT compares whole runs
  - they used whole bilateral ErC; you use an ASHS subdivision, per
    hemisphere, which is far fewer voxels per estimate
  - a smaller ROI makes every per-subject value noisier, so you need more
    subjects to resolve the same effect

Worth noting what they did NOT do, since it bears on the pipeline work: no
slice timing correction, no fieldmap unwarping, no normalisation. Realign,
smooth at 5 mm, done. Preprocessing sophistication is evidently not what
stands between you and this effect.

Finally, check the completeness table at the top. Every bar here describes
the same subjects, which is the only way a difference between sessions can
mean anything - but it also means a subject lost to one bad session is lost
to all of them. If that table is dropping people you expected to see, the
gaps are worth fixing before reading anything below it.
")
