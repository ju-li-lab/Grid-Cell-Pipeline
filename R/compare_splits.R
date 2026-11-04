# ============================================================================
#  compare_splits.R  —  Why do the two trial-splits disagree?
# ----------------------------------------------------------------------------
#  Compares two GridCAT output dirs that differ ONLY in how the grid events
#  of each run were split into GLM1 (estimate orientation) and GLM2 (test):
#     halves   = first half of events  vs  second half
#     oddeven  = odd-numbered events   vs  even-numbered
#
#  Everything else (preprocessing, ROI, pmod, no across-run averaging) is the
#  same, so any difference is caused by the split alone.
#
#  The point of the script is three numbers:
#    1. how far apart the two splits put the grid orientation (phi)
#    2. whether phi reproduces across runs / sessions, which it must if it is
#       an environment-anchored grid orientation (same env in every run)
#    3. how big the split-induced noise is next to the between-subject
#       variance you are actually trying to interpret
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

DIR_HALVES  <- "/Users/juli/Desktop/SOLIS/GLM_output_pmod_avg0_halves_oap_viarest_thr07_smooth4"
DIR_ODDEVEN <- "/Users/juli/Desktop/SOLIS/GLM_output_pmod_avg0_halves_oap_viarest_thr07_oddeven"
HEMI        <- "left"

# 6-fold symmetry => phi lives on a 60 deg circle, so two orientations can be
# at most 30 deg apart and the contrast's sign flips at 15 deg. Everything
# below is in those units.
PERIOD    <- 60
MAX_DIFF  <- PERIOD / 2   # 30 — also the mean |difference| under a uniform null... /2
NULL_MEAN <- PERIOD / 4   # 15 — mean |difference| if phi were random
FLIP      <- PERIOD / 4   # 15 — beyond this the contrast is guaranteed to flip sign

# signed difference between two orientations, wrapped into [-30, 30]
ori_diff <- function(a, b) ((a - b + MAX_DIFF) %% PERIOD) - MAX_DIFF


# ---- read the bits of a GridCAT file we need -------------------------------
#  (import_gridcat.R does the full job; this is the 3-column version, kept
#  local so this script runs on its own)

read_blocks <- function(path) {
  L <- trimws(readLines(path, warn = FALSE))
  L <- sub(";+$", "", L[nzchar(L)])
  starts <- which(startsWith(L, "GRID METRIC;"))
  ends   <- c(starts[-1] - 1L, length(L))
  out <- list()
  for (i in seq_along(starts)) {
    df  <- read.csv(text = paste(L[starts[i]:ends[i]], collapse = "\n"),
                    sep = ";", check.names = FALSE)
    tag <- df[[1]][1]
    key <- if (startsWith(tag, "Magnitude"))    "mag" else
           if (startsWith(tag, "Between"))      "bet" else
           if (startsWith(tag, "Within"))       "wit" else "ori"
    out[[key]] <- df
  }
  out
}

# GridCAT names the run three different ways depending on the block
run_label <- function(x) {
  x <- as.character(x)
  n <- str_match(x, regex("run\\s*(\\d+)", ignore_case = TRUE))[, 2]
  ifelse(!is.na(n), n,
    ifelse(str_detect(x, regex("allRuns|averaged", ignore_case = TRUE)), "avg",
      ifelse(str_detect(x, "^\\d+$"), x, NA_character_)))
}

read_dir <- function(dir, split) {
  files <- list.files(dir, pattern = "_grid_metrics\\.txt$", full.names = TRUE)
  stopifnot(length(files) > 0)
  map_dfr(files, function(p) {
    id <- str_match(basename(p),
                    "^sub-(\\d+[a-z]\\d+)_ses-(\\d+)_[A-Za-z]+-([A-Za-z]+)_")
    b  <- read_blocks(p)
    list(
      tibble(run      = run_label(b$mag[["CONTRAST NAME"]]),
             contrast = as.numeric(b$mag[["MEAN CON-VALUE WITHIN ROI"]]),
             n_voxels = as.integer(b$mag[["VOXELS WITHIN ROI"]]),
             n_nan    = as.integer(b$mag[["NaN VOXELS WITHIN ROI"]])),
      tibble(run      = run_label(b$bet[["VOXEL-WISE GRID ORI IMAGE"]]),
             rayleigh = as.numeric(b$bet[["RAYLEIGH z"]])),
      tibble(run      = run_label(b$ori[["RUN"]]),
             phi      = as.numeric(b$ori[["MEAN GRID ORI IN DEGREES (VOXELS WEIGHTED)"]]))
    ) %>%
      reduce(full_join, by = "run") %>%
      filter(!is.na(run)) %>%
      mutate(subject = id[2], session = as.integer(id[3]), hemisphere = id[4],
             split = split, usable = n_voxels - n_nan, .before = 1)
  })
}

halves  <- read_dir(DIR_HALVES,  "halves")
oddeven <- read_dir(DIR_ODDEVEN, "oddeven")
both    <- bind_rows(halves, oddeven) %>% filter(hemisphere == HEMI)

paired <- full_join(halves, oddeven,
                    by = c("subject", "session", "hemisphere", "run"),
                    suffix = c("_h", "_o")) %>%
  filter(hemisphere == HEMI, run != "avg") %>%
  mutate(d_phi      = ori_diff(phi_h, phi_o),
         will_flip  = abs(d_phi) > FLIP,
         sign_agree = sign(contrast_h) == sign(contrast_o))


# ============================================================================
#  1. The per-run table — where the two folders part company
# ============================================================================

cat("\n== per-run, ", HEMI, " pmErC ===============================\n", sep = "")
paired %>%
  transmute(subject, session, run,
            con_halves  = round(contrast_h, 3),
            con_oddeven = round(contrast_o, 3),
            phi_halves  = round(phi_h, 1),
            phi_oddeven = round(phi_o, 1),
            d_phi       = round(d_phi, 1),
            will_flip, sign_agree, usable = usable_h) %>%
  as.data.frame() %>% print(row.names = FALSE)


# ============================================================================
#  2. Sign flips are explained by the orientation disagreement
# ----------------------------------------------------------------------------
#  The GridCAT contrast is ~ A * cos(6 * (phi_test - phi_train)), so its sign
#  is decided by whether the two orientations land within 15 deg of each other.
#  If that is the mechanism, |d_phi| should predict sign agreement.
# ============================================================================

cat("\n== sign flips vs orientation disagreement =================\n")
cat("sign agreement:                 ",
    sprintf("%.0f%% (%d/%d), vs 50%% p = %.3f",
            100 * mean(paired$sign_agree), sum(paired$sign_agree), nrow(paired),
            binom.test(sum(paired$sign_agree), nrow(paired), 0.5)$p.value), "\n")
cat("runs past the 15 deg flip line: ",
    sprintf("%.0f%%", 100 * mean(paired$will_flip)), "\n")
fit <- glm(sign_agree ~ abs(d_phi), data = paired, family = binomial)
cat("P(agree) ~ |d_phi|:              slope ",
    sprintf("%.3f, p = %.4f", coef(fit)[2], summary(fit)$coefficients[2, 4]), "\n")


# ============================================================================
#  3. THE KEY TEST — on what timescale does phi reproduce?
# ----------------------------------------------------------------------------
#  The environment is identical in every run, so an environment-anchored grid
#  orientation must be the same in run 1, 2 and 3. Compare three timescales
#  against the uniform null (mean |d_phi| = 15 deg if phi is random).
# ============================================================================

pairs_within <- function(df, key) {
  df %>%
    group_by(across(all_of(key))) %>%
    filter(n() > 1) %>%
    summarise(d = list(abs(ori_diff(outer(phi, phi, "-"),
                                    0)[lower.tri(diag(n()))])),
              .groups = "drop") %>%
    unnest(d) %>% pull(d)
}

d_split   <- abs(paired$d_phi)
d_run     <- pairs_within(both %>% filter(run != "avg"),
                          c("split", "subject", "session"))
d_session <- pairs_within(both %>% filter(run == "avg"),
                          c("split", "subject"))

report <- function(label, d) {
  cat(sprintf("%-46s %5.2f deg   p vs chance = %s   %s\n",
              label, mean(d), format.pval(t.test(d, mu = NULL_MEAN)$p.value, 2),
              if (t.test(d, mu = NULL_MEAN)$p.value < 0.05 &&
                  mean(d) < NULL_MEAN) "reproduces" else "AT CHANCE"))
}

cat("\n== does the grid orientation reproduce? ===================\n")
cat("   (uniform null = 15 deg; smaller = more reproducible)\n\n")
report("across SPLITS   (same run, same data)",        d_split)
report("across RUNS     (same env, same session)",     d_run)
report("across SESSIONS (same env, weeks apart)",      d_session)


# ============================================================================
#  4. How big is the split noise next to the effect of interest?
# ============================================================================

subject_means <- paired %>%
  group_by(subject) %>%
  summarise(m = mean(c(contrast_h, contrast_o)), .groups = "drop")

cat("\n== split noise vs the variance you want to interpret =====\n")
cat(sprintf("SD(contrast_oddeven - contrast_halves)   %.3f   <- method noise\n",
            sd(paired$contrast_o - paired$contrast_h)))
cat(sprintf("SD(subject means)                        %.3f   <- the signal\n",
            sd(subject_means$m)))
cat(sprintf("ratio                                    %.1fx\n",
            sd(paired$contrast_o - paired$contrast_h) / sd(subject_means$m)))

cat("\n== is there any grid signal at all? ======================\n")
for (s in c("halves", "oddeven")) {
  x <- if (s == "halves") paired$contrast_h else paired$contrast_o
  cat(sprintf("%-8s mean contrast %+.4f   vs 0: p = %.3f\n",
              s, mean(x), t.test(x, mu = 0)$p.value))
}
# A real 6-fold signal should be clearer where the orientation is more coherent
# across voxels. A negative correlation is the signature of noise.
cat(sprintf("\n|contrast| vs between-voxel coherence:  halves r = %+.2f, oddeven r = %+.2f\n",
            cor(abs(paired$contrast_h), paired$rayleigh_h, use = "complete"),
            cor(abs(paired$contrast_o), paired$rayleigh_o, use = "complete")))
cat("(positive = coherent ROIs give bigger effects; negative = noise)\n")


# ============================================================================
#  5. Orientation jitter shrinks as the ROI grows
# ----------------------------------------------------------------------------
#  If the split disagreement is estimation noise, it should depend on how many
#  voxels the circular mean was taken over.
# ============================================================================

cat("\n== orientation jitter by ROI size =========================\n")
paired %>%
  mutate(tertile = ntile(usable_h, 3)) %>%
  group_by(tertile) %>%
  summarise(usable_voxels = paste(range(usable_h), collapse = "-"),
            mean_d_phi    = round(mean(abs(d_phi)), 2),
            n             = n(), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")
