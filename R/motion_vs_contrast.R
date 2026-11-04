# ============================================================================
#  motion_vs_contrast.R  —  Do two runs that moved apart give different
#                           contrasts?
# ----------------------------------------------------------------------------
#  Troubleshooting, not inference.
#
#  The question is pairwise, not per-run: for every pair of runs in a session,
#  how far apart did the head sit, and how far apart did the contrast come out?
#  The prediction being tested is that a subject whose between-task line is flat
#  (little repositioning) should also show similar contrasts across its runs.
#
#  An earlier version of this script plotted contrast against motion one run at
#  a time. That came out flat on every measure (all |rho| < 0.21), which is what
#  prompted the pairwise view — a run-level scatter cannot see "these two
#  disagree BECAUSE they moved apart". The old view is kept at the bottom as a
#  commented one-liner.
#
#  Needs grid_runs (import_gridcat.R) and motion_between (motion.R). Both are
#  read from disk if they are not already in the session.
#
#  HOW TO RUN (RStudio): source the whole file (Ctrl/Cmd + Shift + S). The
#  plots are wrapped in print() because source() throws away a bare ggplot.
# ============================================================================

library(dplyr)
library(readr)
library(ggplot2)


# ---- settings ---------------------------------------------------------------

HEMI <- "left"          # "bilat", "left" or "right"

# Rotations are in radians and translations in mm, so the six cannot be summed
# as they are; this is the radius that turns a rotation into a displacement.
# motion.R leaves its own value behind when sourced, and that one wins, so the
# two scripts cannot silently disagree.
if (!exists("ROT_RADIUS")) ROT_RADIUS <- 50     # mm


# ---- data -------------------------------------------------------------------

if (!exists("grid_runs"))      grid_runs      <- read_csv("grid_runs.csv",      show_col_types = FALSE)
if (!exists("motion_between")) motion_between <- read_csv("motion_between.csv", show_col_types = FALSE)


# ---- one row per run --------------------------------------------------------
# GridCAT's `run` and the task number in the rp filenames are the same thing,
# so run_num is the key. The "avg" rows are dropped: they have no run_num and
# no single head position belongs to them.

runs <- grid_runs %>%
  filter(!is_average, hemisphere == HEMI) %>%
  select(subject_id, session, task = run_num, contrast) %>%
  inner_join(select(motion_between, subject_id, session, task,
                    x, y, z, pitch, roll, yaw),
             by = c("subject_id", "session", "task"))


# ---- one row per PAIR of runs ------------------------------------------------
# Self-join inside a session, keeping task_a < task_b so each pair appears once
# and not also in reverse. Three runs per session gives three pairs: 1-2, 1-3
# and 2-3.
#
# The movement between two runs is taken as the difference of the six stored
# parameters. Those are each run's position relative to run 1, so differencing
# them is an approximation to composing the two rigid transforms properly —
# exact for the translations, first-order for the rotations. Checked against
# the composed transforms read back from the asub-*.nii headers: worst case
# 0.07 mm across all 84 pairs, mean 0.003 mm, against movements of 0.4-8 mm.
# Not worth re-reading 84 headers to recover.

pairs <- runs %>%
  inner_join(runs, by = c("subject_id", "session"),
             suffix = c("_a", "_b"), relationship = "many-to-many") %>%
  filter(task_a < task_b) %>%
  mutate(
    pair      = paste0("run ", task_a, " vs ", task_b),
    move_mm   = abs(x_a - x_b) + abs(y_a - y_b) + abs(z_a - z_b) +
                ROT_RADIUS * (abs(pitch_a - pitch_b) + abs(roll_a - roll_b) +
                              abs(yaw_a - yaw_b)),
    dcontrast = abs(contrast_a - contrast_b)
  ) %>%
  select(subject_id, session, pair, task_a, task_b, move_mm, dcontrast,
         contrast_a, contrast_b)

cat("hemisphere: ", HEMI, " — ", nrow(pairs), " pairs from ",
    nrow(distinct(pairs, subject_id, session)), " subject-sessions\n", sep = "")


# ---- 1. the plot you asked for ----------------------------------------------
# If repositioning between two runs pushed their contrasts apart, the cloud
# should slope up.

p1 <- ggplot(pairs, aes(x = move_mm, y = dcontrast)) +
  geom_point(aes(colour = pair), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linewidth = 0.5) +
  facet_wrap(~ session, labeller = label_both) +
  labs(title = paste0("Do runs that moved apart disagree? — ", HEMI, " ErC"),
       subtitle = "one point per pair of runs; line is a plain lm",
       x = "movement between the two runs (mm)",
       y = "|difference in contrast|", colour = NULL) +
  theme_minimal()
print(p1)


# ---- 2. the prediction, per subject and session ------------------------------
# "If the between-task line is flat, the contrasts should be similar too."
#
# Flatness of the line is the widest gap between any two runs in the session,
# and the contrast equivalent is the widest gap between any two of its
# contrasts. One point per subject-session, so nothing is counted three times.
# A subject in the bottom-left sat still AND agreed with itself, which is the
# prediction; the bottom-right and top-left corners are where it fails.

spread <- pairs %>%
  group_by(subject_id, session) %>%
  summarise(move_spread     = max(move_mm),
            contrast_spread = max(dcontrast), .groups = "drop")

p2 <- ggplot(spread, aes(x = move_spread, y = contrast_spread)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "grey50", linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(colour = factor(session)), size = 2.5) +
  geom_text(aes(label = sub("^sub-", "", subject_id)),
            size = 2.4, hjust = -0.15, check_overlap = TRUE) +
  labs(title = paste0("Widest spread within a session — ", HEMI, " ErC"),
       subtitle = "bottom-left = sat still and agreed with itself (the prediction)",
       x = "widest movement between any two runs (mm)",
       y = "widest contrast gap between any two runs",
       colour = "session",
       caption = paste("dashed line is least squares: with most subject-sessions",
                       "bunched at low movement,\nit follows whichever few sat",
                       "far to the right. See the drop-the-movers check below.")) +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0, colour = "grey40"))
print(p2)


# ---- 3. rough numbers --------------------------------------------------------
# Spearman, because the movement measure is skewed. Treat as a direction to
# look in: each subject contributes three pairs, so the n below is not a count
# of independent observations and no p-value here means what it normally means.

cat("\npairwise, all sessions together:\n")
cat(sprintf("  rho(movement, |contrast difference|) = %+.3f   (n = %d pairs)\n",
            cor(pairs$move_mm, pairs$dcontrast, method = "spearman"), nrow(pairs)))

pairs %>%
  group_by(session) %>%
  summarise(n = n(),
            rho = cor(move_mm, dcontrast, method = "spearman"),
            .groups = "drop") %>%
  print()

pairs %>%
  group_by(pair) %>%
  summarise(n = n(),
            rho = cor(move_mm, dcontrast, method = "spearman"),
            .groups = "drop") %>%
  print()

cat("\nper subject-session spread:\n")
cat(sprintf("  rho(movement spread, contrast spread) = %+.3f   (n = %d subject-sessions)\n",
            cor(spread$move_spread, spread$contrast_spread, method = "spearman"),
            nrow(spread)))

# Does the upward tilt in p2 survive without the handful that moved a lot?
# Spearman is already rank-based and resists them; the least-squares line in
# the plot is not, so the two can disagree sharply. If dropping a few points
# flips the sign, the tilt was those points, not a trend.
BIG_MOVE <- 4     # mm

still <- filter(spread, move_spread < BIG_MOVE)
cat(sprintf("  pearson r = %+.3f over all %d, but %+.3f over the %d under %g mm\n",
            cor(spread$move_spread, spread$contrast_spread),
            nrow(spread),
            cor(still$move_spread, still$contrast_spread),
            nrow(still), BIG_MOVE))
cat(sprintf("  (%d subject-session(s) at or above %g mm carry the difference)\n",
            nrow(spread) - nrow(still), BIG_MOVE))

# the quartile split is the same question without assuming a straight line:
# do the pairs that moved most actually disagree most?
cat("\n|contrast difference| by movement quartile:\n")
pairs %>%
  mutate(quartile = ntile(move_mm, 4)) %>%
  group_by(quartile) %>%
  summarise(movement_mm       = round(median(move_mm), 2),
            median_dcontrast  = round(median(dcontrast), 3),
            mean_dcontrast    = round(mean(dcontrast), 3),
            n = n(), .groups = "drop") %>%
  print()


# a few more things you might want:

# the same thing signed rather than absolute — does moving push the contrast
# DOWN specifically, rather than just scattering it?
# ggplot(pairs, aes(x = move_mm, y = contrast_b - contrast_a)) + geom_point() +
#   geom_hline(yintercept = 0, linetype = "dashed")

# switch hemisphere: set HEMI above and re-source. To see all three at once,
# drop the hemisphere filter from `runs`, carry it through the self-join
# (add it to the `by`), and facet on it here.

# the old per-run view: contrast against each run's own motion, no pairing
# grid_runs %>% filter(!is_average, hemisphere == HEMI) %>%
#   mutate(task = run_num) %>%
#   inner_join(readr::read_csv("motion_runs.csv"), by = c("subject_id","session","task")) %>%
#   ggplot(aes(x = mean_fd, y = contrast)) + geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x)

# which pairs disagree most, regardless of motion
# pairs %>% arrange(desc(dcontrast)) %>% head(10)

# ggsave("motion_vs_contrast.png", p1, width = 9, height = 4, dpi = 300)
