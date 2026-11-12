# ============================================================================
#  plot_voxel_stability.R — voxel grid orientation stability across runs
# ----------------------------------------------------------------------------
#  Spaghetti plot: one line per VOXEL, x = run, y = grid orientation, one
#  panel per subject x session. A voxel with a stable grid code is a flat
#  line; noise is lines jumping across the whole 0-60 deg range.
#
#  Needs voxel_ori from import_gridcat.R — run that first. The pipeline writes
#  the underlying *_voxel_ori.csv files next to the grid_metrics files.
#
#  Orientations live on a circle (60 deg for 6-fold symmetry), so each voxel's
#  line is unwrapped relative to its first run: a step from 58 to 2 deg is
#  drawn as +4, not -56. The y axis is "change from run 1" in degrees, where
#  the largest possible change is +-30 and random orientations average 15.
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

GRIDCAT_RUN <- NULL      # NULL = first imported run, or name one from run_info
HEMI        <- "bilat"   # "bilat", "left" or "right"
SESSIONS    <- NULL      # NULL = all, or e.g. c(1)
PERIOD      <- 60        # 60 for 6-fold symmetry, 90 for 4-fold

stopifnot(exists("voxel_ori"))
if (is.null(GRIDCAT_RUN)) GRIDCAT_RUN <- unique(voxel_ori$gridcat_run)[1]

d <- voxel_ori %>%
  filter(gridcat_run == GRIDCAT_RUN, hemisphere == HEMI)
if (!is.null(SESSIONS)) d <- d %>% filter(session %in% SESSIONS)
stopifnot(nrow(d) > 0)

# Signed circular difference, in (-PERIOD/2, PERIOD/2]
circ_diff <- function(a, b) ((a - b + PERIOD / 2) %% PERIOD) - PERIOD / 2

d <- d %>%
  group_by(subject_id, session, voxel) %>%
  filter(n() > 1) %>%                       # a voxel seen in one run has no line
  arrange(run, .by_group = TRUE) %>%
  mutate(dev = circ_diff(ori_deg, first(ori_deg))) %>%
  ungroup()

p <- ggplot(d, aes(x = run, y = dev, group = voxel)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_hline(yintercept = c(-15, 15), linetype = "dotted", colour = "grey60") +
  geom_line(alpha = 0.25) +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_grid(session ~ subject_id, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(d$run))) +
  scale_y_continuous(limits = c(-PERIOD / 2, PERIOD / 2)) +
  labs(title = paste0("Voxel grid orientation across runs — ", HEMI,
                      " (", GRIDCAT_RUN, ")"),
       subtitle = paste0("one line per voxel, relative to its run-1 orientation; ",
                         "dotted lines at ±15° (chance level of |Δφ|)"),
       x = "run", y = "Δ orientation vs run 1 (deg)") +
  theme_minimal() +
  theme(legend.position = "none")
print(p)

# The matching number: per subject x session, how many voxels stay within
# +-15 deg of their run-1 orientation in EVERY later run (chance ~ 50% for
# one comparison, lower for more runs).
cat("\n== voxels stable within ±15° of run 1 ===================\n")
d %>%
  filter(run > min(run)) %>%
  group_by(subject_id, session, voxel) %>%
  summarise(stable = all(abs(dev) <= PERIOD / 4), .groups = "drop_last") %>%
  summarise(n_voxels = n(), pct_stable = round(100 * mean(stable), 1),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")
