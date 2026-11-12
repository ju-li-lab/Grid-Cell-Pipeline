# ============================================================================
#  plot_gridcat_run.R — look at ONE GridCAT run
# ----------------------------------------------------------------------------
#  Needs grid_runs from import_gridcat.R — run that first.
#  Plots are wrapped in print() because source() drops bare ggplots.
# ============================================================================

library(dplyr)
library(ggplot2)

GRIDCAT_RUN <- NULL      # NULL = first imported run, or name one from run_info
HEMI        <- "bilat"   # "bilat", "left" or "right" — never pool hemispheres

stopifnot(exists("grid_runs"))
if (is.null(GRIDCAT_RUN)) GRIDCAT_RUN <- unique(grid_runs$gridcat_run)[1]

d <- grid_runs %>%
  filter(gridcat_run == GRIDCAT_RUN, !is_average, hemisphere == HEMI)
stopifnot(nrow(d) > 0)

cat("GridCAT run:", GRIDCAT_RUN, "—", HEMI, "\n")
d %>%
  select(subject_id, session, run, contrast, rayleigh_z, pct_nan) %>%
  print(n = Inf)


# ---- contrast across runs, one line per subject ------------------------------

p_contrast <- ggplot(d, aes(x = run_num, y = contrast,
                            colour = subject_id, group = subject_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ session, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(d$run_num))) +
  labs(title = paste0("Grid contrast across runs — ", HEMI, " (", GRIDCAT_RUN, ")"),
       x = "run", y = "contrast", colour = "subject") +
  theme_minimal()
print(p_contrast)


# ---- Rayleigh z across runs ---------------------------------------------------

p_rayleigh <- ggplot(d, aes(x = run_num, y = rayleigh_z,
                            colour = subject_id, group = subject_id)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ session, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(d$run_num))) +
  labs(title = paste0("Rayleigh z across runs — ", HEMI, " (", GRIDCAT_RUN, ")"),
       subtitle = "dashed line is chance (z = 1)",
       x = "run", y = "Rayleigh z", colour = "subject") +
  theme_minimal()
print(p_rayleigh)


# ---- how much of each mask was usable -----------------------------------------

p_nan <- ggplot(d, aes(x = subject_id, y = pct_nan, fill = factor(run_num))) +
  geom_col(position = "dodge") +
  facet_wrap(~ session, labeller = label_both) +
  labs(title = paste0("Unusable ROI voxels — ", HEMI, " (", GRIDCAT_RUN, ")"),
       x = NULL, y = "% NaN voxels in mask", fill = "run") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_nan)
