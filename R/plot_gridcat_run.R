# ============================================================================
#  plot_gridcat_run.R — look at ONE GridCAT run
# ----------------------------------------------------------------------------
#  Needs grid_runs from import_gridcat.R — run that first. The last section
#  tests the contrast against 0 per hemisphere and session.
#  Plots are wrapped in print() because source() drops bare ggplots.
# ============================================================================

library(dplyr)
library(ggplot2)

GRIDCAT_RUN <- c("results_solis_oap")      # NULL = first imported run, or name one from run_info
HEMI        <- "right"   # "bilat", "left" or "right" — never pool hemispheres

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
# print(p_nan)


# ---- is the contrast different from zero? ---------------------------------------
#  The group-level check that the analysis works as intended: one value per
#  subject and session (the mean of its runs — the "avg" row is a sum),
#  t-tested against 0 for every hemisphere and session, and per group once
#  subject_info.csv has groups.

subject_contrast <- grid_runs %>%
  filter(gridcat_run == GRIDCAT_RUN, !is_average) %>%
  group_by(subject_id, group, hemisphere, session) %>%
  summarise(contrast = mean(contrast, na.rm = TRUE), .groups = "drop")
by_group <- any(!is.na(subject_contrast$group))

contrast_vs_zero <- subject_contrast %>%
  group_by(across(all_of(c("hemisphere", "session", if (by_group) "group")))) %>%
  summarise(n = n(), mean = mean(contrast), sd = sd(contrast), .groups = "drop") %>%
  mutate(t    = mean / (sd / sqrt(n)),
         p    = 2 * pt(-abs(t), n - 1),
         ci95 = qt(0.975, n - 1) * sd / sqrt(n))

cat("\n== contrast against 0 (two-sided one-sample t-test) =====================\n")
contrast_vs_zero %>%
  mutate(across(c(mean, sd, t, ci95), ~ round(.x, 3)), p = signif(p, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)

last_session <- subject_contrast %>%
  group_by(subject_id, hemisphere) %>%
  filter(session == max(session))

p_zero <- ggplot(subject_contrast, aes(x = factor(session), y = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(aes(group = subject_id), colour = "grey70") +
  geom_point(colour = "grey70", size = 1.8) +
  geom_text(data = last_session, aes(label = sub("^sub-", "", subject_id)),
            hjust = -0.3, size = 2.6, colour = "#52514e", check_overlap = TRUE) +
  geom_pointrange(data = contrast_vs_zero,
                  aes(y = mean, ymin = mean - ci95, ymax = mean + ci95),
                  colour = "#0b0b0b", size = 0.5) +
  geom_text(data = contrast_vs_zero, aes(y = Inf, label = sprintf("p = %.2g", p)),
            vjust = 1.5, size = 3) +
  facet_grid(rows = if (by_group) vars(group), cols = vars(hemisphere)) +
  labs(title = paste0("Is the contrast different from 0? (", GRIDCAT_RUN, ")"),
       subtitle = "grey: one line per subject (mean of its runs); black: mean ± 95% CI",
       x = "session", y = "contrast") +
  theme_minimal()
print(p_zero)
