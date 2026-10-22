# Look at one slice of the grid data, per subject.
# Right now: run 1, left hemisphere — contrast and how many NaN voxels.
# Edit the filter line below to look at something else, then re-run.
#
# The plots are wrapped in print() on purpose: RStudio's Source button uses
# source(), which throws away a bare ggplot without drawing it.

library(dplyr)
library(ggplot2)

# uses grid_runs from import_gridcat.R — run that first.
# or, if you saved it to disk (SAVE_CSV <- TRUE):
# grid_runs <- readr::read_csv("grid_runs.csv", show_col_types = FALSE)


# ---- pick what to look at --------------------------------------------------

d <- grid_runs %>% filter(run != "avg", hemisphere == "bilat")

# other slices to try — swap one of these in:
# d <- grid_runs %>% filter(run == "2",   hemisphere == "left")
# d <- grid_runs %>% filter(run == "avg", hemisphere == "bilat")
# d <- grid_runs %>% filter(run == "1",   hemisphere == "left", session == 1)
# d <- grid_runs %>% filter(run != "avg", hemisphere == "left")   # all runs

d %>%
  select(subject_id, group, session, contrast, n_voxels, n_nan_voxels, pct_nan) %>%
  print(n = Inf)


# ---- contrast per subject --------------------------------------------------
# (update the title if you change the filter)

p1 <- ggplot(d, aes(x = subject_id, y = contrast, fill = factor(session))) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0) +
  labs(title = "Grid contrast — run 1, left ErC",
       x = NULL, y = "aligned vs misaligned", fill = "session") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)


# ---- NaN voxels per subject ------------------------------------------------
# how much of the mask was actually usable

p2 <- ggplot(grid_runs, aes(x = subject_id, y = pct_nan, fill = factor(session))) +
  geom_col(position = "dodge") +
  labs(title = "Checkman",
       x = NULL, y = NULL, fill = "session") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p2)


# ---- do the two go together? -----------------------------------------------
# if subjects with more dropout also have bigger contrasts, that's a worry

p3 <- ggplot(d, aes(x = pct_nan, y = contrast)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Write title here",
       x = "X Axis label", y = "Y axis label") +
  theme_minimal()
print(p3)


# ---- contrast across runs, one line per subject ----------------------------
# Spaghetti plot, one panel per session. Does a subject's contrast hold up
# from run to run, or is it carried by a single run?
#
# "avg" is excluded on purpose: it is a summary of these same points, not
# another measurement, so plotting it would bend each line towards its own
# mean. run_num is used on the x axis rather than the `run` factor because it
# is numeric, which keeps the points evenly spaced and lets the line join
# them in run order.
#
# This builds its own slice rather than using `d`, so editing the filter at
# the top of the script cannot quietly put three hemispheres on one line.

HEMI <- "bilat"          # "bilat", "left" or "right"

d_runs <- grid_runs %>% filter(!is_average, hemisphere == HEMI)

p4 <- ggplot(d_runs, aes(x = run_num, y = contrast,
                         colour = subject_id, group = subject_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ session, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(d_runs$run_num))) +
  labs(title = paste0("Grid contrast across runs — ", HEMI, " ErC"),
       x = "run", y = "aligned vs misaligned", colour = "subject") +
  theme_minimal()
print(p4)

# With many subjects the colour legend stops being readable. Two ways out:
#   drop the legend and give each subject its own small panel
# p4 + facet_grid(session ~ subject_id) + theme(legend.position = "none")
#   or keep the spaghetti and add the group mean on top
# p4 + stat_summary(aes(group = 1), fun = mean, geom = "line",
#                   colour = "black", linewidth = 1)


# a few more things you might want:

# colour by group instead of session (once subject_info.csv is filled in)
# ggplot(d, aes(x = subject_id, y = contrast, fill = group)) + geom_col()

# same plot, but one panel per session
# ggplot(d, aes(x = subject_id, y = contrast)) + geom_col() + facet_wrap(~ session)

# NaN voxels as a share of the mask — pct_nan comes straight from the
# importer now, no need to compute it here
# ggplot(d, aes(x = subject_id, y = pct_nan, fill = factor(session))) +
#   geom_col(position = "dodge")

# group means
# d %>% group_by(group) %>% summarise(mean(contrast), sd(contrast), n())

# the plots stay around as p1..p4 — type p1 to draw it again, or save it
# ggsave("contrast_run1_left.png", p1, width = 8, height = 4, dpi = 300)
