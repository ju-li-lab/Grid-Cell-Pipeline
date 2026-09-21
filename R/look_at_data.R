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

d <- grid_runs %>% filter(run == "1", hemisphere == "bilat")

# other slices to try — swap one of these in:
# d <- grid_runs %>% filter(run == "2",   hemisphere == "left")
# d <- grid_runs %>% filter(run == "avg", hemisphere == "bilat")
# d <- grid_runs %>% filter(run == "1",   hemisphere == "left", session == 1)
# d <- grid_runs %>% filter(run != "avg", hemisphere == "left")   # all runs

d %>%
  select(subject_id, group, session, contrast, n_voxels, n_nan_voxels) %>%
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

p2 <- ggplot(d, aes(x = subject_id, y = 100 * n_nan_voxels / n_voxels, fill = factor(session))) +
  geom_col(position = "dodge") +
  labs(title = "NaN voxels in the mask — run 1, left ErC",
       x = NULL, y = NULL, fill = "session") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p2)


# ---- do the two go together? -----------------------------------------------
# if subjects with more dropout also have bigger contrasts, that's a worry

p3 <- ggplot(d, aes(x = 100 * n_nan_voxels / n_voxels, y = contrast)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Write title here",
       x = "X Axis label", y = "Y axis label") +
  theme_minimal()
print(p3)


# a few more things you might want:

# colour by group instead of session (once subject_info.csv is filled in)
# ggplot(d, aes(x = subject_id, y = contrast, fill = group)) + geom_col()

# same plot, but one panel per session
# ggplot(d, aes(x = subject_id, y = contrast)) + geom_col() + facet_wrap(~ session)

# NaN voxels as a share of the mask — fairer when masks differ a lot in size
# d %>% mutate(pct_nan = 100 * n_nan_voxels / n_voxels) %>%
#   ggplot(aes(x = subject_id, y = pct_nan, fill = factor(session))) +
#   geom_col(position = "dodge")

# group means
# d %>% group_by(group) %>% summarise(mean(contrast), sd(contrast), n())

# the plots stay around as p1, p2, p3 — type p1 to draw it again, or save it
# ggsave("contrast_run1_left.png", p1, width = 8, height = 4, dpi = 300)
