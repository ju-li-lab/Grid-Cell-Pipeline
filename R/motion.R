# ============================================================================
#  motion.R  —  How much did each subject move, within and between tasks?
# ----------------------------------------------------------------------------
#  Written to check whether the outliers in the grid metrics are really
#  head motion. Two separate questions, two separate tables:
#
#    motion_runs     one row per subject x session x task
#                    within-task motion: mean/max framewise displacement
#    motion_between  one row per subject x session x task
#                    where the task sat relative to task 1 of the same session
#
#  HOW TO RUN (RStudio): edit PREPROC_DIR below, then source the whole file
#  (Ctrl/Cmd + Shift + S). The plots are wrapped in print() because source()
#  throws away a bare ggplot without drawing it.
# ============================================================================


# ============================================================================
#  SETTINGS  —  edit these
# ============================================================================

PREPROC_DIR <- "/Volumes/train_juli/SOLIS/derivatives/spm-preproc-oap-viarest/"

# Rotations are in radians, translations in mm, so they cannot be summed as
# they are. The usual fix (Power et al. 2012) is to convert a rotation into
# the distance it moves a point on the cortex, i.e. radians x brain radius.
ROT_RADIUS <- 50     # mm

# A volume counts as bad above this framewise displacement.
FD_THRESH <- 0.5     # mm

SAVE_CSV   <- TRUE
OUTPUT_DIR <- "."


# ============================================================================
#  PACKAGES
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(ggplot2)


# ============================================================================
#  1. Find the realignment parameter files
# ============================================================================
#  SPM's Realign & Unwarp writes one rp_*.txt per run, six columns: three
#  translations (mm) and three rotations (rad), each volume relative to the
#  FIRST VOLUME OF ITS OWN RUN. That last point is the whole reason this
#  script has two halves: the rp files start at zero for every run, so they
#  say nothing at all about how the runs sat relative to each other.

rp_files <- list.files(PREPROC_DIR, pattern = "^rp_a.*\\.txt$",
                       recursive = TRUE, full.names = TRUE)

if (length(rp_files) == 0) {
  stop("No rp_*.txt files under '", PREPROC_DIR, "'. Is PREPROC_DIR right?")
}

# sub-01s17_ses-01_task-3_run-2_bold -> subject, session, task
#                       ^^^^^^ the extra run- tag on re-acquired runs is kept
#                       out of the task number on purpose
parse_rp <- function(path) {
  base <- basename(path)
  m <- str_match(base, "^rp_a(sub-[^_]+)_ses-(\\d+)_task-(\\d+)(.*)_bold\\.txt$")
  if (is.na(m[1, 1])) return(NULL)
  tibble(
    subject_id = m[1, 2],
    session    = as.integer(m[1, 3]),
    task       = as.integer(m[1, 4]),
    rp_path    = path,
    # the realigned image the parameters belong to — needed in part 3
    img_path   = file.path(dirname(path),
                           str_replace(base, "^rp_", "") %>%
                             str_replace("\\.txt$", ".nii"))
  )
}

runs <- map(rp_files, parse_rp) %>% list_rbind() %>%
  arrange(subject_id, session, task)

cat("Found ", nrow(runs), " runs from ",
    n_distinct(runs$subject_id), " subjects.\n", sep = "")


# ============================================================================
#  2. Within-task motion
# ============================================================================
#  Framewise displacement: how far the head moved from one volume to the
#  next. Summed absolute differences rather than the vector length, which is
#  what Power et al. use and what most QC thresholds are quoted against.
#
#  Mean FD and max FD answer different questions. A high mean FD is a subject
#  who never sat still; a high max FD with a low mean is one sudden jolt,
#  which is the pattern that produces a single odd run.

read_rp <- function(path) {
  x <- read_table(path, col_names = c("x", "y", "z", "pitch", "roll", "yaw"),
                  col_types = "dddddd")
  x %>% mutate(volume = row_number())
}

motion_vols <- runs %>%
  mutate(rp = map(rp_path, read_rp)) %>%
  select(subject_id, session, task, rp) %>%
  unnest(rp) %>%
  group_by(subject_id, session, task) %>%
  mutate(
    # rotations expressed in mm so the six can be added up
    fd = abs(x - lag(x)) + abs(y - lag(y)) + abs(z - lag(z)) +
         ROT_RADIUS * (abs(pitch - lag(pitch)) + abs(roll - lag(roll)) +
                       abs(yaw - lag(yaw)))
  ) %>%
  ungroup()

motion_runs <- motion_vols %>%
  group_by(subject_id, session, task) %>%
  summarise(
    n_volumes  = n(),
    mean_fd    = mean(fd, na.rm = TRUE),
    max_fd     = max(fd, na.rm = TRUE),
    n_bad      = sum(fd > FD_THRESH, na.rm = TRUE),
    pct_bad    = 100 * n_bad / (n_volumes - 1),
    # drift: the span covered over the run, not frame to frame
    range_trans = max(sqrt(x^2 + y^2 + z^2)) - min(sqrt(x^2 + y^2 + z^2)),
    max_trans   = max(abs(c(x, y, z))),
    max_rot_deg = max(abs(c(pitch, roll, yaw))) * 180 / pi,
    .groups = "drop"
  )

print(motion_runs, n = Inf)


# ============================================================================
#  3. Between-task motion
# ============================================================================
#  Realign & Unwarp ran once per session with one data block per task, so all
#  three tasks were aligned to the first volume of task 1. That alignment is
#  written into the header of each realigned image (the asub-*.nii files):
#  task 1's header is untouched, and task 2's and 3's carry the rigid
#  transform that brought them onto task 1. So the between-task motion is the
#  difference between those headers.
#
#  Sessions were realigned separately, so these numbers are only comparable
#  WITHIN a session. There is no task-1-of-session-1 common reference here.

# The 4x4 voxel-to-world matrix out of a NIfTI-1 header (srow_x/y/z, offset
# 280, twelve float32s), converted to SPM's 1-based voxel indexing so the
# relative transforms come out in the same units as the rp files.
nifti_mat <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  sizeof_hdr <- readBin(con, "integer", n = 1, size = 4)
  endian <- if (sizeof_hdr == 348L) .Platform$endian else
    if (.Platform$endian == "little") "big" else "little"
  seek(con, 280)
  srow <- readBin(con, "double", n = 12, size = 4, endian = endian)
  M <- rbind(matrix(srow, nrow = 3, byrow = TRUE), c(0, 0, 0, 1))
  S <- diag(4); S[1:3, 4] <- -1       # NIfTI counts voxels from 0, SPM from 1
  M %*% S
}

# spm_imatrix for a rigid transform: pull the three translations and three
# rotations back out of a 4x4. Validated against SPM's own rp output.
rigid_params <- function(M) {
  clamp <- function(v) min(max(v, -1), 1)
  R <- M[1:3, 1:3]
  roll  <- asin(clamp(R[1, 3]))
  c_r   <- cos(roll)
  pitch <- atan2(clamp(R[2, 3] / c_r), clamp(R[3, 3] / c_r))
  yaw   <- atan2(clamp(R[1, 2] / c_r), clamp(R[1, 1] / c_r))
  c(M[1:3, 4], pitch, roll, yaw)
}

between_one_session <- function(df) {
  mats <- map(df$img_path, nifti_mat)
  ref  <- mats[[which.min(df$task)]]          # task 1 of this session
  pars <- map(mats, ~ rigid_params(.x %*% solve(ref)))
  df %>%
    mutate(
      x     = map_dbl(pars, 1), y    = map_dbl(pars, 2), z   = map_dbl(pars, 3),
      pitch = map_dbl(pars, 4), roll = map_dbl(pars, 5), yaw = map_dbl(pars, 6)
    )
}

motion_between <- runs %>%
  group_by(subject_id, session) %>%
  group_modify(~ between_one_session(.x)) %>%
  ungroup() %>%
  mutate(
    # same currency as FD above: displacement from task 1, rotations in mm
    disp_from_task1 = abs(x) + abs(y) + abs(z) +
      ROT_RADIUS * (abs(pitch) + abs(roll) + abs(yaw)),
    trans_from_task1 = sqrt(x^2 + y^2 + z^2),
    rot_from_task1_deg = sqrt(pitch^2 + roll^2 + yaw^2) * 180 / pi
  ) %>%
  group_by(subject_id, session) %>%
  # and the step from the previous task, which is the between-task equivalent
  # of a singleframewise displacement
  mutate(
    step_from_prev = abs(x - lag(x)) + abs(y - lag(y)) + abs(z - lag(z)) +
      ROT_RADIUS * (abs(pitch - lag(pitch)) + abs(roll - lag(roll)) +
                    abs(yaw - lag(yaw)))
  ) %>%
  ungroup() %>%
  select(subject_id, session, task, x, y, z, pitch, roll, yaw,
         trans_from_task1, rot_from_task1_deg, disp_from_task1, step_from_prev)

print(motion_between, n = Inf)


# ============================================================================
#  4. Plots  —  one line per subject, one panel per session
# ============================================================================
#  Same shape as the contrast spaghetti in look_at_data.R, so the two can be
#  read side by side: if a subject's odd run is motion, its line should spike
#  in the same place here.

p_within <- ggplot(motion_runs, aes(x = task, y = mean_fd,
                                    colour = subject_id, group = subject_id)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ session, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(motion_runs$task))) +
  labs(title = "Within-task motion",
       subtitle = "mean framewise displacement, rotations as mm at 50 mm radius",
       x = "task", y = "mean FD (mm)", colour = "subject") +
  theme_minimal()
print(p_within)

# the same plot for the worst single moment in each run — this is the one
# that catches a jolt the mean hides. Adding an aes() to a finished plot
# replaces just that one mapping; the data is already motion_runs.
p_within_max <- p_within +
  aes(y = max_fd) +
  geom_hline(yintercept = FD_THRESH, linetype = "dashed", colour = "grey60") +
  labs(title = "Worst single frame per task",
       subtitle = paste0("max framewise displacement; dashed line = FD_THRESH (",
                         FD_THRESH, " mm)"),
       y = "max FD (mm)")
print(p_within_max)

p_between <- ggplot(motion_between, aes(x = task, y = disp_from_task1,
                                        colour = subject_id,
                                        group = subject_id)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ session, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(motion_between$task))) +
  labs(title = "Between-task motion",
       subtitle = "displacement from the start of task 1, same session",
       x = "task", y = "displacement from task 1 (mm)", colour = "subject") +
  theme_minimal()
print(p_between)

# the FD timecourse itself, when a run looks odd and you want to see where
# in the run it went wrong
# p_time <- ggplot(motion_vols, aes(x = volume, y = fd, colour = factor(task))) +
#   geom_line(linewidth = 0.3) +
#   facet_grid(subject_id ~ session, labeller = label_both) +
#   labs(x = "volume", y = "FD (mm)", colour = "task") +
#   theme_minimal()
# print(p_time)


# ============================================================================
#  5. Save
# ============================================================================

if (SAVE_CSV) {
  write_csv(motion_runs,    file.path(OUTPUT_DIR, "motion_runs.csv"))
  write_csv(motion_between, file.path(OUTPUT_DIR, "motion_between.csv"))
  cat("Wrote motion_runs.csv and motion_between.csv to ", OUTPUT_DIR, "\n",
      sep = "")
}

# To line the motion up against the grid metrics, import_gridcat.R first,
# then join on the keys both tables share:
# grid_runs %>%
#   mutate(task = run_num) %>%
#   inner_join(motion_runs, by = c("subject_id", "session", "task")) %>%
#   ggplot(aes(x = mean_fd, y = contrast)) + geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x, se = FALSE)
