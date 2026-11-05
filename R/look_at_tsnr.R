# Look at the tSNR table from tools/tsnr_report.sh, per subject.
# Spaghetti plots: one line per subject, one panel per session, so a run that
# is noisier than its neighbours shows up as a kink rather than disappearing
# into a group mean.
#
# This script imports as well as plots — the tSNR data is one CSV, so there is
# nothing for a separate importer to do. Column names are made to match
# import_gridcat.R (subject_id, session, hemisphere, run_num) so the two can
# be joined; the last section does that.
#
# The plots are wrapped in print() on purpose: RStudio's Source button uses
# source(), which throws away a bare ggplot without drawing it.

library(dplyr)
library(ggplot2)


# ---- settings --------------------------------------------------------------

TSNR_CSV <- "/Users/juli/Desktop/SOLIS/tsnr_test2.tsv"     # from tools/tsnr_report.sh
HEMI     <- "left"                 # "left" or "right"; see the print() below
METRIC   <- "roi_tsnr"             # any numeric column: roi_tsnr, brain_tsnr,
                                   # roi_mean, n_vox, roi_over_mthresh

# Found whether the working directory is the repo root or R/ itself, because
# RStudio's Session > Set Working Directory > To Source File Location puts you
# in R/ and a path written for one is wrong in the other.
find_csv <- function(name) {
  hit <- Filter(file.exists, c(name, file.path("R", name),
                               file.path("..", "R", name)))
  if (!length(hit)) stop("Cannot find ", name, ". Run tools/tsnr_report.sh ",
                         "first, or set TSNR_CSV to a full path.")
  hit[1]
}

# The tSNR CSV leaves the mask name unparsed, because ROI naming changes
# between studies. These pull the BIDS-style tags back out when they are
# there and give NA when they are not, so a differently-named mask still
# imports — it just lands in `roi` under its own filename.
tag <- function(x, key) {
  pat <- paste0(".*", key, "-([A-Za-z0-9]+).*")
  ifelse(grepl(pat, x), sub(pat, "\\1", x), NA_character_)
}
first_int <- function(x) {
  suppressWarnings(as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1", x)))
}


# ---- import ----------------------------------------------------------------

tsnr <- readr::read_csv(find_csv(TSNR_CSV), show_col_types = FALSE) %>%
  rename(subject_id = subject) %>%
  mutate(
    session    = first_int(session),            # "ses-01" -> 1
    run_num    = first_int(run),                # "3_run-2" -> 3, "rest" -> NA
    hemisphere = tag(mask, "hemi"),
    roi        = coalesce(tag(mask, "label"), mask)
  ) %>%
  arrange(subject_id, session, hemisphere, run_num)

# What is actually in the file. If HEMI below finds nothing, this is why.
cat("subjects:    ", paste(unique(tsnr$subject_id), collapse = ", "), "\n",
    "sessions:    ", paste(sort(unique(tsnr$session)), collapse = ", "), "\n",
    "hemispheres: ", paste(unique(tsnr$hemisphere), collapse = ", "), "\n",
    "ROIs:        ", paste(unique(tsnr$roi), collapse = ", "), "\n",
    "images:      ", paste(unique(tsnr$image), collapse = ", "), "\n", sep = "")


# ---- pick what to look at --------------------------------------------------

d <- tsnr %>% filter(hemisphere == HEMI, !is.na(run_num))

if (nrow(d) == 0) {
  stop("No rows for hemisphere '", HEMI, "'. Available: ",
       paste(unique(tsnr$hemisphere), collapse = ", "))
}

d %>%
  select(subject_id, session, run, n_vox, roi_tsnr, brain_tsnr,
         roi_mean, spm_global, roi_over_mthresh) %>%
  print(n = Inf)


# ---- the spaghetti plot ----------------------------------------------------
# run_num on the x axis rather than the `run` factor, because it is numeric:
# that keeps the points evenly spaced and lets the line join them in run
# order. The resting-state run is dropped above (run_num is NA for it) — it
# is a different acquisition, not a fourth task run.

p1 <- ggplot(d, aes(x = run_num, y = .data[[METRIC]],
                    colour = subject_id, group = subject_id)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ session, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(d$run_num))) +
  labs(title = paste0(METRIC, " across runs — ", HEMI, " ", unique(d$roi)[1]),
       x = "run", y = METRIC, colour = "subject") +
  theme_minimal()
print(p1)


# ---- the same plot for the one metric with a meaningful line ---------------
# roi_over_mthresh is roi_mean / (0.8 * spm_global): below 1 the ROI sits
# under SPM's implicit-mask threshold and the GLM discards those voxels.
# That makes 1.0 a real boundary rather than a reference, so it gets its own
# plot instead of riding on METRIC above.

p2 <- ggplot(d, aes(x = run_num, y = roi_over_mthresh,
                    colour = subject_id, group = subject_id)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ session, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(d$run_num))) +
  labs(title = paste0("ROI mean vs SPM's masking threshold — ", HEMI),
       subtitle = "below the dashed line the GLM drops the voxel",
       x = "run", y = "roi_mean / (0.8 x global)", colour = "subject") +
  theme_minimal()
print(p2)


# ---- ROI against whole brain ----------------------------------------------
# Does a low ROI tSNR mean a bad region, or a bad run? If the two move
# together the run was noisy overall; if only the ROI drops, it is the region.

# p3 <- ggplot(d, aes(x = brain_tsnr, y = roi_tsnr, colour = subject_id)) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
#   geom_point(size = 2) +
#   facet_wrap(~ session, labeller = label_both) +
#   labs(title = paste0("ROI vs whole-brain tSNR — ", HEMI),
#        subtitle = "dashed line is ROI = brain",
#        x = "whole-brain tSNR", y = "ROI tSNR", colour = "subject") +
#   theme_minimal()
# print(p3)

# both hemispheres at once, one row each
# tsnr %>% filter(!is.na(run_num)) %>%
#   ggplot(aes(run_num, roi_tsnr, colour = subject_id, group = subject_id)) +
#   geom_line(alpha = 0.6) + geom_point() +
#   facet_grid(hemisphere ~ session, labeller = label_both) + theme_minimal()

# one panel per subject, when the colour legend stops being readable
# p1 + facet_grid(session ~ subject_id) + theme(legend.position = "none")

# the spaghetti with the group mean drawn over it
# p1 + stat_summary(aes(group = 1), fun = mean, geom = "line",
#                   colour = "black", linewidth = 1)


# ---- does tSNR explain the contrast? ---------------------------------------
# Only meaningful if BOTH CSVs come from the same derivatives tree: the masks
# are re-coregistered per tree, so joining tsnr_viarest.csv to a grid_runs.csv
# exported from the oap tree compares different voxels. Check the paths first.

GRID_RUNS <- "grid_runs.csv"       # from import_gridcat.R
found    <- Filter(file.exists, c(GRID_RUNS, file.path("R", GRID_RUNS),
                                  file.path("..", "R", GRID_RUNS)))

if (length(found)) {
  gr <- readr::read_csv(found[1], show_col_types = FALSE) %>%
    filter(!is_average) %>%
    select(subject_id, session, hemisphere, run_num, contrast, rayleigh_z)

  both <- inner_join(d, gr, by = c("subject_id", "session",
                                   "hemisphere", "run_num"))
  cat("\njoined rows:", nrow(both), "of", nrow(d), "tSNR rows\n")

  p4 <- ggplot(both, aes(x = roi_tsnr, y = contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(aes(colour = subject_id), size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "black", linewidth = 0.5) +
    labs(title = paste0("Contrast against ROI tSNR — ", HEMI),
         subtitle = "one point per run; a slope here means tSNR limits the effect",
         x = "ROI tSNR", y = "aligned vs misaligned", colour = "subject") +
    theme_minimal()
  print(p4)
} else {
  cat("\n", GRID_RUNS, " not found — skipping the contrast join.\n", sep = "")
}

# the plots stay around as p1, p2, p4 — type p1 to draw it again, or save it
# ggsave("tsnr_left.png", p1, width = 8, height = 4, dpi = 300)
