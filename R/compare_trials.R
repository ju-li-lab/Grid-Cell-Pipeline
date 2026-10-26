# Compare preprocessing trials from grid_consistency.csv
#
# Each run of tools/grid_consistency/ appends rows tagged with its --trial
# label. This script asks whether one trial produced more reproducible grid
# orientations than another, while keeping the voxel-count confound visible
# rather than hidden.
#
# Edit the settings block, then Source.
#
# The plots are wrapped in print() on purpose: RStudio's Source button uses
# source(), which throws away a bare ggplot without drawing it.

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)


# ---- settings --------------------------------------------------------------

CSV <- "/Volumes/train_juli/SOLIS/Grid-Cell-Pipeline/grid_consistency.csv"

# Give the trials readable names. Run the script once with this empty to see
# the labels that are actually in the file, then fill it in.
#   names  = the trial label in the CSV
#   values = what to call it in the output
TRIAL_NAMES <- c(
  "20261026_114007" = "blip +1, thr 0.6",
  "20261026_123509" = "blip +1, thr 0.8",
  "20261026_131221" = "blip -1"
)

# Trial to compare the others against.
REFERENCE <- "blip +1, thr 0.8"

# Units with very few voxels make every consistency statistic unstable, and
# they are not comparable between trials that kept different numbers of
# voxels. Units below this in ANY trial are dropped from the matched
# comparison (they are still shown in the unmatched one).
MIN_VOXELS <- 10

HEMI <- NULL          # NULL = all hemispheres, or "bilat" / "left" / "right"


# ---- read ------------------------------------------------------------------

cons <- read_csv(CSV, show_col_types = FALSE)

if (length(TRIAL_NAMES) == 0 || !any(cons$trial %in% names(TRIAL_NAMES))) {
  cat("\nTrial labels present in", CSV, ":\n")
  print(cons %>% count(trial))
  stop("Fill in TRIAL_NAMES above with these labels, then re-run.")
}

cons <- cons %>%
  filter(trial %in% names(TRIAL_NAMES)) %>%
  mutate(trial = factor(unname(TRIAL_NAMES[trial]), levels = unname(TRIAL_NAMES)),
         unit  = paste(subject, session, hemi, comparison, sep = "/"))
if (!is.null(HEMI)) cons <- cons %>% filter(hemi == HEMI)

CHANCE <- unique(cons$chance_pct)[1]


# ---- which statistic to compare on -----------------------------------------
# Not z_diff, and not a resultant length. Both depend on how many voxels
# survived, and these trials differ a lot in that:
#
#   z_diff = n * r^2          grows with n for the SAME true consistency
#   rho    = sqrt((z-1)/(n-1))  corrects that, but explodes when n is tiny,
#                               so a trial that lost voxels looks better
#
# pct_stable has the same expectation at every n (CHANCE under the null), so
# it is the one that survives a change in voxel count. Pool it across voxels
# rather than averaging ROIs of wildly different size: a 3-voxel ROI should
# not carry the same weight as a 29-voxel one.

pooled <- function(d) {
  # Units with no usable voxels carry no information and would turn the
  # weighted mean into NA. Drop them, but count them - a trial that produced
  # empty ROIs is telling you something.
  d %>% summarise(
    units   = n(),
    empty   = sum(n_valid == 0 | is.na(pct_stable)),
    voxels  = sum(n_valid),
    pct     = 100 * sum(pct_stable / 100 * n_valid, na.rm = TRUE) /
                    sum(n_valid[!is.na(pct_stable)]),
    se      = 100 * sqrt(0.25 / sum(n_valid[!is.na(pct_stable)])),
    .groups = "drop")
}


# ---- 1. all units ----------------------------------------------------------

cat("\n=== all units (voxel counts differ between trials - see n_valid) ===\n")
cons %>% group_by(trial) %>% pooled() %>%
  mutate(vs_chance = sprintf("%+.1f", pct - CHANCE)) %>% print()

cat("\n=== how many voxels each trial had to work with ===\n")
cons %>% group_by(trial) %>%
  summarise(mask = sum(n_mask), valid = sum(n_valid),
            retained = sprintf("%.0f%%", 100 * sum(n_valid) / sum(n_mask)),
            `units n=0` = sum(n_valid == 0),
            `units n<5` = sum(n_valid < 5), `units n<10` = sum(n_valid < 10),
            .groups = "drop") %>% print()


# ---- 2. matched units ------------------------------------------------------
# The honest comparison available from this file: units that carry enough
# voxels in EVERY trial. A proper fixed voxel set needs the WriteValidMask /
# CommonMaskDir workflow in tools/grid_consistency/README.md.

ok_units <- cons %>% group_by(unit) %>%
  filter(n() == nlevels(trial), all(n_valid >= MIN_VOXELS)) %>%
  pull(unit) %>% unique()

matched <- cons %>% filter(unit %in% ok_units)

cat(sprintf("\n=== matched units only (n_valid >= %d in every trial: %d of %d) ===\n",
            MIN_VOXELS, length(ok_units), n_distinct(cons$unit)))
matched %>% group_by(trial) %>% pooled() %>%
  mutate(vs_chance = sprintf("%+.1f", pct - CHANCE)) %>% print()


# ---- 3. paired differences against the reference ---------------------------

ref <- matched %>% filter(trial == REFERENCE) %>%
  select(unit, ref_pct = pct_stable, ref_n = n_valid)

paired <- matched %>% filter(trial != REFERENCE) %>%
  inner_join(ref, by = "unit") %>%
  mutate(d_pct = pct_stable - ref_pct, d_n = n_valid - ref_n)

cat(sprintf("\n=== paired differences vs '%s' ===\n", REFERENCE))
paired %>% group_by(trial) %>%
  summarise(units = n(),
            median_d_pct = median(d_pct),
            better = sum(d_pct > 0), worse = sum(d_pct < 0),
            median_d_voxels = median(d_n),
            p_wilcox = tryCatch(wilcox.test(d_pct)$p.value, error = function(e) NA_real_),
            .groups = "drop") %>% print()


# ---- 4. plots --------------------------------------------------------------

p1 <- cons %>% group_by(trial) %>% pooled() %>%
  bind_rows(matched %>% group_by(trial) %>% pooled() %>% mutate(set = "matched")) %>%
  mutate(set = ifelse(is.na(set), "all units", set)) %>%
  ggplot(aes(x = trial, y = pct, fill = set)) +
  geom_hline(yintercept = CHANCE, linetype = "dashed", colour = "grey40") +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = pct - 1.96*se, ymax = pct + 1.96*se),
                position = position_dodge(0.8), width = 0.2) +
  annotate("text", x = 0.6, y = CHANCE, label = "chance", vjust = -0.5,
           hjust = 0, size = 3, colour = "grey40") +
  labs(title = "Cross-run grid orientation stability",
       subtitle = "pooled over voxels; error bars are 95% CI",
       x = NULL, y = "% stable voxels", fill = NULL) +
  theme_minimal()
print(p1)

# Does an apparent gain just track losing voxels? If the points slope, the
# trial that kept fewer voxels is being flattered by selection, not by
# better preprocessing.
p2 <- ggplot(paired, aes(x = d_n, y = d_pct, colour = trial)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
  labs(title = paste("Change in stability vs change in voxel count, vs", REFERENCE),
       subtitle = "a slope means the 'gain' is selection, not preprocessing",
       x = "change in valid voxels", y = "change in % stable", colour = NULL) +
  theme_minimal()
print(p2)

p3 <- ggplot(matched, aes(x = trial, y = pct_stable, group = unit)) +
  geom_hline(yintercept = CHANCE, linetype = "dashed", colour = "grey40") +
  geom_line(alpha = 0.3) +
  geom_point(aes(size = n_valid), alpha = 0.6) +
  scale_size_continuous(range = c(1, 4)) +
  labs(title = "Per unit, matched set", x = NULL, y = "% stable voxels",
       size = "voxels") +
  theme_minimal()
print(p3)


# ---- 5. what to conclude ---------------------------------------------------

cat("\n---- reading this ----\n")
cat("A trial is better if it is ABOVE chance by more, on the matched set,\n")
cat("without having thrown away voxels to get there. Check three things:\n")
cat("  1. the matched-set table, not the all-units one\n")
cat("  2. the retained-voxel column: a trial that kept far fewer voxels is\n")
cat("     suspect even if its percentage looks good\n")
cat("  3. p2: if the points slope, the difference is selection\n")
cat("If no trial clears chance by more than its CI, the honest answer is\n")
cat("that none of them has reproducible grid signal, and the setting you are\n")
cat("varying is not the thing holding you back.\n")
