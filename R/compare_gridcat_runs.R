# ============================================================================
#  compare_gridcat_runs.R — compare two or more GridCAT runs
# ----------------------------------------------------------------------------
#  Import several runs first: set GRIDCAT_DIRS in import_gridcat.R to all the
#  output directories you want to compare, then source it, then this.
#
#  Three views:
#    1. which settings differ between the runs (from run_info)
#    2. the same measurement under each run, paired per subject/session/run
#    3. does phi (the grid orientation) agree between runs?
#
#  Plots are wrapped in print() because source() drops bare ggplots.
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

HEMI   <- "bilat"      # never pool hemispheres: each has its own orientation
METRIC <- "contrast"   # any numeric column of grid_runs: contrast, rayleigh_z, ...

# 6-fold symmetry: orientations live on a 60 deg circle, so the largest
# possible difference is 30 deg and chance is 15 deg on average.
PERIOD   <- 60
ori_diff <- function(a, b) abs(((a - b + PERIOD / 2) %% PERIOD) - PERIOD / 2)

stopifnot(exists("grid_runs"), exists("run_info"))
runs_here <- unique(grid_runs$gridcat_run)
if (length(runs_here) < 2) {
  stop("Only one GridCAT run imported (", runs_here,
       "). Add more directories to GRIDCAT_DIRS in import_gridcat.R.")
}


# ---- 1. which settings differ -------------------------------------------------

differing <- run_info %>%
  select(-any_of("created")) %>%
  summarise(across(everything(), ~ n_distinct(.x, na.rm = FALSE))) %>%
  pivot_longer(everything()) %>%
  filter(value > 1 | name == "gridcat_run") %>%
  pull(name)

cat("\n== settings that differ between the runs ==================\n")
run_info %>% select(all_of(differing)) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n== settings shared by all runs ============================\n")
run_info %>%
  select(-all_of(differing), -any_of("created")) %>%
  distinct() %>%
  pivot_longer(everything()) %>%
  as.data.frame() %>% print(row.names = FALSE)


# ---- 2. the same measurement under each run -----------------------------------

d <- grid_runs %>%
  filter(!is_average, hemisphere == HEMI) %>%
  mutate(unit = paste(subject_id, session, run, sep = "/"),
         value = .data[[METRIC]])

cat("\n== ", METRIC, " per run (", HEMI, ") =====================\n", sep = "")
d %>%
  group_by(gridcat_run) %>%
  summarise(n = sum(!is.na(value)),
            mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# Paired lines: each unit (subject/session/run) under every GridCAT run.
# Only units present in ALL runs are drawn, so the comparison stays paired.
complete_units <- d %>%
  filter(!is.na(value)) %>%
  count(unit) %>%
  filter(n == length(runs_here)) %>%
  pull(unit)

p_paired <- d %>%
  filter(unit %in% complete_units) %>%
  ggplot(aes(x = gridcat_run, y = value, group = unit)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(alpha = 0.3) +
  geom_point(aes(colour = subject_id), size = 1.5, alpha = 0.7) +
  facet_wrap(~ session, labeller = label_both) +
  labs(title = paste0(METRIC, " under each GridCAT run — ", HEMI),
       subtitle = "one line per subject/session/run; only units present in every run",
       x = NULL, y = METRIC, colour = "subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_paired)

# Pairwise scatter against the first run.
ref <- runs_here[1]
wide <- d %>%
  filter(unit %in% complete_units) %>%
  select(gridcat_run, unit, subject_id, session, value) %>%
  pivot_wider(names_from = gridcat_run, values_from = value)

for (other in setdiff(runs_here, ref)) {
  r <- cor(wide[[ref]], wide[[other]], use = "complete.obs")
  p_sc <- ggplot(wide, aes(x = .data[[ref]], y = .data[[other]])) +
    geom_abline(linetype = "dashed", colour = "grey60") +
    geom_point(aes(colour = subject_id), size = 2, alpha = 0.8) +
    labs(title = paste0(METRIC, ": ", other, " vs ", ref),
         subtitle = sprintf("r = %.2f over %d paired units", r, nrow(wide)),
         colour = "subject") +
    theme_minimal()
  print(p_sc)
}


# ---- 3. does phi agree between runs? -------------------------------------------
#  Chance is 15 deg; agreement well below that means both runs recover the
#  same orientation, not just a similar contrast.

phi <- grid_runs %>%
  filter(!is_average, hemisphere == HEMI, !is.na(ori_weighted)) %>%
  mutate(unit = paste(subject_id, session, run, sep = "/")) %>%
  select(gridcat_run, unit, phi = ori_weighted) %>%
  pivot_wider(names_from = gridcat_run, values_from = phi)

cat("\n== |delta phi| between runs (chance = 15 deg) ==============\n")
for (other in setdiff(runs_here, ref)) {
  dphi <- ori_diff(phi[[ref]], phi[[other]])
  ok <- !is.na(dphi)
  if (!any(ok)) next
  tt <- t.test(dphi[ok], mu = PERIOD / 4)
  cat(sprintf("%-30s vs %-30s  n=%3d  mean |dphi| = %5.2f deg  p vs chance = %.3f\n",
              other, ref, sum(ok), mean(dphi[ok]), tt$p.value))
}
cat("\n")
