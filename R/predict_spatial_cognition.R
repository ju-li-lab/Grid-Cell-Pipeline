# ============================================================================
#  predict_spatial_cognition.R — does a grid measure predict spatial
#  cognition? ONE GridCAT run at a time.
# ----------------------------------------------------------------------------
#  Needs grid_runs + grid_stability (import_gridcat.R) and spatial_cog
#  (import_spatial_cognition.R) — source those first. Groups are read fresh
#  from subject_info.csv.
#
#  For every grid measure (A) x hemisphere x spatial-cognition score (B):
#    overall  do subjects with a higher A have a higher B? Correlation of
#             each subject's mean over the MRI sessions.
#    change   does a drop in A go with a drop in B? Correlation of the
#             changes from the first MRI session, one test per later session.
#  r > 0: A and B rise and fall together. Every test stands on its own and
#  is not corrected; the summary says how many p < .05 the same tests give
#  with the subjects shuffled, i.e. from noise alone.
#
#  You get the tables `overall` and `change` (also saved to
#  results/<GridCAT run>/), a heatmap of every r, and the scatter plots
#  behind the tests with p < .05.
#  Plots are wrapped in print() because source() drops bare ggplots.
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(ggplot2)

GRIDCAT_RUN   <- NULL      # NULL = first imported run, or name one from run_info
PREDICTORS    <- c("contrast", "temporal_stability", "spatial_stability",
                   "ori_consistency", "usable_voxels")
MEASURES      <- "Score"   # NULL = every measure, or e.g. c("Score", "RT")
PER_CONDITION <- FALSE     # TRUE = every condition too, not only their mean
BY_GROUP      <- FALSE     # TRUE = each group of subject_info.csv on its own
MIN_SUBJECTS  <- 5         # fewer subjects -> no test
N_SHUFFLES    <- 200       # subject shuffles for the chance level; 0 = skip
PERIOD        <- 60        # 60 for 6-fold symmetry, 90 for 4-fold

SUBJECT_INFO_CSV <- "subject_info.csv"
SAVE_CSV <- TRUE

stopifnot(exists("grid_runs"), exists("grid_stability"), exists("spatial_cog"))
if (is.null(GRIDCAT_RUN)) GRIDCAT_RUN <- unique(grid_runs$gridcat_run)[1]


# ---- A: grid measures, one row per subject x session x hemisphere ------------
#  contrast            mean of the runs' contrasts (the "avg" row is a sum)
#  temporal_stability  mean % stable voxels over the run-vs-run pairs; the
#                      "avg_vs" pairs are circular and left out
#  spatial_stability   mean between-voxel Rayleigh z of the runs
#  ori_consistency     how well the runs' ROI orientations agree: resultant
#                      length on the PERIOD circle, 1 = identical runs,
#                      about 0.5 for three random orientations
#  usable_voxels       ROI voxels SPM's mask left. Not a grid measure: a
#                      yardstick for what data quality alone predicts

grid_measures <- grid_runs %>%
  filter(gridcat_run == GRIDCAT_RUN, !is_average) %>%
  group_by(subject_id, session, hemisphere) %>%
  summarise(
    contrast          = mean(contrast, na.rm = TRUE),
    spatial_stability = mean(rayleigh_z, na.rm = TRUE),
    ori_consistency   = if (sum(!is.na(ori_weighted)) > 1) {
      Mod(mean(exp(2i * pi * ori_weighted / PERIOD), na.rm = TRUE))
    } else NA_real_,
    usable_voxels     = max(n_voxels - n_nan_voxels),
    .groups = "drop") %>%
  left_join(
    grid_stability %>%
      filter(gridcat_run == GRIDCAT_RUN, run_a != "avg", run_b != "avg") %>%
      group_by(subject_id, session, hemisphere) %>%
      summarise(temporal_stability = mean(pct_stable, na.rm = TRUE),
                .groups = "drop"),
    by = c("subject_id", "session", "hemisphere"))
stopifnot(nrow(grid_measures) > 0)

if (file.exists(SUBJECT_INFO_CSV)) {
  grid_measures <- left_join(grid_measures,
                             read_csv(SUBJECT_INFO_CSV, show_col_types = FALSE),
                             by = "subject_id")
}
if (!"group" %in% names(grid_measures)) grid_measures$group <- NA_character_


# ---- B: spatial-cognition scores, one row per subject x MRI session x score ---
#  File sessions that stand for the same MRI session are averaged; a test
#  with several conditions is scored by their mean ("all").

sc_scores <- spatial_cog %>%
  filter(!is.na(session), subject_id %in% grid_measures$subject_id,
         is.null(MEASURES) | measure %in% MEASURES) %>%
  group_by(subject_id, session, test, condition, measure) %>%
  summarise(score = mean(value), .groups = "drop")
condition_means <- sc_scores %>%
  group_by(test) %>%
  filter(n_distinct(condition) > 1) %>%
  group_by(subject_id, session, test, measure) %>%
  summarise(condition = "all", score = mean(score), .groups = "drop")
sc_scores <- bind_rows(sc_scores, condition_means)
if (!PER_CONDITION) sc_scores <- filter(sc_scores, condition == "all")

# Named by test alone ("SUT") unless conditions or several measures are in play
several_measures <- n_distinct(sc_scores$measure) > 1
sc_scores <- sc_scores %>%
  arrange(match(test, unique(spatial_cog$test)), condition != "all", condition, measure) %>%
  mutate(score_name = str_squish(paste(test,
                                       if_else(condition == "all", "", condition),
                                       if (several_measures) measure else "")),
         score_name = factor(score_name, levels = unique(score_name)))


# ---- the tests --------------------------------------------------------------------

keys       <- c(if (BY_GROUP) "group", "hemisphere", "predictor", "score_name")
group_keys <- if (BY_GROUP) "group" else character(0)

pair_up <- function(gm) {
  gm %>%
    pivot_longer(all_of(PREDICTORS), names_to = "predictor", values_to = "A") %>%
    inner_join(select(sc_scores, subject_id, session, score_name, B = score),
               by = c("subject_id", "session"), relationship = "many-to-many") %>%
    filter(!is.na(A), !is.na(B))
}

# overall: each subject's mean over the sessions with both A and B
overall_values <- function(pairs) {
  pairs %>%
    group_by(across(all_of(c(keys, "subject_id")))) %>%
    summarise(A = mean(A), B = mean(B), .groups = "drop")
}

# change: every later session minus the first one
change_values <- function(pairs) {
  first <- min(pairs$session)
  pairs %>%
    group_by(across(all_of(c(keys, "subject_id")))) %>%
    filter(any(session == first)) %>%
    mutate(A = A - A[session == first], B = B - B[session == first]) %>%
    ungroup() %>%
    filter(session != first) %>%
    mutate(change = paste(session, "-", first))
}

# Pearson r across subjects; same p as lm(B ~ A)
correlate <- function(values, by) {
  values %>%
    group_by(across(all_of(by))) %>%
    filter(n() >= MIN_SUBJECTS) %>%
    summarise(r = suppressWarnings(cor(A, B)), n = n(), .groups = "drop") %>%
    filter(!is.na(r)) %>%
    mutate(p = 2 * pt(-abs(r) * sqrt((n - 2) / (1 - r^2)), n - 2))
}

pairs   <- pair_up(grid_measures)
overall <- correlate(overall_values(pairs), keys)
change  <- correlate(change_values(pairs), c(keys, "change"))


# ---- how many p < .05 would noise give? -------------------------------------------
#  Shuffle which subject's grid measures go with which subject's scores and
#  count again. The shuffles keep every overlap between the tests (left and
#  right inside bilat, the scores of one person), so their spread is the
#  honest chance level.

chance <- NULL
if (N_SHUFFLES > 0) {
  set.seed(1)
  chance <- map(seq_len(N_SHUFFLES), function(i) {
    shuffled <- grid_measures %>%
      group_by(across(all_of(group_keys))) %>%
      mutate(subject_id = {
        ids <- unique(subject_id)
        setNames(ids[sample.int(length(ids))], ids)[subject_id]
      }) %>%
      ungroup() %>%
      pair_up()
    tibble(overall = sum(correlate(overall_values(shuffled), keys)$p < 0.05),
           change  = sum(correlate(change_values(shuffled), c(keys, "change"))$p < 0.05))
  }) %>%
    list_rbind()
}

summary_line <- function(label, tests, shuffled_hits) {
  line <- sprintf("%-8s %3d of %3d tests with p < .05", label,
                  sum(tests$p < 0.05), nrow(tests))
  if (length(shuffled_hits)) {
    line <- paste0(line, sprintf("   shuffled subjects: %.1f on average, 95%% at or below %d",
                                 mean(shuffled_hits),
                                 as.integer(quantile(shuffled_hits, 0.95, type = 1))))
  }
  cat(line, "\n")
}

cat("\n== ", GRIDCAT_RUN, ": does a grid measure predict spatial cognition? ==\n", sep = "")
cat("scores: ", paste(levels(sc_scores$score_name), collapse = ", "), "\n", sep = "")
cat("subjects: ", n_distinct(pairs$subject_id), " overall, ",
    n_distinct(change_values(pairs)$subject_id), " with two sessions\n\n", sep = "")
summary_line("overall", overall, chance$overall)
summary_line("change", change, chance$change)

results <- bind_rows(
  mutate(overall, question = "overall"),
  change %>% mutate(question = paste("change", change)) %>% select(-change)
) %>%
  mutate(question = factor(question, levels = unique(question)))

cat("\ntests with p < .05 (uncorrected):\n")
if (any(results$p < 0.05)) {
  results %>%
    filter(p < 0.05) %>%
    arrange(p) %>%
    select(question, all_of(keys), r, p, n) %>%
    mutate(r = round(r, 2), p = signif(p, 2)) %>%
    as.data.frame() %>% print(row.names = FALSE)
} else {
  cat("  none\n")
}
cat("\n")

if (isTRUE(SAVE_CSV)) {
  out_dir <- file.path("results", GRIDCAT_RUN)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(overall, file.path(out_dir, "overall.csv"))
  write_csv(change,  file.path(out_dir, "change.csv"))
  message("Saved overall.csv and change.csv to ", normalizePath(out_dir))
}


# ---- every r at a glance ------------------------------------------------------------
#  One heatmap; with more than 20 scores, one per test.

heatmap <- function(res, title) {
  ggplot(res, aes(x = predictor, y = score_name, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(sprintf("%.2f", r), if_else(p < 0.05, "*", ""))),
              size = 3, colour = "#0b0b0b") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948",
                         midpoint = 0, limits = c(-1, 1)) +
    scale_x_discrete(limits = PREDICTORS) +
    scale_y_discrete(limits = rev) +
    facet_grid(rows = if (BY_GROUP) vars(question, group) else vars(question),
               cols = vars(hemisphere), scales = "free_y", space = "free_y") +
    labs(title = title, subtitle = "r across subjects; * p < .05, uncorrected",
         x = NULL, y = NULL, fill = "r") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

title <- paste0("Does a grid measure predict spatial cognition? (", GRIDCAT_RUN, ")")
if (nlevels(droplevels(results$score_name)) <= 20) {
  print(heatmap(results, title))
} else {
  results$test <- word(as.character(results$score_name), 1)
  for (t in unique(results$test)) {
    print(heatmap(filter(results, test == t), paste0(title, " — ", t)))
  }
  results$test <- NULL
}


# ---- the data behind the tests with p < .05 ----------------------------------------
#  overall: one point per subject (its mean over the sessions);
#  change: one point per subject (later minus first session).

hits <- results %>% filter(p < 0.05) %>% arrange(p) %>% head(12)
if (nrow(hits)) {
  values <- bind_rows(
    mutate(overall_values(pairs), question = "overall"),
    change_values(pairs) %>% mutate(question = paste("change", change))
  )
  p_hits <- hits %>%
    mutate(panel = sprintf("%s: %s %s vs %s\nr = %.2f, p = %.2g, n = %d",
                           question, hemisphere, predictor, score_name, r, p, n)) %>%
    select(all_of(keys), question, panel) %>%
    inner_join(values, by = c(keys, "question")) %>%
    ggplot(aes(x = A, y = B)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#2a78d6", linewidth = 0.7) +
    geom_point(colour = "#2a78d6", size = 2.5) +
    geom_text(aes(label = str_remove(subject_id, "^sub-")),
              size = 2.5, vjust = -0.9, colour = "#52514e") +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    facet_wrap(~ panel, scales = "free") +
    labs(title = paste0("Tests with p < .05 (", GRIDCAT_RUN, ")"),
         x = "grid measure (A)", y = "spatial-cognition score (B)") +
    theme_minimal() +
    theme(strip.text = element_text(size = 8))
  print(p_hits)
}

rm(condition_means, several_measures, pairs, chance, hits, title)


# ---- what you get --------------------------------------------------------------------
#  grid_measures: one row per subject x session x hemisphere of this run, with
#    the grid measures above and subject_info.
#  sc_scores: one row per subject x MRI session x score.
#  overall: one row per (group x) hemisphere x grid measure x score: r, n, p.
#  change: the same per later session ("2 - 1" = session 2 minus session 1).
