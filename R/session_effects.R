# ============================================================================
#  Session differences, resolved per run
# ----------------------------------------------------------------------------
#  Built to answer one question properly: does the grid code differ BETWEEN
#  sessions, and is that difference the same in every run?
#
#  The thing this script exists to avoid: concluding that sessions differ
#  because one is significant and another is not. Those are two separate
#  tests against zero, and the gap between "p = 0.04" and "p = 0.09" can sit
#  on top of almost identical effects. A difference needs its own test, and
#  because the completeness filter guarantees the same subjects in every
#  session, that test can be paired - which is the most sensitive version
#  available to you.
#
#  Needs grid_runs and grid_stability from import_gridcat.R - run that first.
#
#  The plots are wrapped in print() on purpose: RStudio's Source button uses
#  source(), which throws away a bare ggplot without drawing it.
# ============================================================================

library(dplyr)
library(ggplot2)

# complete_subjects() lives in its own file so that editing or restoring
# import_gridcat.R cannot take it with it. Sourced automatically.
if (!exists("complete_subjects")) {
  .cs <- c("complete_subjects.R", "R/complete_subjects.R",
           "../R/complete_subjects.R")
  .cs <- .cs[file.exists(.cs)]
  if (!length(.cs)) stop(
    "Cannot find complete_subjects.R. In RStudio use\n",
    "  Session > Set Working Directory > To Source File Location,\n",
    "then Source this script again.")
  source(.cs[1]); rm(.cs)
}


# ---- settings --------------------------------------------------------------

HEMI <- "bilat"        # "left", "right" or "bilat"

METRIC <- "contrast"  # per-run metric to test. "contrast" is the grid code
                      # magnitude; "rayleigh_z" is spatial stability.

BASELINE <- 1         # the session everything else is compared against,
                      # normally the pre-intervention one.

REQUIRE_COMPLETE <- TRUE   # same rule as group_summary.R. Leave this on for
                      # anything comparing sessions: an unbalanced panel
                      # makes a session difference out of who turned up.


# ---- data ------------------------------------------------------------------

if (!exists("grid_runs")) stop("Run import_gridcat.R first - grid_runs is missing.")
if (!METRIC %in% names(grid_runs)) stop("No column called '", METRIC, "'.")

runs <- grid_runs %>% filter(hemisphere == HEMI, !is_average)

stab <- NULL
if (exists("grid_stability") && !is.null(grid_stability) && nrow(grid_stability)) {
  stab <- grid_stability %>%
    filter(hemisphere == HEMI, run_a != "avg", run_b != "avg")
}

if (REQUIRE_COMPLETE) {
  keep <- complete_subjects(runs, stab, value_col = METRIC)
  runs <- runs %>% filter(subject_id %in% keep)
  if (!is.null(stab)) stab <- stab %>% filter(subject_id %in% keep)
}

d <- runs %>%
  mutate(ses = factor(session), v = .data[[METRIC]]) %>%
  select(subject_id, ses, session, run, run_num, v,
         any_of(c("group", "n_voxels", "n_nan_voxels")))

sessions <- sort(unique(d$session))
others   <- setdiff(sessions, BASELINE)
if (!BASELINE %in% sessions) stop("Session ", BASELINE, " is not in the data.")
if (!length(others))         stop("Only one session present - nothing to compare.")

cat(sprintf("\n%s, %s hemisphere: %d subjects, sessions %s, runs %s\n",
            METRIC, HEMI, dplyr::n_distinct(d$subject_id),
            paste(sessions, collapse = "/"),
            paste(sort(unique(d$run_num)), collapse = "/")))


# ---- 1. the headline test --------------------------------------------------
#  One value per subject per session (averaged over runs), then a paired test
#  of each session against the baseline. This is the comparison your claim
#  about the intervention actually rests on.

subj <- d %>%
  group_by(subject_id, session) %>%
  summarise(v = mean(v, na.rm = TRUE), .groups = "drop")

paired_test <- function(a, b) {
  tt <- tryCatch(t.test(a, b, paired = TRUE), error = function(e) NULL)
  if (is.null(tt)) return(data.frame(n = NA, diff = NA, ci_lo = NA, ci_hi = NA,
                                     t = NA, p = NA, dz = NA))
  dif <- a - b
  data.frame(n = sum(!is.na(dif)), diff = unname(tt$estimate),
             ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
             t = unname(tt$statistic), p = tt$p.value,
             dz = mean(dif, na.rm = TRUE) / sd(dif, na.rm = TRUE))
}

wide <- subj %>% tidyr::pivot_wider(names_from = session, values_from = v,
                                    names_prefix = "s")

head_tbl <- lapply(others, function(s) {
  cbind(comparison = sprintf("session %d - session %d", BASELINE, s),
        paired_test(wide[[paste0("s", BASELINE)]], wide[[paste0("s", s)]]))
}) %>% bind_rows()

cat("\n=== session difference, averaged over runs (paired) ===\n")
head_tbl %>%
  mutate(across(c(diff, ci_lo, ci_hi, dz), ~round(.x, 4)),
         t = round(t, 2), p = signif(p, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nper-session means (same subjects throughout):\n")
subj %>% group_by(session) %>%
  summarise(n = n(), mean = round(mean(v, na.rm = TRUE), 4),
            sem = round(sd(v, na.rm = TRUE) / sqrt(n()), 4), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)


# ---- 2. the same difference, run by run ------------------------------------
#  Is the session effect spread evenly across runs, or carried by one of
#  them? A difference that lives in run 1 alone is more likely an order or
#  novelty effect than an intervention effect.
#
#  These are 2 sessions x 3 runs = 6 uncorrected tests. Treat them as a
#  description of where the overall effect sits, not as six findings.

wide_run <- d %>%
  select(subject_id, session, run_num, v) %>%
  tidyr::pivot_wider(names_from = session, values_from = v, names_prefix = "s")

run_tbl <- lapply(others, function(s) {
  lapply(sort(unique(wide_run$run_num)), function(r) {
    w <- wide_run %>% filter(run_num == r)
    cbind(run = r, comparison = sprintf("%d - %d", BASELINE, s),
          paired_test(w[[paste0("s", BASELINE)]], w[[paste0("s", s)]]))
  }) %>% bind_rows()
}) %>% bind_rows()

cat("\n=== session difference per run (paired, uncorrected) ===\n")
run_tbl %>%
  mutate(across(c(diff, ci_lo, ci_hi, dz), ~round(.x, 4)),
         t = round(t, 2), p = signif(p, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)


# ---- 3. does the session effect depend on run? -----------------------------
#  The omnibus version of panel 2. A session x run interaction says the
#  session difference is not the same in every run; no interaction says one
#  number describes it.

if (requireNamespace("lme4", quietly = TRUE)) {
  dm <- d %>% mutate(run_f = factor(run_num)) %>% filter(!is.na(v))
  fit <- function(f) suppressMessages(suppressWarnings(
    lme4::lmer(f, data = dm, REML = FALSE)))
  ok <- tryCatch({
    m_full  <- fit(v ~ ses * run_f + (1 | subject_id))
    m_add   <- fit(v ~ ses + run_f + (1 | subject_id))
    m_norun <- fit(v ~ ses + (1 | subject_id))
    m_noses <- fit(v ~ run_f + (1 | subject_id))
    cat("\n=== mixed model: v ~ session * run + (1 | subject) ===\n")
    cat("(likelihood-ratio tests; random intercept per subject)\n")
    lr <- function(a, b, what) {
      an <- anova(a, b)
      cat(sprintf("  %-22s chisq = %6.2f, df = %d, p = %s\n", what,
                  an$Chisq[2], an$Df[2], signif(an$`Pr(>Chisq)`[2], 3)))
    }
    lr(m_add,   m_full, "session x run")
    lr(m_noses, m_add,  "session (main)")
    lr(m_norun, m_add,  "run (main)")
    TRUE
  }, error = function(e) { cat("\n(mixed model skipped: ", conditionMessage(e), ")\n", sep=""); FALSE })
} else {
  cat("\n(lme4 not installed - skipping the interaction test)\n")
}


# ---- 4. is this an intervention effect, or just time? ----------------------
#  A session difference that shows up in EVERY group is not evidence for the
#  intervention: repetition, familiarity with the task, and more direct
#  routes through the environment all push in the same direction as a
#  session number. Only a group x session difference isolates the
#  intervention, and that needs a group that did not receive it.

has_groups <- "group" %in% names(d) &&
  dplyr::n_distinct(d$group[!is.na(d$group)]) >= 2

if (has_groups) {
  cat("\n=== session difference within each group ===\n")
  d %>%
    filter(!is.na(group)) %>%
    group_by(group, subject_id, session) %>%
    summarise(v = mean(v, na.rm = TRUE), .groups = "drop") %>%
    group_by(group, session) %>%
    summarise(n = n(), mean = round(mean(v, na.rm = TRUE), 4),
              sem = round(sd(v, na.rm = TRUE) / sqrt(n()), 4), .groups = "drop") %>%
    as.data.frame() %>% print(row.names = FALSE)
  cat("\nThe intervention claim lives in the DIFFERENCE between these two\n")
  cat("rows of differences, not in either one on its own.\n")
} else {
  cat("\n(no group column filled in - see note 4 in the closing text)\n")
}


# ---- 5. plots --------------------------------------------------------------

sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# P1: each session's run profile.
s1 <- d %>% group_by(ses, run_num) %>%
  summarise(mean = mean(v, na.rm = TRUE), se = sem(v), .groups = "drop")

p1 <- ggplot(s1, aes(x = run_num, y = mean, colour = ses, group = ses)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.08) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = sort(unique(s1$run_num))) +
  labs(title = paste0(METRIC, " by run and session - ", HEMI, " ErC"),
       subtitle = "group mean +- SEM, same subjects in every session",
       x = "run", y = METRIC, colour = "session") +
  theme_minimal()
print(p1)

# P2: the paired difference itself, per run. This is the plot that carries
# the claim: if the ribbon crosses zero, that run shows no session effect.
p2 <- ggplot(run_tbl, aes(x = run, y = diff, colour = comparison,
                          group = comparison)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.08,
                position = position_dodge(0.15)) +
  geom_line(position = position_dodge(0.15)) +
  geom_point(size = 2.5, position = position_dodge(0.15)) +
  scale_x_continuous(breaks = sort(unique(run_tbl$run))) +
  labs(title = "Within-subject session difference, per run",
       subtitle = "95% CI of the paired difference; crossing zero = no effect in that run",
       x = "run", y = paste("difference in", METRIC), colour = NULL) +
  theme_minimal()
print(p2)

# P3: every subject's own trajectory. At this sample size a group difference
# can be two people, and this is where you see that.
p3 <- ggplot(d, aes(x = ses, y = v, group = subject_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(alpha = 0.35) +
  geom_point(size = 1.4, alpha = 0.6) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               colour = "firebrick", linewidth = 1.1) +
  stat_summary(aes(group = 1), fun = mean, geom = "point",
               colour = "firebrick", size = 2.5) +
  facet_wrap(~ run_num, labeller = label_both) +
  labs(title = "Per subject across sessions, one panel per run",
       subtitle = "red is the group mean",
       x = "session", y = METRIC) +
  theme_minimal()
print(p3)

# P4: the second metric. Stability is defined on PAIRS of runs, so it cannot
# go on a per-run axis - it gets its own panel, against its 50% chance line.
if (!is.null(stab) && nrow(stab)) {
  s4 <- stab %>%
    group_by(ses = factor(session), comparison) %>%
    summarise(mean = mean(pct_stable, na.rm = TRUE), se = sem(pct_stable),
              .groups = "drop")

  p4 <- ggplot(s4, aes(x = comparison, y = mean, colour = ses, group = ses)) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1,
                  position = position_dodge(0.2)) +
    geom_point(size = 2.5, position = position_dodge(0.2)) +
    labs(title = paste0("Temporal stability by run-pair and session - ", HEMI, " ErC"),
         subtitle = "dashed line is the 50% chance level",
         x = "run pair", y = "% stable voxels", colour = "session") +
    theme_minimal()
  print(p4)
}


# ---- how to read this ------------------------------------------------------

cat("
---- reading this ----

1. Table 1 is the test. 'Session 1 was significant and session 2 was not' is
   not a session difference - it is two separate tests against zero. What
   supports an intervention claim is a paired difference whose confidence
   interval excludes zero, and that is what table 1 reports.

2. Table 2 and plot P2 say WHERE that difference sits. An effect spread
   evenly over runs is the one that behaves like a session effect. One that
   lives entirely in run 1 is more likely novelty or order, and one that
   grows across runs within a session looks like fatigue. The mixed model's
   'session x run' line is the formal version of the same question.

3. Check P3 before believing any of it. At this sample size a group mean can
   be moved by two subjects, and the paired tests above will not tell you
   that - they will just report a small p.

4. The confound that matters most here: session number is also time, task
   familiarity, and whatever subjects learned about the environment. All of
   those push in the same direction as an intervention. A session difference
   found in a single group cannot separate them. If you have subjects who
   did not receive the intervention, fill in the group column of
   subject_info.csv - section 4 then gives you the comparison that can.

5. On the ROI change you are running now: a whole-ErC mask instead of the
   posterior subdivision is a good move on its own terms - it is what the
   benchmark paper used, and more voxels means a less noisy per-subject
   estimate, which is the main thing limiting you. But it is a different
   analysis, not a confirmation of this one. Decide which mask is the
   primary one before you see both results, and report the other as what it
   is. Running the same comparison two ways and keeping the better number is
   how a p-value stops meaning anything.
")
