# ============================================================================
#  paths.R  —  What route did each run actually travel?
# ----------------------------------------------------------------------------
#  The navigation is passive, so the trajectory is set by the experiment, not by
#  the subject. This reads the EventData files and answers three things:
#
#    1. how many distinct routes there are, and which run used which
#    2. whether each route samples headings evenly, overall and within the two
#       GLM1/GLM2 splits — if it does not, the recovered grid orientation is
#       partly a readout of where the route went
#    3. whether the grid orientation repeats when the SAME subject walks the
#       SAME route in a different session, which is the cleanest test this
#       design allows
#
#  Interactive version: https://claude.ai/artifact/3R6N4nqr4YLSaax9Un4VRH
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

EVENT_DIR <- "/Volumes/train_juli/SOLIS/GLM_runauto/EventData"
GRID_RUNS <- "R/grid_runs.csv" # from import_gridcat.R, for the section 3 test
HEMI      <- "left"

# 6-fold symmetry: orientations live on a 60 deg circle, so the largest possible
# difference is 30 deg and two random orientations differ by 15 deg on average.
PERIOD    <- 60
NULL_MEAN <- PERIOD / 4
ori_diff  <- function(a, b) abs(((a - b + PERIOD/2) %% PERIOD) - PERIOD/2)


# ---- read the event files --------------------------------------------------
#  Three event types; only `translate` carries a heading, and those are the
#  grid events. The macOS `._` resource forks are not data.

read_events <- function(path) {
  read.csv(path, header = FALSE, col.names = c("type", "onset", "dur", "dir"),
           fill = TRUE, stringsAsFactors = FALSE) %>%
    mutate(dir = suppressWarnings(as.numeric(dir)))
}

files <- list.files(EVENT_DIR, pattern = "_EventData\\.txt$", full.names = TRUE)
files <- files[!startsWith(basename(files), "._")]
stopifnot(length(files) > 0)

events <- map_dfr(files, function(p) {
  id <- str_match(basename(p), "^sub-(\\d+[a-z]\\d+)_ses-(\\d+)_task-(\\d+)_")
  read_events(p) %>%
    mutate(subject = id[2], session = as.integer(id[3]), run = as.integer(id[4]))
})

translates <- events %>% filter(type == "translate")


# ============================================================================
#  1. How many distinct routes?
# ----------------------------------------------------------------------------
#  A route is the ordered sequence of translate headings. Two runs share a route
#  if that whole sequence matches.
# ============================================================================

routes <- translates %>%
  group_by(subject, session, run) %>%
  summarise(seq = paste(dir, collapse = ","),
            n_legs = n(),
            path_seconds = sum(dur), .groups = "drop")

# Most-used route first, ties broken on the sequence itself so the names are
# stable between runs of this script. One sequence is another plus a single
# trailing event (the scanner stopping mid-leg); that is the same route, so it
# is folded in as "P3+" rather than counted as a fourth.
route_id <- routes %>%
  count(seq, name = "n_runs") %>%
  arrange(desc(n_runs), seq) %>%
  mutate(route = paste0("P", row_number()))

is_tail_of <- function(a, b) startsWith(a, paste0(b, ","))
for (i in seq_len(nrow(route_id))) {
  parent <- which(map_lgl(route_id$seq, ~ is_tail_of(route_id$seq[i], .x)))
  parent <- setdiff(parent, i)
  if (length(parent)) route_id$route[i] <- paste0(route_id$route[parent[1]], "+")
}

routes <- left_join(routes, route_id, by = "seq")

cat("\n== distinct routes ========================================\n")
route_id %>%
  mutate(n_legs = lengths(strsplit(seq, ","))) %>%
  select(route, n_runs, n_legs) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n", nrow(routes), " runs ran ",
    sum(!str_detect(route_id$route, "\\+")), " distinct routes",
    if (any(str_detect(route_id$route, "\\+")))
      " (a \"+\" route is the same route with a truncated tail)" else "",
    ".\n", sep = "")

cat("\n== which route each run used ==============================\n")
routes %>%
  select(subject, session, run, route) %>%
  pivot_wider(names_from = run, values_from = route, names_prefix = "run") %>%
  arrange(subject, session) %>%
  as.data.frame() %>% print(row.names = FALSE)

# How precisely is a route replayed? Compare leg durations across the runs that
# share it. Anything near zero means the trajectory is literally replayed.
cat("\n== replay precision ======================================\n")
translates %>%
  left_join(routes %>% select(subject, session, run, route),
            by = c("subject", "session", "run")) %>%
  group_by(route, subject, session, run) %>%
  mutate(leg = row_number()) %>%
  group_by(route, leg) %>%
  summarise(spread = max(dur) - min(dur), .groups = "drop_last") %>%
  summarise(worst_leg_spread_s = round(max(spread), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)


# Every run's LAST translate event is cut short by the scanner stopping. A
# zero-duration event is degenerate in an SPM design and a 0.06 s one carries
# no signal, yet both still enter the parametric modulator.
cat("\n== final translate event of each run ======================\n")
last_legs <- translates %>%
  group_by(subject, session, run) %>%
  slice_tail(n = 1) %>% ungroup()
cat("  runs ending in a leg under 0.5 s: ", sum(last_legs$dur < 0.5),
    " of ", nrow(last_legs), "\n", sep = "")
cat("  of those, exactly zero seconds:   ", sum(last_legs$dur == 0), "\n", sep = "")
cat("  shortest NON-final translate:     ",
    round(min(translates %>% group_by(subject, session, run) %>%
              slice_head(n = -1) %>% pull(dur)), 2), " s\n", sep = "")


# ============================================================================
#  2. Does each route sample headings evenly?
# ----------------------------------------------------------------------------
#  The modulators are cos(6*theta) and sin(6*theta). They only behave if the
#  headings are spread evenly over the 60 deg wedge. Two checks:
#    resultant length  — 0 is perfectly even; random sampling gives ~1/sqrt(n)
#    corr(cos, sin)    — should be ~0, or the two regressors trade off
# ============================================================================

resultant <- function(dir) {
  x <- mean(cos(6 * dir * pi / 180)); y <- mean(sin(6 * dir * pi / 180))
  sqrt(x^2 + y^2)
}
collinearity <- function(dir) {
  suppressWarnings(cor(cos(6 * dir * pi / 180), sin(6 * dir * pi / 180)))
}

one_run_per_route <- routes %>% group_by(route) %>% slice(1) %>% ungroup()

cat("\n== heading sampling, whole run ===========================\n")
one_run_per_route %>%
  left_join(translates, by = c("subject", "session", "run")) %>%
  group_by(route) %>%
  summarise(n          = n(),
            resultant  = round(resultant(dir), 4),
            chance     = round(1 / sqrt(n()), 4),
            cor_cos_sin= round(collinearity(dir), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n== heading sampling, within each split's half of the events =\n")
one_run_per_route %>%
  left_join(translates, by = c("subject", "session", "run")) %>%
  group_by(route) %>%
  mutate(i = row_number(), half = if_else(i <= n() / 2, "first half", "second half"),
         parity = if_else(i %% 2 == 1, "odd events", "even events")) %>%
  pivot_longer(c(half, parity), values_to = "subset") %>%
  group_by(route, subset) %>%
  summarise(n = n(), resultant = round(resultant(dir), 4),
            chance = round(1 / sqrt(n()), 4),
            cor_cos_sin = round(collinearity(dir), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)


# ============================================================================
#  3. Does the grid orientation repeat when the route is held fixed?
# ----------------------------------------------------------------------------
#  Same subject, same route, different session: the person, the room, the
#  headings, their order and their timing are all identical. A real
#  environment-anchored grid orientation has to survive that.
# ============================================================================

if (!file.exists(GRID_RUNS)) {
  cat("\n(", GRID_RUNS, " not found — run import_gridcat.R first for section 3)\n",
      sep = "")
} else {
  grid <- readr::read_csv(GRID_RUNS, show_col_types = FALSE) %>%
    filter(hemisphere == HEMI, run != "avg") %>%
    transmute(subject = str_remove(subject_id, "^sub-"),
              session, run = as.integer(as.character(run)), phi = ori_weighted)

  d <- routes %>%
    select(subject, session, run, route) %>%
    inner_join(grid, by = c("subject", "session", "run"))

  pairs_of <- function(df, keys) {
    df %>% group_by(across(all_of(keys))) %>% filter(n() > 1) %>%
      summarise(p = list(as.vector(ori_diff(outer(phi, phi, "-"), 0)[lower.tri(diag(n()))])),
                .groups = "drop") %>%
      unnest(p) %>% pull(p)
  }

  report <- function(label, v) {
    if (!length(v)) return(invisible())
    tt <- t.test(v, mu = NULL_MEAN)
    cat(sprintf("%-46s n=%3d  mean |dphi| = %5.2f deg  p vs chance = %.2f\n",
                label, length(v), mean(v), tt$p.value))
  }

  cat("\n== does the orientation repeat? (chance = 15 deg) =========\n")
  report("same subject, same route, other session",
         pairs_of(d, c("subject", "route")))
  report("same subject, same session, other route",
         pairs_of(d, c("subject", "session")))
  report("same route, different subjects",
         pairs_of(d, c("route")))
}
cat("\n")
