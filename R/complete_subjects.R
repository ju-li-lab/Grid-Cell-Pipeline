# ============================================================================
#  complete_subjects()  —  balanced subject sets
# ----------------------------------------------------------------------------
#  Returns the subjects that have a value in EVERY session x run cell, and
#  (if a stability frame is given) every run-pair as well.
#
#  Why this matters: a subject who ran well in session 1 but dropped out of
#  session 3 makes session 1 look better than session 3 purely by being in
#  it. Any comparison BETWEEN sessions needs the same people in each.
#
#  Pass frames that are already filtered to one hemisphere and to whatever
#  sessions you care about, with the "avg" rows removed. The grid a subject
#  has to fill is built from what you pass in, so looking at one session
#  asks only for that session's runs.
#
#  A row that exists but carries NA counts as missing - an ROI where every
#  voxel was NaN tells you nothing and must not be averaged as a zero.
#
#  This lives in its own file, rather than inside import_gridcat.R, because
#  import_gridcat.R gets edited every time DATA_DIR changes and restoring an
#  older copy of it would silently take this function with it. Nothing here
#  needs editing.
#
#  group_summary.R and session_effects.R source this themselves - you do not
#  need to run it by hand.
# ============================================================================

library(dplyr)

complete_subjects <- function(runs, stab = NULL, value_col = "contrast",
                              verbose = TRUE) {

  if (!value_col %in% names(runs)) stop("No column called '", value_col, "'.")

  cells_want <- runs %>% distinct(session, run) %>% nrow()
  have_runs  <- runs %>%
    filter(!is.na(.data[[value_col]])) %>%
    distinct(subject_id, session, run) %>%
    count(subject_id, name = "cells")

  keep <- have_runs %>% filter(cells == cells_want) %>% pull(subject_id)

  pairs_want <- 0
  have_pairs <- NULL
  if (!is.null(stab) && nrow(stab)) {
    pairs_want <- stab %>% distinct(session, comparison) %>% nrow()
    have_pairs <- stab %>%
      filter(!is.na(pct_stable)) %>%
      distinct(subject_id, session, comparison) %>%
      count(subject_id, name = "pairs")
    keep <- intersect(keep,
                      have_pairs %>% filter(pairs == pairs_want) %>% pull(subject_id))
  }

  all_subs <- sort(unique(runs$subject_id))
  dropped  <- setdiff(all_subs, keep)

  if (verbose) {
    cat(sprintf("\n=== completeness ===\neach subject needs %d session x run cells",
                cells_want))
    if (pairs_want > 0) cat(sprintf(" and %d run-pairs", pairs_want))
    cat(sprintf("\nkept %d of %d subjects\n", length(keep), length(all_subs)))

    if (length(dropped)) {
      fill0 <- function(x) { x[is.na(x)] <- 0L; x }
      tbl <- data.frame(
        subject_id = dropped,
        cells      = fill0(have_runs$cells[match(dropped, have_runs$subject_id)]),
        of_cells   = cells_want,
        stringsAsFactors = FALSE)
      if (!is.null(have_pairs)) {
        tbl$pairs    <- fill0(have_pairs$pairs[match(dropped, have_pairs$subject_id)])
        tbl$of_pairs <- pairs_want
      }
      cat("\ndropped:\n")
      print(tbl, row.names = FALSE)
    }
  }

  if (!length(keep)) {
    stop("No subject has a complete set. Look at the table above, then ",
         "either fix the gaps or set REQUIRE_COMPLETE <- FALSE.")
  }
  keep
}
