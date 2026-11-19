# ============================================================================
#  import_spatial_cognition.R — import spatial-cognition test results (.xlsx)
# ----------------------------------------------------------------------------
#  A workbook has one sheet per test (SUT, 4MTN, PtO, Navis), one row per
#  subject x session x condition, and the performance measures as columns.
#  Its subject codes lack the study code (sub-01, not sub-01c06) and its
#  session numbers are not the MRI sessions, so every import names both:
#
#    import_spatial_cog(files, study = "c06", sessions = c("0" = 1, "9" = 2))
#
#  `study` is appended to each subject code. `sessions` reads
#  "file session = MRI session"; file sessions you leave out stay in the
#  table with session = NA, and several file sessions may point at the same
#  MRI session (the analysis then averages them). All files of one study can
#  go into one call; stack one call per study at the bottom of this file.
#
#  Result: spatial_cog, one row per subject x file session x test x
#  condition x measure.
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(readxl)

SAVE_CSV   <- TRUE
OUTPUT_DIR <- "."

# Columns that are neither a condition nor a performance measure.
ID_COLUMNS <- c("Subject", "Session", "Version", "Date", "Time", "Protocol")

# Columns that are conditions even when they hold numbers. Text columns
# (Difficulty, Axis, Mode, ...) are conditions anyway; every other numeric
# column is a measure.
CONDITION_COLUMNS <- c("Condition", "MovementProb")


# ---- one sheet = one test -----------------------------------------------------

is_number_column <- function(x) {
  all(is.na(x) | !is.na(suppressWarnings(as.numeric(as.character(x)))))
}

read_test_sheet <- function(path, sheet) {
  d <- read_excel(path, sheet = sheet)
  other <- setdiff(names(d), ID_COLUMNS)
  conditions <- other[other %in% CONDITION_COLUMNS |
                        !map_lgl(d[other], is_number_column)]
  measures <- setdiff(other, conditions)

  # "Easy Pitch", "MovementProb=1", ...; "all" when the test has no conditions
  d <- d %>%
    mutate(across(all_of(conditions) & where(is.numeric),
                  ~ if_else(is.na(.x), NA_character_, paste0(cur_column(), "=", .x))),
           across(all_of(conditions), as.character))
  if (length(conditions)) {
    d <- unite(d, "condition", all_of(conditions), sep = " ", na.rm = TRUE)
  }
  if (!"condition" %in% names(d)) d$condition <- ""

  d %>%
    transmute(
      subject    = str_trim(as.character(Subject)),
      sc_session = as.integer(parse_number(as.character(Session))),
      date       = as.Date(Date),
      test       = sheet,
      condition  = if_else(condition == "", "all", condition),
      across(all_of(measures), ~ suppressWarnings(as.numeric(as.character(.x))))
    ) %>%
    pivot_longer(all_of(measures), names_to = "measure", values_to = "value") %>%
    # Navis writes -1 when no target was revisited: that is not a time
    mutate(value = if_else(measure == "TimeToFirstRevisit" & value < 0, NA, value)) %>%
    filter(!is.na(value))
}


# ---- one study ----------------------------------------------------------------

import_spatial_cog <- function(files, study, sessions) {
  files <- path.expand(files)
  d <- map(files, function(f) {
    map(excel_sheets(f), ~ read_test_sheet(f, .x)) %>%
      list_rbind() %>%
      mutate(file = basename(f))
  }) %>%
    list_rbind()

  d %>%
    mutate(
      subject_id = paste0(str_replace(subject, "^sub-(\\d)$", "sub-0\\1"), study),
      study      = study,
      session    = as.integer(unname(sessions[as.character(sc_session)]))
    ) %>%
    select(study, subject_id, sc_session, session, date, everything(), -subject)
}


# ---- run the import -------------------------------------------------------------
#  One call per study. CHECK the session mappings against the MRI dates; the
#  table printed below lists every file session with its dates.

spatial_cog <- bind_rows(
  import_spatial_cog("~/Downloads/SPCGN1_SIRIUSCONTROL_summary_2026-05-11.xlsx",
                     study = "c06", sessions = c("0" = 1, "9" = 2)),
  import_spatial_cog("~/Downloads/SPCGN1_SIRIUSCREW_summary_2026-05-12_1.xlsx",
                     study = "s06", sessions = c("0" = 1, "9" = 2)),
  import_spatial_cog("~/Downloads/SOLIS100_SPCGN1_summary_2026-08-25.xlsx",
                     study = "s17", sessions = c("BDC-5" = 1, "MD90" = 2)),
)

cat("\n== file sessions and the MRI session they stand for =====================\n")
spatial_cog %>%
  group_by(study, sc_session, session) %>%
  summarise(subjects = n_distinct(subject_id),
            from = min(date), to = max(date), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n== tests, conditions and measures ===========================================\n")
for (t in unique(spatial_cog$test)) {
  d <- spatial_cog[spatial_cog$test == t, ]
  cat(sprintf("%-6s conditions: %s\n       measures:   %s\n", t,
              paste(unique(d$condition), collapse = ", "),
              paste(unique(d$measure), collapse = ", ")))
}

if (exists("grid_runs")) {
  mapped  <- unique(spatial_cog$subject_id[!is.na(spatial_cog$session)])
  no_grid <- setdiff(mapped, grid_runs$subject_id)
  no_sc   <- setdiff(unique(grid_runs$subject_id), mapped)
  if (length(no_grid)) cat("\nNo GridCAT data for:", paste(no_grid, collapse = ", "), "\n")
  if (length(no_sc))   cat("No spatial-cognition data for:", paste(no_sc, collapse = ", "), "\n")
}
cat("\n")

if (isTRUE(SAVE_CSV)) {
  write_csv(spatial_cog, file.path(OUTPUT_DIR, "spatial_cog.csv"))
  message("Saved spatial_cog.csv to ", normalizePath(OUTPUT_DIR))
}


# ---- what you get -----------------------------------------------------------------
#  spatial_cog: study, subject_id (sub-01c06), sc_session (as in the file),
#    session (MRI session, NA if not mapped), date, test (sheet name),
#    condition ("all" if the test has none), measure, value, file.
