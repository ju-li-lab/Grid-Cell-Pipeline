# ============================================================================
#  import_gridcat.R  —  Mass-import GridCAT grid-metric files into R
# ----------------------------------------------------------------------------
#  Reads every *_grid_metrics.txt file in ONE GridCAT output directory and
#  returns tidy data frames organised by subject / session / run / hemisphere.
#
#  HOW TO RUN (RStudio):
#    1. Edit the SETTINGS block below (mainly DATA_DIR).
#    2. Source the whole file  (Ctrl/Cmd + Shift + S).
#    3. You get two data frames in your environment:
#         grid_runs       one row per subject x session x hemisphere x run
#         grid_stability  one row per run-PAIR (within-voxel stability)
#
#  Nothing below the SETTINGS block needs editing for normal use.
# ============================================================================


# ============================================================================
#  SETTINGS  —  edit these
# ============================================================================

# >>> 1. Which directory to import (one GLM output variant at a time).
DATA_DIR <- "/Volumes/juli_ssd/Archive/GLM_output_7mm_6f_am_keep1_avg0"

# >>> 2. Subject-level info (group, age, ...).
#     A CSV with one row per subject. If it does not exist yet, the script
#     writes a TEMPLATE with an empty `group` column for you to fill in.
SUBJECT_INFO_CSV <- "subject_info.csv"

# >>> 3. Where to write the imported CSVs. Set SAVE_CSV <- FALSE to skip.
#     Note: the filenames are fixed, so importing a second variant overwrites
#     the first. Add a prefix here if you want to keep several side by side.
SAVE_CSV   <- FALSE
OUTPUT_DIR <- "."

# >>> 4. Filename pattern.
#     Matches e.g.  sub-01c06_ses-01_ErC-bilat_grid_metrics.txt
#                        ^^ ^^^     ^^   ^^^ ^^^^^
#                        |  |       |    |   +-- hemisphere
#                        |  |       |    +------ roi label
#                        |  |       +----------- session
#                        |  +------------------- subject code + cohort
#                        +---------------------- subject number
FILE_PATTERN <- "^sub-(\\d+)([a-zA-Z])(\\d+)_ses-(\\d+)_([A-Za-z]+)-([A-Za-z]+)_grid_metrics\\.txt$"


# ============================================================================
#  PACKAGES
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)


# ============================================================================
#  1. Parse a filename into its identifying fields
# ============================================================================

parse_filename <- function(filename) {
  m <- str_match(filename, FILE_PATTERN)
  if (is.na(m[1, 1])) return(NULL)

  tibble(
    subject_id   = paste0("sub-", m[1, 2], m[1, 3], m[1, 4]),  # sub-01c06
    subject_num  = as.integer(m[1, 2]),                        # 1
    subject_code = tolower(m[1, 3]),                           # "c" / "s"
    cohort       = m[1, 4],                                    # "06" / "13"
    session      = as.integer(m[1, 5]),                        # 1, 2, 3, ...
    roi          = m[1, 6],                                    # "ErC"
    hemisphere   = tolower(m[1, 7])                            # bilat/left/right
  )
}


# ============================================================================
#  2. Split one GridCAT file into its four blocks
# ----------------------------------------------------------------------------
#  Each file holds four ";"-separated tables, each introduced by a header
#  line starting with "GRID METRIC;". They are returned as a named list.
# ============================================================================

read_grid_metrics_file <- function(path) {

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]      # drop blank separator lines
  lines <- sub(";+$", "", lines)     # drop trailing ";" (would add a NA column)

  starts <- which(startsWith(lines, "GRID METRIC;"))
  if (length(starts) == 0) return(list())
  ends <- c(starts[-1] - 1L, length(lines))

  blocks <- list()
  for (i in seq_along(starts)) {
    body <- lines[(starts[i] + 1L):ends[i]]
    if (length(body) == 0) next

    df <- utils::read.csv(
      text             = paste(c(lines[starts[i]], body), collapse = "\n"),
      sep              = ";",
      check.names      = FALSE,
      stringsAsFactors = FALSE
    )

    tag <- df[["GRID METRIC"]][1]
    key <- case_when(
      str_starts(tag, "Magnitude")              ~ "magnitude",
      str_starts(tag, "Between-voxel")          ~ "between_voxel",
      str_starts(tag, "Within-voxel")           ~ "within_voxel",
      str_starts(tag, "Mean grid orientation")  ~ "mean_ori",
      .default = NA_character_
    )
    if (!is.na(key)) blocks[[key]] <- df
  }
  blocks
}


# ============================================================================
#  3. Normalise run labels
# ----------------------------------------------------------------------------
#  GridCAT writes the run three different ways depending on the block:
#     "Run2_translate-..."                  -> "2"
#     "voxelwiseGridOri-translate-run2-deg" -> "2"
#     RUN column: 2 / "averaged across runs" -> "2" / "avg"
#  Everything becomes "1", "2", "3", ... or "avg".
# ============================================================================

normalise_run <- function(x) {
  x <- as.character(x)
  n <- str_match(x, regex("run\\s*(\\d+)", ignore_case = TRUE))[, 2]
  out <- ifelse(!is.na(n), n, NA_character_)
  is_avg <- str_detect(x, regex("allRuns|averaged", ignore_case = TRUE))
  out[is_avg] <- "avg"
  # plain numeric RUN column ("1", "2", ...)
  plain <- is.na(out) & str_detect(x, "^\\d+$")
  out[plain] <- x[plain]
  out
}


# ============================================================================
#  4. Turn the blocks of ONE file into tidy rows
# ============================================================================

tidy_one_file <- function(path) {

  ids <- parse_filename(basename(path))
  if (is.null(ids)) return(NULL)

  blocks <- read_grid_metrics_file(path)
  if (length(blocks) == 0) return(NULL)

  # ---- 4a. Contrast magnitude: one row per run + "avg" --------------------
  magnitude <- NULL
  if (!is.null(blocks$magnitude)) {
    b <- blocks$magnitude
    magnitude <- tibble(
      run           = normalise_run(b[["CONTRAST NAME"]]),
      contrast      = as.numeric(b[["MEAN CON-VALUE WITHIN ROI"]]),
      n_voxels      = as.integer(b[["VOXELS WITHIN ROI"]]),
      n_nan_voxels  = as.integer(b[["NaN VOXELS WITHIN ROI"]])
    )
  }

  # ---- 4b. Between-voxel coherence: Rayleigh z / p ------------------------
  between <- NULL
  if (!is.null(blocks$between_voxel)) {
    b <- blocks$between_voxel
    between <- tibble(
      run        = normalise_run(b[["VOXEL-WISE GRID ORI IMAGE"]]),
      rayleigh_z = as.numeric(b[["RAYLEIGH z"]]),
      rayleigh_p = as.numeric(b[["RAYLEIGH p"]])
    )
  }

  # ---- 4c. Mean grid orientation ------------------------------------------
  mean_ori <- NULL
  if (!is.null(blocks$mean_ori)) {
    b <- blocks$mean_ori
    mean_ori <- tibble(
      run            = normalise_run(b[["RUN"]]),
      grid_event     = b[["GRID EVENT"]],
      ori_weighted   = as.numeric(b[["MEAN GRID ORI IN DEGREES (VOXELS WEIGHTED)"]]),
      ori_unweighted = as.numeric(b[["MEAN GRID ORI IN DEGREES (VOXELS NOT WEIGHTED)"]])
    )
  }

  # ---- 4d. Join the three run-level blocks --------------------------------
  run_level <- list(magnitude, between, mean_ori) %>%
    compact() %>%
    reduce(full_join, by = "run") %>%
    filter(!is.na(run))

  run_level <- bind_cols(ids[rep(1, nrow(run_level)), ], run_level)

  # ---- 4e. Within-voxel stability: one row per PAIR of runs ---------------
  stability <- NULL
  if (!is.null(blocks$within_voxel)) {
    b <- blocks$within_voxel
    stability <- tibble(
      run_a      = normalise_run(b[["VOXEL-WISE GRID ORI IMAGE 1"]]),
      run_b      = normalise_run(b[["VOXEL-WISE GRID ORI IMAGE 2"]]),
      pct_stable = as.numeric(b[["% STABLE VOXELS WITHIN ROI"]]),
      n_voxels   = as.integer(b[["IMAGE 1 VOXELS WITHIN ROI"]])
    ) %>%
      filter(!is.na(run_a), !is.na(run_b)) %>%
      mutate(comparison = paste(run_a, run_b, sep = "_vs_"))

    stability <- bind_cols(ids[rep(1, nrow(stability)), ], stability)
  }

  list(runs = run_level, stability = stability)
}


# ============================================================================
#  5. Import the whole directory
# ============================================================================

import_gridcat <- function(dir, verbose = TRUE) {

  if (!dir.exists(dir)) stop("Directory not found: ", dir)

  files <- list.files(dir, pattern = "_grid_metrics\\.txt$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No *_grid_metrics.txt files found in: ", dir)
  }
  if (verbose) message("Reading ", length(files), " files from ", dir, " ...")

  parsed <- map(files, function(p) {
    tryCatch(tidy_one_file(p), error = function(e) {
      warning("Could not read ", basename(p), ": ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
  }) %>% compact()

  skipped <- length(files) - length(parsed)
  if (verbose && skipped > 0) message("  ", skipped, " file(s) skipped.")

  runs      <- map(parsed, "runs")      %>% compact() %>% list_rbind()
  stability <- map(parsed, "stability") %>% compact() %>% list_rbind()

  # Consistent, sortable factor levels: 1, 2, 3, ... then "avg"
  all_runs <- c(runs$run, stability$run_a, stability$run_b)
  nums     <- sort(unique(suppressWarnings(as.integer(all_runs[all_runs != "avg"]))))
  lv       <- c(as.character(nums), "avg")

  runs <- runs %>%
    mutate(
      run        = factor(run, levels = lv),
      run_num    = suppressWarnings(as.integer(as.character(run))), # NA for "avg"
      is_average = run == "avg",
      hemisphere = factor(hemisphere, levels = c("left", "right", "bilat"))
    ) %>%
    arrange(subject_id, session, hemisphere, run)

  if (nrow(stability)) {
    stability <- stability %>%
      mutate(
        run_a      = factor(run_a, levels = lv),
        run_b      = factor(run_b, levels = lv),
        hemisphere = factor(hemisphere, levels = c("left", "right", "bilat"))
      ) %>%
      arrange(subject_id, session, hemisphere, run_a, run_b)
  }

  list(runs = runs, stability = stability)
}


# ============================================================================
#  6. Subject-level info (group, and anything else you want to add later)
# ----------------------------------------------------------------------------
#  The importer deliberately knows nothing about groups. Group membership
#  lives in ONE editable CSV that is joined on `subject_id`. Add as many
#  extra columns as you like (age, sex, dropout, ...) — they all come along.
# ============================================================================

make_subject_info_template <- function(runs, path) {
  tpl <- runs %>%
    distinct(subject_id, subject_num, subject_code, cohort) %>%
    arrange(subject_num, subject_code, cohort) %>%
    mutate(group = NA_character_)        # <- you fill this in

  write_csv(tpl, path, na = "")
  message("Wrote subject-info template: ", normalizePath(path))
  message("  -> Open it, fill in the `group` column, add any other columns ",
          "you need, then re-run this script.")
  tpl
}

attach_subject_info <- function(df, path) {
  if (!file.exists(path)) return(df)
  if (!"subject_id" %in% names(df)) return(df)   # nothing was imported

  info <- read_csv(path, show_col_types = FALSE)
  if (!"subject_id" %in% names(info)) {
    warning("`", path, "` has no `subject_id` column — not joined.",
            call. = FALSE)
    return(df)
  }

  # Keep subject_id plus any column the imported data does not already have,
  # so the CSV is free to carry group, age, sex, ... without clashing.
  info <- info %>% select(subject_id, !any_of(setdiff(names(df), "subject_id")))

  unknown <- setdiff(unique(df$subject_id), info$subject_id)
  if (length(unknown)) {
    warning("No row in ", path, " for: ", paste(unknown, collapse = ", "),
            " — those rows get NA.", call. = FALSE)
  }

  left_join(df, info, by = "subject_id") %>%
    relocate(names(info)[names(info) != "subject_id"], .after = subject_id)
}


# ============================================================================
#  7. RUN THE IMPORT
# ============================================================================

imported       <- import_gridcat(DATA_DIR)
grid_runs      <- imported$runs
grid_stability <- imported$stability

# Create the subject-info template on first run, then join it.
if (!file.exists(SUBJECT_INFO_CSV)) {
  make_subject_info_template(grid_runs, SUBJECT_INFO_CSV)
}
grid_runs      <- attach_subject_info(grid_runs,      SUBJECT_INFO_CSV)
grid_stability <- attach_subject_info(grid_stability, SUBJECT_INFO_CSV)


# ---- Summary ---------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("  directory:      ", basename(sub("/+$", "", DATA_DIR)), "\n", sep = "")
cat("  grid_runs:      ", nrow(grid_runs), " rows\n", sep = "")
cat("  grid_stability: ", nrow(grid_stability), " rows\n", sep = "")
cat("  subjects:       ", dplyr::n_distinct(grid_runs$subject_id), "\n", sep = "")
cat("  sessions:       ", paste(sort(unique(grid_runs$session)), collapse = ", "), "\n", sep = "")
cat("  hemispheres:    ", paste(levels(droplevels(grid_runs$hemisphere)), collapse = ", "), "\n", sep = "")
cat("  runs:           ", paste(levels(droplevels(grid_runs$run)), collapse = ", "), "\n", sep = "")
if ("group" %in% names(grid_runs)) {
  if (all(is.na(grid_runs$group))) {
    cat("  groups:         (not filled in yet — edit ", SUBJECT_INFO_CSV, ")\n", sep = "")
  } else {
    cat("  groups:         ",
        paste(names(table(grid_runs$group)), collapse = ", "), "\n", sep = "")
  }
}
cat("------------------------------------------------------------\n\n")

# How many runs does each subject/session have? Flags aborted or short sessions.
one_hemi <- levels(droplevels(grid_runs$hemisphere))[1]
run_counts <- grid_runs %>%
  filter(!is_average, hemisphere == one_hemi) %>%
  count(subject_id, session, name = "n_runs")
if (dplyr::n_distinct(run_counts$n_runs) > 1) {
  cat("NOTE — these subject/sessions do not have the usual number of runs:\n")
  print(as.data.frame(run_counts %>% filter(n_runs != max(n_runs))))
  cat("\n")
}


# ---- Save ------------------------------------------------------------------

if (isTRUE(SAVE_CSV)) {
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  write_csv(grid_runs,      file.path(OUTPUT_DIR, "grid_runs.csv"))
  write_csv(grid_stability, file.path(OUTPUT_DIR, "grid_stability.csv"))
  message("Saved grid_runs.csv and grid_stability.csv to ",
          normalizePath(OUTPUT_DIR))
}


# ============================================================================
#  8. WHAT YOU GET  /  EXAMPLES
# ----------------------------------------------------------------------------
#
#  grid_runs — one row per subject x session x hemisphere x run
#    subject_id      "sub-01c06"        subject_num   1
#    group, ...      whatever you put in subject_info.csv
#    subject_code    "c" / "s"          cohort        "06" / "13"
#    session         integer: 1, 2, 3, ... (any number of sessions)
#    roi             "ErC"              hemisphere    left / right / bilat
#    run             "1","2","3",...,"avg"            run_num  1,2,3, NA for "avg"
#    is_average      TRUE for the across-runs row
#    contrast        aligned vs misaligned mean con-value   <- main measure
#    n_voxels, n_nan_voxels
#    rayleigh_z, rayleigh_p        between-voxel orientation coherence
#    ori_weighted, ori_unweighted  mean grid orientation in degrees
#    grid_event      "translate"
#
#  grid_stability — one row per PAIR of runs (a different grain, so it is a
#  separate table): run_a, run_b, comparison ("1_vs_2"), pct_stable, n_voxels
#
#  ---- typical uses --------------------------------------------------------
#
#  # per-run data only, bilateral ErC
#  d <- grid_runs %>% filter(!is_average, hemisphere == "bilat")
#
#  # the across-runs summary value per subject/session
#  d_avg <- grid_runs %>% filter(is_average, hemisphere == "bilat")
#
#  # average the two unilateral masks instead of using bilat
#  d_uni <- grid_runs %>%
#    filter(!is_average, hemisphere %in% c("left", "right")) %>%
#    group_by(subject_id, group, session, run) %>%
#    summarise(contrast = mean(contrast), .groups = "drop")
#
#  # wide: one column per session
#  grid_runs %>%
#    filter(is_average, hemisphere == "bilat") %>%
#    select(subject_id, group, session, contrast) %>%
#    pivot_wider(names_from = session, values_from = contrast,
#                names_prefix = "ses")
#
#  # group means per session
#  d %>% group_by(group, session) %>%
#    summarise(m = mean(contrast), sd = sd(contrast), n = n(), .groups = "drop")
#
#  # mixed model (needs lme4 / lmerTest). `session` is an integer, so use
#  # factor(session) for categorical sessions, or leave it numeric for a
#  # linear trend across sessions.
#  # lmerTest::lmer(contrast ~ group * factor(session) + (1 | subject_id), data = d)
#
#  ---- ADDING GROUP INFO ---------------------------------------------------
#  Everything subject-level lives in subject_info.csv. The first run writes a
#  template with an empty `group` column:
#
#     subject_id,subject_num,subject_code,cohort,group
#     sub-01c06,1,c,06,
#     sub-01s06,1,s,06,
#
#  Fill in `group`, add any other columns (age, sex, exclude, ...) and re-run
#  this script — every column is joined onto both tables automatically.
#  Subjects missing from the CSV get NA and a warning, so nothing is dropped
#  silently. The template is never overwritten once it exists.
# ============================================================================
