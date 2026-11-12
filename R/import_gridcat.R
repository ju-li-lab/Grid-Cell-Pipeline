# ============================================================================
#  import_gridcat.R — import one or more GridCAT output directories
# ----------------------------------------------------------------------------
#  Each GridCAT output directory (GLM_output_<variant>) holds:
#    *_grid_metrics.txt   GridCAT's grid metrics per subject x session x ROI
#    *_voxel_ori.csv      per-voxel grid orientation per run
#    run_info.tsv         the settings that produced the run (one row)
#
#  Edit GRIDCAT_DIRS, then source the whole file. You get four data frames:
#    run_info        one row per GridCAT run: every setting
#    grid_runs       one row per gridcat_run x subject x session x hemisphere x run
#    grid_stability  one row per run-PAIR (within-voxel stability)
#    voxel_ori       one row per voxel x run (for the stability spaghetti)
#
#  All tables carry `gridcat_run`, so a single run can be analysed alone and
#  several runs can be compared against each other. Everything subject-level
#  (group, age, sex, ...) lives in subject_info.csv and is joined on
#  subject_id; session-level scores can later be joined the same way on
#  subject_id + session.
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)


# ---- settings ---------------------------------------------------------------

# One or more GridCAT output directories. Names become the gridcat_run label;
# unnamed entries are labelled by their folder name.
GRIDCAT_DIRS <- c(
  "/Users/juli/Desktop/SOLIS/GLM_output_6f_t5_am_keep1_avg_1"
)

SUBJECT_INFO_CSV <- "subject_info.csv"

SAVE_CSV   <- TRUE
OUTPUT_DIR <- "."

# Matches e.g. sub-01c06_ses-01_ErC-bilat_grid_metrics.txt
FILE_PATTERN <- "^sub-(\\d+)([a-zA-Z])(\\d+)_ses-(\\d+)_([A-Za-z]+)-([A-Za-z]+)_grid_metrics\\.txt$"


# ---- filename parsing -------------------------------------------------------

parse_filename <- function(filename) {
  m <- str_match(filename, FILE_PATTERN)
  if (is.na(m[1, 1])) return(NULL)
  tibble(
    subject_id   = paste0("sub-", m[1, 2], m[1, 3], m[1, 4]),
    subject_num  = as.integer(m[1, 2]),
    subject_code = tolower(m[1, 3]),
    cohort       = m[1, 4],
    session      = as.integer(m[1, 5]),
    roi          = m[1, 6],
    hemisphere   = tolower(m[1, 7])
  )
}


# ---- grid_metrics files -----------------------------------------------------
#  Four ";"-separated blocks per file, each introduced by "GRID METRIC;".

read_grid_metrics_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- sub(";+$", "", lines)

  starts <- which(startsWith(lines, "GRID METRIC;"))
  if (length(starts) == 0) return(list())
  ends <- c(starts[-1] - 1L, length(lines))

  blocks <- list()
  for (i in seq_along(starts)) {
    body <- lines[(starts[i] + 1L):ends[i]]
    if (length(body) == 0) next
    df <- utils::read.csv(
      text = paste(c(lines[starts[i]], body), collapse = "\n"),
      sep = ";", check.names = FALSE, stringsAsFactors = FALSE
    )
    tag <- df[["GRID METRIC"]][1]
    key <- case_when(
      str_starts(tag, "Magnitude")             ~ "magnitude",
      str_starts(tag, "Between-voxel")         ~ "between_voxel",
      str_starts(tag, "Within-voxel")          ~ "within_voxel",
      str_starts(tag, "Mean grid orientation") ~ "mean_ori",
      .default = NA_character_
    )
    if (!is.na(key)) blocks[[key]] <- df
  }
  blocks
}

# GridCAT writes the run three different ways; everything becomes "1","2",...,"avg".
normalise_run <- function(x) {
  x <- as.character(x)
  n <- str_match(x, regex("run\\s*(\\d+)", ignore_case = TRUE))[, 2]
  out <- ifelse(!is.na(n), n, NA_character_)
  out[str_detect(x, regex("allRuns|averaged", ignore_case = TRUE))] <- "avg"
  plain <- is.na(out) & str_detect(x, "^\\d+$")
  out[plain] <- x[plain]
  out
}

tidy_one_file <- function(path) {
  ids <- parse_filename(basename(path))
  if (is.null(ids)) return(NULL)
  blocks <- read_grid_metrics_file(path)
  if (length(blocks) == 0) return(NULL)

  magnitude <- NULL
  if (!is.null(blocks$magnitude)) {
    b <- blocks$magnitude
    magnitude <- tibble(
      run          = normalise_run(b[["CONTRAST NAME"]]),
      contrast     = as.numeric(b[["MEAN CON-VALUE WITHIN ROI"]]),
      n_voxels     = as.integer(b[["VOXELS WITHIN ROI"]]),
      n_nan_voxels = as.integer(b[["NaN VOXELS WITHIN ROI"]])
    ) %>%
      mutate(pct_nan = if_else(!is.na(n_voxels) & n_voxels > 0,
                               100 * n_nan_voxels / n_voxels, NA_real_))
  }

  between <- NULL
  if (!is.null(blocks$between_voxel)) {
    b <- blocks$between_voxel
    between <- tibble(
      run        = normalise_run(b[["VOXEL-WISE GRID ORI IMAGE"]]),
      rayleigh_z = as.numeric(b[["RAYLEIGH z"]]),
      rayleigh_p = as.numeric(b[["RAYLEIGH p"]])
    )
  }

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

  run_level <- list(magnitude, between, mean_ori) %>%
    compact() %>%
    reduce(full_join, by = "run") %>%
    filter(!is.na(run))
  run_level <- bind_cols(ids[rep(1, nrow(run_level)), ], run_level)

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


# ---- one GridCAT output directory -------------------------------------------

read_run_info <- function(dir) {
  path <- file.path(dir, "run_info.tsv")
  if (!file.exists(path)) {
    warning("No run_info.tsv in ", dir, " — settings columns will be NA.",
            call. = FALSE)
    return(tibble(gridcat_run = basename(sub("/+$", "", dir))))
  }
  read_tsv(path, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
    mutate(across(c(any_of(c("xfold_symmetry", "event_usage_glm1",
                             "event_usage_glm2", "masking_threshold",
                             "keep_unused_grid_events", "avg_ori_across_runs",
                             "use_voxel_weighting", "hpf_cutoff", "smooth_fwhm",
                             "slice_timing", "roi_tsnr_filter"))),
                  as.numeric))
}

read_voxel_ori <- function(dir, label) {
  files <- list.files(dir, pattern = "_voxel_ori\\.csv$", full.names = TRUE)
  files <- files[!startsWith(basename(files), "._")]
  if (length(files) == 0) return(NULL)
  map(files, ~ read_csv(.x, show_col_types = FALSE)) %>%
    list_rbind() %>%
    mutate(
      gridcat_run = label,
      subject_id  = subject,
      session     = as.integer(str_remove(session, "^ses-")),
      hemisphere  = tolower(str_extract(roi, "(?<=-)[A-Za-z]+$")),
      roi         = str_remove(roi, "-[A-Za-z]+$"),
      run         = as.integer(run),
      .keep = "unused"
    )
}

import_one_dir <- function(dir, label, verbose = TRUE) {
  if (!dir.exists(dir)) stop("Directory not found: ", dir)
  files <- list.files(dir, pattern = "_grid_metrics\\.txt$", full.names = TRUE)
  if (length(files) == 0) stop("No *_grid_metrics.txt files found in: ", dir)
  if (verbose) message("[", label, "] reading ", length(files), " files ...")

  parsed <- map(files, function(p) {
    tryCatch(tidy_one_file(p), error = function(e) {
      warning("Could not read ", basename(p), ": ", conditionMessage(e), call. = FALSE)
      NULL
    })
  }) %>% compact()

  runs      <- map(parsed, "runs")      %>% compact() %>% list_rbind()
  stability <- map(parsed, "stability") %>% compact() %>% list_rbind()
  runs$gridcat_run <- label
  if (nrow(stability)) stability$gridcat_run <- label

  list(
    info      = read_run_info(dir) %>% mutate(gridcat_run = label),
    runs      = runs,
    stability = stability,
    voxel_ori = read_voxel_ori(dir, label)
  )
}

import_gridcat <- function(dirs, verbose = TRUE) {
  labels <- names(dirs)
  if (is.null(labels)) labels <- rep("", length(dirs))
  labels[labels == ""] <- basename(sub("/+$", "", dirs[labels == ""]))

  parts <- map2(dirs, labels, import_one_dir, verbose = verbose)

  runs      <- map(parts, "runs")      %>% list_rbind()
  stability <- map(parts, "stability") %>% compact() %>% list_rbind()
  voxel_ori <- map(parts, "voxel_ori") %>% compact() %>% list_rbind()
  info      <- map(parts, "info")      %>% list_rbind() %>%
    relocate(gridcat_run)

  # Sortable run factor: 1, 2, 3, ..., then "avg"
  all_runs <- c(runs$run, stability$run_a, stability$run_b)
  nums <- sort(unique(suppressWarnings(as.integer(all_runs[all_runs != "avg"]))))
  lv   <- c(as.character(nums), "avg")

  runs <- runs %>%
    mutate(
      run        = factor(run, levels = lv),
      run_num    = suppressWarnings(as.integer(as.character(run))),
      is_average = run == "avg",
      hemisphere = factor(hemisphere, levels = c("left", "right", "bilat"))
    ) %>%
    relocate(gridcat_run) %>%
    arrange(gridcat_run, subject_id, session, hemisphere, run)

  if (nrow(stability)) {
    stability <- stability %>%
      mutate(
        run_a      = factor(run_a, levels = lv),
        run_b      = factor(run_b, levels = lv),
        hemisphere = factor(hemisphere, levels = c("left", "right", "bilat"))
      ) %>%
      relocate(gridcat_run) %>%
      arrange(gridcat_run, subject_id, session, hemisphere, run_a, run_b)
  }

  list(info = info, runs = runs, stability = stability, voxel_ori = voxel_ori)
}


# ---- subject-level info -------------------------------------------------------
#  Group and anything else subject-level lives in ONE editable CSV, joined on
#  subject_id. Add columns (age, sex, exclude, ...) and they all come along.

make_subject_info_template <- function(runs, path) {
  tpl <- runs %>%
    distinct(subject_id, subject_num, subject_code, cohort) %>%
    arrange(subject_num, subject_code, cohort) %>%
    mutate(group = NA_character_)
  write_csv(tpl, path, na = "")
  message("Wrote subject-info template: ", normalizePath(path),
          "\n  -> fill in `group` (and any other columns), then re-run.")
  tpl
}

attach_subject_info <- function(df, path) {
  if (is.null(df) || !file.exists(path)) return(df)
  if (!"subject_id" %in% names(df)) return(df)

  info <- read_csv(path, show_col_types = FALSE)
  if (!"subject_id" %in% names(info)) {
    warning("`", path, "` has no `subject_id` column — not joined.", call. = FALSE)
    return(df)
  }
  info <- info %>% select(subject_id, !any_of(setdiff(names(df), "subject_id")))

  unknown <- setdiff(unique(df$subject_id), info$subject_id)
  if (length(unknown)) {
    warning("No row in ", path, " for: ", paste(unknown, collapse = ", "),
            " — those rows get NA.", call. = FALSE)
  }
  left_join(df, info, by = "subject_id") %>%
    relocate(names(info)[names(info) != "subject_id"], .after = subject_id)
}


# ---- run the import -----------------------------------------------------------

imported       <- import_gridcat(GRIDCAT_DIRS)
run_info       <- imported$info
grid_runs      <- imported$runs
grid_stability <- imported$stability
voxel_ori      <- imported$voxel_ori

if (!file.exists(SUBJECT_INFO_CSV)) {
  make_subject_info_template(grid_runs, SUBJECT_INFO_CSV)
}
grid_runs      <- attach_subject_info(grid_runs,      SUBJECT_INFO_CSV)
grid_stability <- attach_subject_info(grid_stability, SUBJECT_INFO_CSV)
voxel_ori      <- attach_subject_info(voxel_ori,      SUBJECT_INFO_CSV)

cat("\n------------------------------------------------------------\n")
cat("  gridcat runs:   ", paste(run_info$gridcat_run, collapse = ", "), "\n", sep = "")
cat("  grid_runs:      ", nrow(grid_runs), " rows\n", sep = "")
cat("  grid_stability: ", nrow(grid_stability), " rows\n", sep = "")
cat("  voxel_ori:      ", if (is.null(voxel_ori)) 0 else nrow(voxel_ori), " rows\n", sep = "")
cat("  subjects:       ", n_distinct(grid_runs$subject_id), "\n", sep = "")
cat("  sessions:       ", paste(sort(unique(grid_runs$session)), collapse = ", "), "\n", sep = "")
cat("------------------------------------------------------------\n\n")

# Sessions with an unusual number of runs are worth knowing about.
one_hemi <- levels(droplevels(grid_runs$hemisphere))[1]
run_counts <- grid_runs %>%
  filter(!is_average, hemisphere == one_hemi) %>%
  count(gridcat_run, subject_id, session, name = "n_runs")
if (n_distinct(run_counts$n_runs) > 1) {
  cat("NOTE — these subject/sessions do not have the usual number of runs:\n")
  print(as.data.frame(run_counts %>% filter(n_runs != max(n_runs))))
  cat("\n")
}

if (isTRUE(SAVE_CSV)) {
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  write_csv(run_info,       file.path(OUTPUT_DIR, "run_info.csv"))
  write_csv(grid_runs,      file.path(OUTPUT_DIR, "grid_runs.csv"))
  write_csv(grid_stability, file.path(OUTPUT_DIR, "grid_stability.csv"))
  if (!is.null(voxel_ori)) write_csv(voxel_ori, file.path(OUTPUT_DIR, "voxel_ori.csv"))
  message("Saved CSVs to ", normalizePath(OUTPUT_DIR))
}


# ---- what you get -------------------------------------------------------------
#  grid_runs: gridcat_run, subject_id (+ subject_info columns), session, roi,
#    hemisphere, run ("1","2",...,"avg"), run_num, is_average, contrast,
#    n_voxels, n_nan_voxels, pct_nan, rayleigh_z, rayleigh_p, ori_weighted,
#    ori_unweighted, grid_event.
#  grid_stability: same keys per run PAIR: run_a, run_b, comparison, pct_stable.
#  voxel_ori: gridcat_run, subject_id, session, roi, hemisphere, run, voxel,
#    x, y, z, ori_deg — one row per voxel per run.
#  run_info: one row per gridcat_run with every setting; join it onto any
#    table with left_join(..., by = "gridcat_run").
