# ============================================================================
#  rayleigh.R  —  What the Rayleigh z in the GridCAT output actually measures
# ----------------------------------------------------------------------------
#  GridCAT reports it under "Between-voxel grid orientation coherence within
#  ROI", computed on the voxelwiseGridOri-*-runN-deg image. One angle per
#  voxel; the time series is already collapsed by then.
#
#  The statistic is CircStat's circ_rtest:
#
#       rbar = | sum_v exp(i * 6 * phi_v) | / n        mean resultant length
#       z    = n * rbar^2                              <- what GridCAT prints
#       p    = exp( sqrt(1 + 4n + 4(n^2 - R^2)) - (1 + 2n) ),   R = n * rbar
#
#  so z runs from 0 (voxels point everywhere) to n (all voxels identical), and
#  n is the number of USABLE voxels, i.e. VOXELS WITHIN ROI minus NaN. This
#  script verifies that reading against the p-values in the files, then asks
#  whether a high z means anything for the contrast.
# ============================================================================

library(dplyr); library(tidyr); library(stringr); library(purrr)

DIRS <- c(halves  = "/Users/juli/Desktop/SOLIS/GLM_output_pmod_avg0_halves_oap_viarest_thr07_smooth4",
          oddeven = "/Users/juli/Desktop/SOLIS/GLM_output_pmod_avg0_halves_oap_viarest_thr07_oddeven")
HEMI <- "left"
VOXEL_MM3    <- 2 * 2 * 3.45   # the EPI grid these ROIs live on
SMOOTH_FWHM  <- 5              # config.sh; the folder name says 4, so check

read_blocks <- function(path) {
  L <- trimws(readLines(path, warn = FALSE)); L <- sub(";+$", "", L[nzchar(L)])
  st <- which(startsWith(L, "GRID METRIC;")); en <- c(st[-1] - 1L, length(L))
  out <- list()
  for (i in seq_along(st)) {
    df  <- read.csv(text = paste(L[st[i]:en[i]], collapse = "\n"), sep = ";",
                    check.names = FALSE)
    tag <- df[[1]][1]
    key <- if (startsWith(tag, "Magnitude")) "mag" else
           if (startsWith(tag, "Between"))   "bet" else
           if (startsWith(tag, "Within"))    "wit" else "ori"
    out[[key]] <- df
  }
  out
}
run_label <- function(x) {
  x <- as.character(x); n <- str_match(x, regex("run\\s*(\\d+)", TRUE))[, 2]
  ifelse(!is.na(n), n,
    ifelse(str_detect(x, regex("allRuns|averaged", TRUE)), "avg",
      ifelse(str_detect(x, "^\\d+$"), x, NA_character_)))
}
read_dir <- function(dir, split) {
  map_dfr(list.files(dir, pattern = "_grid_metrics\\.txt$", full.names = TRUE),
    function(p) {
      id <- str_match(basename(p), "^sub-(\\d+[a-z]\\d+)_ses-(\\d+)_[A-Za-z]+-([A-Za-z]+)_")
      b  <- read_blocks(p)
      list(
        tibble(run = run_label(b$bet[["VOXEL-WISE GRID ORI IMAGE"]]),
               z   = as.numeric(b$bet[["RAYLEIGH z"]]),
               p   = as.numeric(b$bet[["RAYLEIGH p"]]),
               n_voxels = as.integer(b$bet[["VOXELS WITHIN ROI"]]),
               n_nan    = as.integer(b$bet[["NaN VOXELS WITHIN ROI"]])),
        tibble(run = run_label(b$mag[["CONTRAST NAME"]]),
               contrast = as.numeric(b$mag[["MEAN CON-VALUE WITHIN ROI"]])),
        tibble(run = run_label(b$ori[["RUN"]]),
               phi = as.numeric(b$ori[["MEAN GRID ORI IN DEGREES (VOXELS WEIGHTED)"]]))
      ) %>% reduce(full_join, by = "run") %>% filter(!is.na(run)) %>%
        mutate(subject = id[2], session = as.integer(id[3]), hemisphere = id[4],
               split = split, usable = n_voxels - n_nan, .before = 1)
    })
}

A <- imap_dfr(DIRS, ~ read_dir(.x, .y)) %>% filter(run != "avg")
L <- A %>% filter(hemisphere == HEMI)

# Zar's approximation, the one circ_rtest uses.
zar_p <- function(R, n) exp(sqrt(1 + 4*n + 4*(n^2 - R^2)) - (1 + 2*n))
rbar  <- function(z, n) sqrt(z / n)


# ============================================================================
#  1. Confirm the formula and the n it uses
# ============================================================================

cat("\n== is n the usable voxels, or all of them? ================\n")
chk <- A %>% filter(!is.na(z), !is.na(p), usable > 0) %>%
  mutate(err_usable = abs(zar_p(sqrt(z * usable),   usable)   - p),
         err_all    = abs(zar_p(sqrt(z * n_voxels), n_voxels) - p))
cat(sprintf("  n = usable : median error vs the reported p = %.5f\n", median(chk$err_usable)))
cat(sprintf("  n = all    : median error vs the reported p = %.5f\n", median(chk$err_all)))
cat(sprintf("  rows checked %d. (residual is the files' 3-decimal rounding)\n", nrow(chk)))
cat(sprintf("  any z above its ceiling n? %s (largest z/n = %.3f)\n",
            if (any(chk$z > chk$usable + 1e-6)) "YES" else "no", max(chk$z / chk$usable)))


# ============================================================================
#  2. The ceiling, per ROI size
# ----------------------------------------------------------------------------
#  z has no fixed scale: it tops out at n. So it is NOT comparable between runs
#  whose ROIs differ in size. rbar = sqrt(z/n) is.
# ============================================================================

cat("\n== z cannot exceed the voxel count =======================\n")
L %>% group_by(usable) %>%
  summarise(runs = n(), z_ceiling = first(usable),
            z_mean = round(mean(z), 2), z_max = round(max(z), 2),
            rbar_mean = round(mean(rbar(z, usable)), 3),
            rbar_if_independent = round(1 / sqrt(first(usable)), 3),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)


# ============================================================================
#  3. Is the test's null defensible here?
# ----------------------------------------------------------------------------
#  Rayleigh assumes n INDEPENDENT angles. Smoothing makes neighbouring voxels
#  near-copies, so the effective n is far below the voxel count and z is
#  inflated with no grid code present.
# ============================================================================

cat("\n== against the null (z ~ Exp(1), mean 1) =================\n")
cat(sprintf("  observed mean z   %.2f\n", mean(L$z)))
cat(sprintf("  runs at p < .05   %.0f%%   (null: 5%%)\n",  100 * mean(L$p < .05)))
cat(sprintf("  runs at p < .001  %.0f%%   (null: 0.1%%)\n", 100 * mean(L$p < .001)))
m <- lm(z ~ usable, L)
cat(sprintf("  z grows %.2f per extra voxel (p = %.3g) — independent voxels would give 0\n",
            coef(m)[2], summary(m)$coefficients[2, 4]))
resel_vox <- SMOOTH_FWHM^3 / VOXEL_MM3
cat(sprintf("\n  at FWHM %g mm one resel = %.1f voxels of %.1f mm3,\n",
            SMOOTH_FWHM, resel_vox, VOXEL_MM3))
cat(sprintf("  so these ROIs span %.2f to %.2f independent elements, not %d to %d.\n",
            min(L$usable)/resel_vox, max(L$usable)/resel_vox,
            min(L$usable), max(L$usable)))


# ============================================================================
#  4. Does a high z buy a usable orientation, or a contrast?
# ============================================================================

W <- L %>% select(subject, session, run, split, z, contrast, phi, usable) %>%
  pivot_wider(names_from = split, values_from = c(z, contrast, phi, usable)) %>%
  mutate(d_phi = abs(((phi_halves - phi_oddeven + 30) %% 60) - 30),
         z_min = pmin(z_halves, z_oddeven))

cat("\n== high z vs the orientation's stability =================\n")
cat(sprintf("  cor(z, |d_phi| between the two splits) = %+.2f, p = %.3f\n",
            cor(W$z_min, W$d_phi), cor.test(W$z_min, W$d_phi)$p.value))
W %>% mutate(tertile = ntile(z_min, 3)) %>% group_by(tertile) %>%
  summarise(z_range = paste(round(range(z_min), 1), collapse = " - "),
            mean_d_phi = round(mean(d_phi), 1), n = n(), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("  (chance |d_phi| is 15 deg, and 15 deg is also where the contrast flips sign)\n")

cat("\n== high z vs the contrast ================================\n")
cat(sprintf("  cor(z, contrast)   = %+.2f\n", cor(L$z, L$contrast)))
cat(sprintf("  cor(z, |contrast|) = %+.2f   <- a real effect would make this positive\n",
            cor(L$z, abs(L$contrast))))
cat("\n  coherent ROIs with no contrast to show for it:\n")
L %>% filter(z > 3, abs(contrast) < 0.1) %>% arrange(desc(z)) %>% head(6) %>%
  transmute(run = paste0("sub-", subject, " ses-", session, " run-", run), split,
            usable, z = round(z, 2), rbar = round(rbar(z, usable), 3),
            voxel_SD_deg = round(sqrt(-2 * log(rbar(z, usable))) * 180/pi / 6, 1),
            contrast = round(contrast, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")
