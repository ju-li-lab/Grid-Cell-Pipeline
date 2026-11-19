# ============================================================================
#  plot_spatial_cognition.R — every subject's spatial-cognition scores over
#  all sessions
# ----------------------------------------------------------------------------
#  Needs spatial_cog from import_spatial_cognition.R. One plot per test x
#  measure, one panel per condition (plus their mean, "all"), one line per
#  subject across every session in the file. Filled points are the sessions
#  that stand for an MRI session. The plots go to the Plots pane and, one
#  page each, to results/spatial_cognition_<study>.pdf.
#  Plots are wrapped in print() because source() drops bare ggplots.
# ============================================================================

library(dplyr)
library(stringr)
library(ggplot2)

STUDY    <- "s17"   # NULL = first imported study, or e.g. "s06"
MEASURES <- NULL   # NULL = every measure, or e.g. c("Score", "RT")

stopifnot(exists("spatial_cog"))
if (is.null(STUDY)) STUDY <- unique(spatial_cog$study)[1]

sc <- spatial_cog %>%
  filter(study == STUDY, is.null(MEASURES) | measure %in% MEASURES) %>%
  group_by(subject_id, sc_session, session, test, condition, measure) %>%
  summarise(value = mean(value), .groups = "drop")
stopifnot(nrow(sc) > 0)

# A test with several conditions also gets their mean, as the analysis uses it
condition_means <- sc %>%
  group_by(test) %>%
  filter(n_distinct(condition) > 1) %>%
  group_by(subject_id, sc_session, session, test, measure) %>%
  summarise(condition = "all", value = mean(value), .groups = "drop")
sc <- bind_rows(sc, condition_means) %>%
  mutate(test      = factor(test, levels = unique(spatial_cog$test)),
         condition = factor(condition,
                            levels = c(setdiff(unique(spatial_cog$condition), "all"), "all")),
         mri       = !is.na(session),
         subject   = str_extract(subject_id, "(?<=^sub-)\\d+"))

spaghetti <- function(d, title) {
  line_ends <- d %>%
    group_by(condition, subject_id) %>%
    filter(sc_session == max(sc_session))
  ggplot(d, aes(x = sc_session, y = value, group = subject_id)) +
    geom_line(colour = "#2a78d6", alpha = 0.55, linewidth = 0.6) +
    geom_point(aes(shape = mri), colour = "#2a78d6", size = 2) +
    geom_text(data = line_ends, aes(label = subject), hjust = -0.4, size = 2.6,
              colour = "#52514e", check_overlap = TRUE) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                       labels = c(`FALSE` = "no", `TRUE` = "yes"),
                       name = "stands for an MRI session") +
    scale_x_continuous(breaks = sort(unique(d$sc_session)),
                       expand = expansion(mult = c(0.03, 0.1))) +
    facet_wrap(~ condition) +
    labs(title = title, x = "session, as numbered in the file", y = NULL) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

plots <- sc %>%
  group_by(test, measure) %>%
  group_map(~ spaghetti(.x, paste0(STUDY, " ", .y$test, ": ", .y$measure)))

for (p in plots) print(p)

dir.create("results", showWarnings = FALSE)
pdf_file <- file.path("results", paste0("spatial_cognition_", STUDY, ".pdf"))
pdf(pdf_file, width = 11, height = 7)
for (p in plots) print(p)
invisible(dev.off())
message(length(plots), " plots, also saved to ", normalizePath(pdf_file))

rm(condition_means, p)
