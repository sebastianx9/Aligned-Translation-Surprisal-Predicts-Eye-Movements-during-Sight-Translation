# ── RQ3: descriptive pairplot of raw variables ─────────────────────────────
#
# Pure exploratory/descriptive figure, analogous to rq2_eda_pairplot.R:
# c_nmt, c_mono, log(GD), and log(RRT), diagonal densities + lower-triangle
# scatter + upper-triangle correlations. No modelling here. Translate-stage
# only (RQ3 does not compare conditions), so there is no colour-split
# grouping variable, unlike the RQ2 version.
#
# log(GD) is defined for all fixated translate-stage words (n=5,240);
# log(RRT) is defined only for words with at least one regression-in
# (n=2,536) and is NA otherwise. GGally computes each pairwise panel on
# whatever rows have both variables defined, so GD-involving panels use
# the full n and RRT-involving panels use the smaller, conditional-on-
# regression subset -- matching how the two outcomes are actually
# analysed separately in RQ3 (see rq_twostage_loo.R).
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT_DIR  <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(dplyr); library(ggplot2); library(GGally)})

em   <- read.csv(file.path(DATA_DIR, "eye_measures_word.csv"),          stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)

predictors <- nmt %>%
  left_join(mono %>% select(sentence_id, word_index, mono_surprisal=surprisal_sum),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, nmt_surprisal=surprisal_soft, mono_surprisal)

df <- em %>% filter(stage == "translate") %>%
  left_join(predictors, by=c("sentence_id","word_index")) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
plot_df <- df %>%
  transmute(
    `c[nmt]`    = z(nmt_surprisal),
    `c[mono]`   = z(mono_surprisal),
    `log(GD)`   = ifelse(!is.na(gd_ms) & gd_ms > 0, log(gd_ms), NA_real_),
    `log(RRT)`  = ifelse(!is.na(regress_in) & regress_in == 1 & !is.na(rrt_ms) & rrt_ms > 0,
                          log(rrt_ms), NA_real_)
  )

point_colour <- "#0072B2"

p <- ggpairs(
  plot_df, columns = 1:4,
  upper = list(continuous = wrap("cor", size = 3.5, colour = "grey30")),
  lower = list(continuous = wrap("points", size = 0.4, alpha = 0.15, colour = point_colour)),
  diag  = list(continuous = function(data, mapping, ...) {
    ggplot(data, mapping) +
      geom_density(colour = point_colour, fill = point_colour, alpha = 0.2, linewidth = 0.8) +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank())
  }),
  labeller = label_parsed
) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10))

ggsave(file.path(OUT_DIR, "rq3_eda_pairplot.pdf"), p, width = 7.5, height = 7, device = "pdf")
cat("Saved rq3_eda_pairplot.pdf\n")
