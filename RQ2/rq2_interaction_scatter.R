# ── RQ2: full-data scatterplot visualising the condition x predictor interaction ──
#
# Every word-participant observation is plotted (n=11,193), not binned means,
# using heavy point transparency to handle overplotting. The y-axis is a
# partial residual (controls, condition, and the other predictor's effect
# removed), not raw log(TFD): the raw/unadjusted relationship is dominated by
# word length/frequency confounds shared across both conditions, which masks
# the differential, condition-specific effect the joint model isolates (see
# rq2_eda_pairplot.R for the raw/unadjusted relationships instead).
#
# Same data pipeline as rq2_stage_slopes_check.R / rq2_brm_joint.R: pooled
# z-scoring across both conditions, same stoplight exclusion.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT_DIR  <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(dplyr); library(ggplot2); library(lme4); library(patchwork)})

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep="\t",
                   header=TRUE, stringsAsFactors=FALSE, quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))

sent_lengths <- nmt %>% group_by(sentence_id) %>%
  summarise(sent_len=max(word_index)+1, .groups="drop")

predictors <- nmt %>% left_join(sent_lengths, by="sentence_id") %>%
  mutate(word_length=nchar(word), word_position=word_index/(sent_len-1),
         word_lower=tolower(trimws(word))) %>%
  left_join(freq, by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id, word_index, mono_surprisal=surprisal_sum),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal=surprisal_soft, mono_surprisal, log10_freq)

df <- fix %>% filter(stage %in% c("translate","read")) %>%
  left_join(predictors, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity),
         condition=factor(ifelse(stage=="translate","Translation","Reading"),
                           levels=c("Reading","Translation"))) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
                     c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))

cat(sprintf("N = %d (translate: %d, read: %d)\n",
            nrow(df), sum(df$condition=="Translation"), sum(df$condition=="Reading")))

# ── Partial residuals: remove controls + condition + the OTHER predictor's
# effect before binning, so the plotted relationship isolates each predictor's
# own effect (matching what the joint model's coefficients report) rather than
# the raw marginal relationship, which is dominated by word length/frequency
# confounds shared across both conditions and masks the differential effect. ──

ctrl <- lmerControl(optimizer = "bobyqa")

# Correct component-plus-residual construction: adding back only the
# intercept (as an earlier version of this script did) leaves the residual
# exactly orthogonal to keep_var by construction (keep_var is already in
# the model via condition*keep_var), so a naive intercept-only "partial
# residual" has a within-condition slope against keep_var of ~0 (verified
# numerically), not the model's coefficient. Adding back keep_var's own
# fitted contribution (its main effect + its interaction with condition)
# is what reproduces the model's actual per-condition slope in the plot.
resid_for <- function(data, drop_var) {
  keep_var <- setdiff(c("c_nmt", "c_mono"), drop_var)
  f <- as.formula(paste0(
    "log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity + condition * ", keep_var,
    " + (1|participant) + (1|sentence_id)"))
  m <- lmer(f, data = data, REML = FALSE, control = ctrl)
  mm <- model.matrix(m)
  b  <- fixef(m)
  keep_cols <- grep(keep_var, colnames(mm), fixed = TRUE, value = TRUE)
  component <- as.numeric(mm[, keep_cols, drop = FALSE] %*% b[keep_cols])
  data$partial_resid <- residuals(m) + component + b["(Intercept)"]
  data
}

df_nmt_resid  <- resid_for(df, "c_nmt")  %>% transmute(condition, x = c_nmt,  y = partial_resid, predictor = "c[nmt]")
df_mono_resid <- resid_for(df, "c_mono") %>% transmute(condition, x = c_mono, y = partial_resid, predictor = "c[mono]")
plot_data <- bind_rows(df_nmt_resid, df_mono_resid)

# ── Binned summary: fixed-width bins (not equal-count deciles), so bin
# width reflects the actual data scale rather than being squeezed by the
# right-skewed distribution. Point size encodes n per bin, giving an honest
# view of where the data mass is (big bubbles near the centre, small ones
# in the sparse tails) instead of showing 11k raw overplotted points. Bins
# with n<3 are dropped as too noisy to represent a mean. ───────────────────
binwidth <- 0.25
bin_summary <- plot_data %>%
  mutate(x_bin = round(x / binwidth) * binwidth) %>%
  group_by(predictor, condition, x_bin) %>%
  summarise(y_mean = mean(y), n = n(), .groups = "drop") %>%
  filter(n >= 3)

# ── Colours: Okabe-Ito orange/blue, validated colourblind-safe ──────────────
cond_colours <- c("Reading" = "#D55E00", "Translation" = "#0072B2")

common_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    strip.text = element_text(face = "italic", size = 13),
    legend.position = "top",
    axis.line = element_line(colour = "black", linewidth = 0.3),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  )

# ── Top row: every raw observation (n=11,193), heavy transparency to handle
# overplotting. Shows the true noise level directly, at the cost of making
# the trend hard to see by eye. ─────────────────────────────────────────────
p_raw <- ggplot(plot_data, aes(x = x, y = y, colour = condition)) +
  geom_point(size = 0.5, alpha = 0.12) +
  geom_smooth(aes(group = condition, fill = condition), method = "lm",
              se = TRUE, linewidth = 0.8, alpha = 0.15) +
  facet_wrap(~predictor, labeller = label_parsed) +
  scale_colour_manual(values = cond_colours, name = NULL) +
  scale_fill_manual(values = cond_colours, guide = "none") +
  labs(x = "Predictor value (SD units)", y = "Partial residual of log(TFD)") +
  common_theme +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(legend.position = "top")

# ── Bottom row: fixed-width bin means, bubble size = n per bin. Solid
# filled circles (shape 19, coloured not filled) rather than shape 21 with
# colour=NA -- ggplot2 silently drops every row of a layer when colour is
# set to NA as a fixed parameter, which is what produced empty bubbles in
# an earlier version of this figure. Soft alpha lets overlapping bubbles
# blend instead of looking cluttered. ───────────────────────────────────
p_bin <- ggplot(plot_data, aes(x = x, y = y, colour = condition)) +
  geom_smooth(aes(group = condition, fill = condition), method = "lm",
              se = TRUE, linewidth = 0.8, alpha = 0.15) +
  geom_point(data = bin_summary,
             aes(x = x_bin, y = y_mean, size = n),
             shape = 19, alpha = 0.36) +
  facet_wrap(~predictor, labeller = label_parsed) +
  scale_colour_manual(values = cond_colours, name = NULL) +
  scale_fill_manual(values = cond_colours, guide = "none") +
  scale_size_area(max_size = 9, name = "N per bin") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Predictor value (SD units)", y = "Partial residual of log(TFD)",
       title = "Binned means (bubble size = N per bin)") +
  common_theme +
  theme(plot.title = element_text(size = 11, face = "plain", colour = "grey40"))

ggsave(file.path(OUT_DIR, "rq2_interaction_scatter.pdf"), p_raw, width = 7, height = 4.2, device = "pdf")
cat("Saved rq2_interaction_scatter.pdf\n")
