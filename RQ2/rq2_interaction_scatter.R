# ── RQ2: binned scatterplot visualising the condition x predictor interaction ──
#
# Raw scatterplots at n=11,193 would be an unreadable overplotted blob, so this
# bins each predictor into within-condition deciles and plots mean log(TFD)
# +/- SE per bin. This shows the diverging slopes (the interaction) directly in
# the data, complementing rq2_crossover.pdf (which shows model-based posterior
# slopes, an abstraction) with what the underlying data pattern looks like.
#
# Same data pipeline as rq2_stage_slopes_check.R / rq2_brm_joint.R: pooled
# z-scoring across both conditions, same stoplight exclusion.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT_DIR  <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(dplyr); library(ggplot2); library(lme4)})

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

resid_for <- function(data, drop_var) {
  keep_var <- setdiff(c("c_nmt", "c_mono"), drop_var)
  f <- as.formula(paste0(
    "log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity + condition * ", keep_var,
    " + (1|participant) + (1|sentence_id)"))
  m <- lmer(f, data = data, REML = FALSE, control = ctrl)
  data$partial_resid <- residuals(m) + fixef(m)["(Intercept)"]
  data
}

df_nmt_resid  <- resid_for(df, "c_nmt")
df_mono_resid <- resid_for(df, "c_mono")

bin_summary <- function(data, predictor_col, n_bins = 10) {
  data %>%
    group_by(condition) %>%
    mutate(bin = ntile(.data[[predictor_col]], n_bins)) %>%
    group_by(condition, bin) %>%
    summarise(
      x_mean   = mean(.data[[predictor_col]]),
      y_mean   = mean(partial_resid),
      y_se     = sd(partial_resid) / sqrt(n()),
      .groups = "drop"
    )
}

nmt_bins  <- bin_summary(df_nmt_resid,  "c_nmt")  %>% mutate(predictor = "c[nmt]")
mono_bins <- bin_summary(df_mono_resid, "c_mono") %>% mutate(predictor = "c[mono]")
plot_data <- bind_rows(nmt_bins, mono_bins)

# ── Colours: Okabe-Ito orange/blue, validated colourblind-safe ──────────────
cond_colours <- c("Reading" = "#D55E00", "Translation" = "#0072B2")

p <- ggplot(plot_data, aes(x = x_mean, y = y_mean, colour = condition)) +
  geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se), width = 0, linewidth = 0.5) +
  geom_point(size = 2) +
  geom_smooth(aes(group = condition), method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~predictor, labeller = label_parsed) +
  scale_colour_manual(values = cond_colours, name = NULL) +
  labs(x = "Predictor value (within-condition decile mean, SD units)",
       y = "Partial residual of log(TFD)") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    strip.text = element_text(face = "italic", size = 13),
    legend.position = "top",
    axis.line = element_line(colour = "black", linewidth = 0.3),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  )

ggsave(file.path(OUT_DIR, "rq2_interaction_scatter.pdf"), p, width = 7, height = 4.2, device = "pdf")
cat("Saved rq2_interaction_scatter.pdf\n")
