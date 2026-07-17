# ── RQ2: descriptive pairplot of raw variables ────────────────────────────────
#
# Pure exploratory/descriptive figure: raw (unadjusted) c_nmt, c_mono, and
# log(TFD), with univariate histograms on the diagonal and pairwise
# relationships off-diagonal, split by condition. No modelling here.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT_DIR  <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(dplyr); library(ggplot2); library(GGally)})

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
  mutate(log_tfd=log(total_fixation_duration_ms),
         condition=factor(ifelse(stage=="translate","Translation","Reading"),
                           levels=c("Reading","Translation"))) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

# Standardised (z-scored), matching the c_nmt / c_mono notation used
# throughout the rest of the dissertation, so this figure's axes are the
# same variables the models are fitted on, not a different raw quantity.
z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
plot_df <- df %>%
  transmute(condition,
            `c[nmt]`  = z(nmt_surprisal),
            `c[mono]` = z(mono_surprisal),
            `log(TFD)` = log_tfd)

cond_colours <- c("Reading" = "#D55E00", "Translation" = "#0072B2")

p <- ggpairs(
  plot_df, columns = 2:4, aes(colour = condition),
  upper = list(continuous = wrap("cor", size = 3.5)),
  lower = list(continuous = wrap("points", size = 0.4, alpha = 0.15)),
  diag  = list(continuous = function(data, mapping, ...) {
    ggplot(data, mapping) +
      geom_density(aes(colour = condition, fill = condition), alpha = 0.15, linewidth = 0.8) +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank())
  }),
  labeller = label_parsed,
  legend = c(2, 1)
) +
  scale_colour_manual(values = cond_colours, name = NULL) +
  scale_fill_manual(values = cond_colours, name = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "rq2_eda_pairplot.pdf"), p, width = 7.5, height = 7, device = "pdf")
cat("Saved rq2_eda_pairplot.pdf\n")
