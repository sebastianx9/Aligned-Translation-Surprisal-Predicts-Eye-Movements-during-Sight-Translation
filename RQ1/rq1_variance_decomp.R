# ── Variance decomposition of word-level LPD differences ──────────────────────
#
# Purpose: diagnose why Var(delta_s) is large in sentence-level LOO-CV.
# Is it primarily between-sentence heterogeneity (real signal differences),
# or within-sentence sampling noise (too few words per sentence)?
#
# Method: fit full-data M0 and M_nmt, compute word-level LPD_diff, then:
#   Model 1: lmer(LPD_diff ~ 1            + (1|sentence_id))
#   Model 2: lmer(LPD_diff ~ ambiguity    + (1|sentence_id))
#
# Key quantities:
#   sigma2_between  = VarCorr(...)$sentence_id  (between-sentence variance)
#   sigma2_within   = sigma(...)^2              (within-sentence variance)
#   sigma2_within / n_words_per_sentence        (sampling noise component)
#
# If sampling noise << sigma2_between:
#   → high Var(delta_s) is genuine sentence heterogeneity, not measurement noise
#   → adding more words per sentence would NOT help
#   → need to find sentence-level variables that explain between-sentence variance
#
# Note: in-sample LPD (not held-out) — for diagnostic purposes only.
# ──────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"

library(lme4)
library(dplyr)

# ── 1. Load and prepare data (same pipeline as main analysis) ─────────────────

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors = FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors = FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors = FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep = "\t",
                   header = TRUE, stringsAsFactors = FALSE, quote = "") %>%
  select(word = Word, log10_freq = Lg10WF) %>%
  mutate(word = tolower(trimws(word)))

sent_lengths <- nmt %>%
  group_by(sentence_id) %>%
  summarise(sent_len = max(word_index) + 1, .groups = "drop")

predictors <- nmt %>%
  left_join(sent_lengths, by = "sentence_id") %>%
  mutate(word_length   = nchar(word),
         word_position = word_index / (sent_len - 1),
         word_lower    = tolower(trimws(word))) %>%
  left_join(freq, by = c("word_lower" = "word")) %>%
  left_join(mono %>% select(sentence_id, word_index,
                             mono_surprisal = surprisal_sum),
            by = c("sentence_id", "word_index")) %>%
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal = surprisal_soft, mono_surprisal, log10_freq)

df <- fix %>%
  left_join(predictors, by = c("sentence_id", "word_index")) %>%
  mutate(log_tfd   = log(total_fixation_duration_ms),
         ambiguity = factor(ambiguity))

df_trans <- df %>%
  filter(stage == "translate",
         !is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq))

z <- function(x, ref) (x - mean(ref, na.rm = TRUE)) / sd(ref, na.rm = TRUE)
df_trans <- df_trans %>%
  mutate(c_nmt  = z(nmt_surprisal,  nmt_surprisal),
         c_mono = z(mono_surprisal, mono_surprisal),
         c_wlen = z(word_length,    word_length),
         c_wpos = z(word_position,  word_position),
         c_freq = z(log10_freq,     log10_freq))

n_words_per_sent <- nrow(df_trans) / n_distinct(df_trans$sentence_id)
cat(sprintf("N obs: %d   N sentences: %d   Mean words/sentence: %.1f\n",
            nrow(df_trans), n_distinct(df_trans$sentence_id), n_words_per_sent))

# ── 2. Fit full-data M0 and M_nmt (in-sample) ────────────────────────────────

ctrl   <- lmerControl(optimizer = "bobyqa")
base_f <- log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity +
  (1 | participant) + (1 | sentence_id)
nmt_f  <- update(base_f, . ~ . + c_nmt)

cat("\nFitting full-data models...\n")
m0   <- lmer(base_f, data = df_trans, REML = FALSE, control = ctrl)
m_nt <- lmer(nmt_f,  data = df_trans, REML = FALSE, control = ctrl)

# Word-level LPD difference (in-sample)
lpd0  <- dnorm(df_trans$log_tfd, fitted(m0),   sigma(m0),   log = TRUE)
lpd_nt <- dnorm(df_trans$log_tfd, fitted(m_nt), sigma(m_nt), log = TRUE)
df_trans$lpd_diff <- lpd_nt - lpd0

cat(sprintf("Mean word-level LPD_diff: %+.6f\n", mean(df_trans$lpd_diff)))

# ── 3. Variance decomposition ─────────────────────────────────────────────────

cat("\nFitting variance decomposition models...\n")
vc1 <- lmer(lpd_diff ~ 1         + (1 | sentence_id), data = df_trans, REML = TRUE)
vc2 <- lmer(lpd_diff ~ ambiguity + (1 | sentence_id), data = df_trans, REML = TRUE)

extract_vc <- function(model, label) {
  sb2 <- as.numeric(VarCorr(model)$sentence_id)   # sigma2_between
  sw2 <- sigma(model)^2                            # sigma2_within
  noise <- sw2 / n_words_per_sent                  # sampling noise per sentence
  total <- sb2 + noise
  list(label = label, sb2 = sb2, sw2 = sw2,
       noise = noise, total = total,
       pct_noise = 100 * noise / total,
       icc = sb2 / (sb2 + sw2))
}

r1 <- extract_vc(vc1, "No covariate")
r2 <- extract_vc(vc2, "+ ambiguity")

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("Variance decomposition of word-level LPD_diff  (c_nmt vs baseline)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

for (r in list(r1, r2)) {
  cat(sprintf("[ %s ]\n", r$label))
  cat(sprintf("  σ²_between (sentence)       = %10.7f\n", r$sb2))
  cat(sprintf("  σ²_within  (word-level)     = %10.7f\n", r$sw2))
  cat(sprintf("  σ²_within / %.0f (noise/sent) = %10.7f\n", n_words_per_sent, r$noise))
  cat(sprintf("  Var(δ_s) ≈ σ²_between + noise = %10.7f\n", r$total))
  cat(sprintf("  Sampling noise %% of Var(δ_s) = %6.1f%%\n", r$pct_noise))
  cat(sprintf("  ICC                          = %10.4f\n\n", r$icc))
}

cat(sprintf("Ambiguity explains %.1f%% of σ²_between\n",
            100 * (r1$sb2 - r2$sb2) / r1$sb2))

cat("\n── Interpretation ────────────────────────────────────────────────────\n")
if (r1$pct_noise < 20) {
  cat("Sampling noise is a MINOR component of Var(δ_s).\n")
  cat("→ More words per sentence would not meaningfully reduce SE.\n")
  cat("→ Between-sentence heterogeneity is the dominant source of variance.\n")
  cat("→ Look for sentence-level variables (length, syntactic complexity,\n")
  cat("  ambiguity, congruency) to explain why some sentences show\n")
  cat("  stronger c_nmt effects than others.\n")
} else {
  cat("Sampling noise is a SUBSTANTIAL component of Var(δ_s).\n")
  cat("→ More words per sentence (or longer sentences) would help.\n")
}
