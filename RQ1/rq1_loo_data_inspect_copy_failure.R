# ── RQ1: Sentence-level LOO-CV with ambiguity covariate ───────────────────────
#
# 200-fold leave-one-sentence-out CV. No random seed for fold assignment.
# Baseline includes `ambiguity` (A/U) to reduce residual variance and produce
# cleaner per-sentence LPD estimates.
#
# SE = sd(delta_s) / sqrt(200)
# ──────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT_DIR  <- "/Users/sebastianx/Dissertation RQ1"

library(lme4)
library(dplyr)

# ── 1. Load data ──────────────────────────────────────────────────────────────

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

# Copy-failure exclusion: corpus scan (copy_failures_full_scan.csv) found 26
# verbatim-copy positions across 200 sentences. Of these, 15 are content words
# (proper nouns excluded). 14 involve established Czech loanwords whose c_nmt
# falls within the normal range (<=1.32 SD); only "stoplight" (S003) has no
# Czech borrowed form and produces artefactually inflated c_nmt (+11.81 SD).
# Only "stoplight" is excluded; other copy-failure words have valid c_nmt.
copy_failures <- data.frame(
  sentence_id = "S003",
  word_index  = 3L,
  stringsAsFactors = FALSE
)

df_trans <- df_trans %>%
  anti_join(copy_failures, by = c("sentence_id", "word_index"))

z <- function(x, ref) (x - mean(ref, na.rm = TRUE)) / sd(ref, na.rm = TRUE)
df_trans <- df_trans %>%
  mutate(c_nmt  = z(nmt_surprisal,  nmt_surprisal),
         c_mono = z(mono_surprisal, mono_surprisal),
         c_wlen = z(word_length,    word_length),
         c_wpos = z(word_position,  word_position),
         c_freq = z(log10_freq,     log10_freq))

cat(sprintf("Excluded: 'stoplight' (S003) — 3 word-participant observations\n"))
cat(sprintf("Translate-stage observations after exclusion: %d\n", nrow(df_trans)))
cat(sprintf("Sentences: %d  (A: %d, U: %d)\n",
            n_distinct(df_trans$sentence_id),
            sum(df_trans$ambiguity == "A" & !duplicated(df_trans$sentence_id)),
            sum(df_trans$ambiguity == "U" & !duplicated(df_trans$sentence_id))))

# ── 2. Model formulae ─────────────────────────────────────────────────────────

ctrl   <- lmerControl(optimizer = "bobyqa")
base_f <- log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity +
  (1 | participant) + (1 | sentence_id)
nmt_f  <- update(base_f, . ~ . + c_nmt)
mono_f <- update(base_f, . ~ . + c_mono)

held_out_llh <- function(model, test_data) {
  mu  <- predict(model, newdata = test_data,
                 re.form = ~(1 | participant), allow.new.levels = TRUE)
  sig <- sigma(model)
  mean(dnorm(test_data$log_tfd, mean = mu, sd = sig, log = TRUE))
}

# ── 3. 200-fold LOO-CV ────────────────────────────────────────────────────────

sids  <- unique(df_trans$sentence_id)
n_s   <- length(sids)

delta_nmt  <- numeric(n_s)
delta_mono <- numeric(n_s)

cat(sprintf("\nRunning %d-fold sentence LOO-CV (with ambiguity covariate)...\n", n_s))
t_start <- proc.time()

for (i in seq_along(sids)) {
  if (i %% 50 == 0 || i == 1) {
    elapsed <- (proc.time() - t_start)["elapsed"]
    eta     <- if (i > 1) elapsed / (i - 1) * (n_s - i + 1) else NA
    cat(sprintf("  [%3d/%d]  elapsed: %.0fs  ETA: %.0fs\n", i, n_s, elapsed, eta))
  }

  s     <- sids[i]
  train <- df_trans %>% filter(sentence_id != s)
  test  <- df_trans %>% filter(sentence_id == s)

  m0   <- suppressMessages(lmer(base_f,  data = train, REML = FALSE, control = ctrl))
  m_nt <- suppressMessages(lmer(nmt_f,   data = train, REML = FALSE, control = ctrl))
  m_mn <- suppressMessages(lmer(mono_f,  data = train, REML = FALSE, control = ctrl))

  llh0          <- held_out_llh(m0,   test)
  delta_nmt[i]  <- held_out_llh(m_nt, test) - llh0
  delta_mono[i] <- held_out_llh(m_mn, test) - llh0
}

cat(sprintf("\nDone. Total time: %.1f min\n", (proc.time() - t_start)["elapsed"] / 60))

# ── 4. Inference ──────────────────────────────────────────────────────────────

se_fn <- function(x) sd(x) / sqrt(length(x))

perm_p_fn <- function(x, n_perm = 1000) {
  obs <- mean(x)
  mean(replicate(n_perm,
    mean(x * sample(c(-1, 1), length(x), replace = TRUE))) >= obs)
}

set.seed(42)
p_nmt  <- perm_p_fn(delta_nmt)
p_mono <- perm_p_fn(delta_mono)

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("RQ1: 200-fold sentence LOO-CV (ambiguity covariate in baseline)\n")
cat("SE = sd(delta_s) / sqrt(200)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

for (nm in c("c_nmt", "c_mono")) {
  d  <- if (nm == "c_nmt") delta_nmt else delta_mono
  p  <- if (nm == "c_nmt") p_nmt     else p_mono
  ci <- mean(d) + c(-1, 1) * qt(0.975, df = n_s - 1) * se_fn(d)
  cat(sprintf("%-8s  Δllh = %+.5f  SE = %.5f  95%%CI = [%+.5f, %+.5f]  p_perm = %.4f%s\n",
              nm, mean(d), se_fn(d), ci[1], ci[2], p,
              ifelse(p < .01, " **", ifelse(p < .05, " *", ""))))
}

# ── 5. Save results ───────────────────────────────────────────────────────────

results_df <- data.frame(
  sentence_id = sids,
  delta_nmt   = delta_nmt,
  delta_mono  = delta_mono
)
write.csv(results_df, file.path(OUT_DIR, "rq1_loo_ambiguity_results.csv"), row.names = FALSE)
cat("\nSaved: rq1_loo_ambiguity_results.csv\n")
