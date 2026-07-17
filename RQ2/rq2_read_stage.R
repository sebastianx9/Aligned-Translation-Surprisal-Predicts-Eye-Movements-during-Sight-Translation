# ── Attention features predicting TFD in reading-aloud stage ──────────────────
#
# Compare:
#   (A) attention features (normalised) vs baseline  — reading-aloud stage
#   (B) c_nmt vs baseline                            — reading-aloud stage (reference)
#
# If attention features predict reading-aloud but not translation:
#   → they capture source reading difficulty, overridden by translation load
# If neither:
#   → simply underpowered / features not useful in this dataset
# If c_nmt also predicts reading-aloud:
#   → c_nmt captures general difficulty, not translation-specific
# ──────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"

library(lme4)
library(dplyr)

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
attn <- read.csv(file.path(DATA_DIR, "attention_features_6_norm.csv"),  stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep="\t",
                   header=TRUE, stringsAsFactors=FALSE, quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>%
  mutate(word=tolower(trimws(word)))

sent_lengths <- nmt %>%
  group_by(sentence_id) %>%
  summarise(sent_len=max(word_index)+1, .groups="drop")

predictors <- nmt %>%
  left_join(sent_lengths, by="sentence_id") %>%
  mutate(word_length   = nchar(word),
         word_position = word_index / (sent_len - 1),
         word_lower    = tolower(trimws(word))) %>%
  left_join(freq, by=c("word_lower"="word")) %>%
  left_join(attn %>% select(sentence_id, word_index,
                             H_e=attn_entropy, f_e=attn_context,
                             f_eos=attn_eos, f_recv=attn_recv, f_cross=attn_cross),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal=surprisal_soft, log10_freq, H_e, f_e, f_eos, f_recv, f_cross)

df <- fix %>%
  left_join(predictors, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity))

# ── Reading-aloud stage ───────────────────────────────────────────────────────

df_read <- df %>%
  filter(stage=="read",
         !is.na(nmt_surprisal), !is.na(log10_freq),
         !is.na(H_e), !is.na(f_e), !is.na(f_eos), !is.na(f_recv), !is.na(f_cross)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
df_read <- df_read %>%
  mutate(c_nmt    = z(nmt_surprisal), c_wlen=z(word_length),
         c_wpos   = z(word_position),  c_freq=z(log10_freq),
         c_He     = z(H_e),   c_fe    = z(f_e),
         c_feos   = z(f_eos), c_frecv = z(f_recv), c_fcross=z(f_cross))

cat(sprintf("Reading-aloud: N obs = %d  |  N sentences = %d\n\n",
            nrow(df_read), n_distinct(df_read$sentence_id)))

# ── LOO-CV ────────────────────────────────────────────────────────────────────

ctrl <- lmerControl(optimizer="bobyqa")
base_f <- log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity +
            (1|participant) + (1|sentence_id)
nmt_f  <- update(base_f, .~.+c_nmt)

held_out_llh <- function(model, test_data) {
  mu <- predict(model, newdata=test_data,
                re.form=~(1|participant), allow.new.levels=TRUE)
  mean(dnorm(test_data$log_tfd, mean=mu, sd=sigma(model), log=TRUE))
}
se_fn  <- function(x) sd(x)/sqrt(length(x))
perm_p <- function(x, n_perm=1000) {
  obs <- mean(x); set.seed(42)
  mean(replicate(n_perm, mean(x*sample(c(-1,1),length(x),replace=TRUE))) >= obs)
}

feat_vars   <- c("c_He","c_fe","c_feos","c_frecv","c_fcross")
feat_labels <- c("H_e","f_e","f_eos","f_recv","f_cross")
sids <- unique(df_read$sentence_id)
n_s  <- length(sids)

delta_nmt   <- numeric(n_s)
delta_feats <- matrix(0, n_s, length(feat_vars), dimnames=list(NULL,feat_vars))

cat("── LOO-CV (reading-aloud stage) ─────────────────────────────────────────\n")
t0 <- proc.time()
for (i in seq_along(sids)) {
  s     <- sids[i]
  train <- df_read %>% filter(sentence_id != s)
  test  <- df_read %>% filter(sentence_id == s)

  m0  <- suppressMessages(lmer(base_f, data=train, REML=FALSE, control=ctrl))
  m_n <- suppressMessages(lmer(nmt_f,  data=train, REML=FALSE, control=ctrl))
  llh0         <- held_out_llh(m0,  test)
  delta_nmt[i] <- held_out_llh(m_n, test) - llh0

  for (fv in feat_vars) {
    f1 <- update(base_f, as.formula(paste(".~.+", fv)))
    m1 <- suppressMessages(lmer(f1, data=train, REML=FALSE, control=ctrl))
    delta_feats[i, fv] <- held_out_llh(m1, test) - llh0
  }
  if (i %% 50 == 0 || i == 1)
    cat(sprintf("  [%3d/200]  %.0fs\n", i, (proc.time()-t0)["elapsed"]))
}
cat(sprintf("\nDone: %.1f min\n", (proc.time()-t0)["elapsed"]/60))

# ── Results ───────────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("Reading-aloud stage — 200-fold LOO-CV (normalised features, 1000 perms)\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")

p_nmt <- perm_p(delta_nmt)
cat(sprintf("c_nmt    Δllh=%+.5f  SE=%.5f  p=%.3f%s\n",
            mean(delta_nmt), se_fn(delta_nmt), p_nmt,
            ifelse(p_nmt<.05,"  *","")))

cat("\nAttention features vs baseline:\n")
cat(sprintf("  %-8s  %-9s %-8s %-6s\n","feature","Δllh","SE","p"))
cat(rep("─",40),"\n",sep="")
for (k in seq_along(feat_vars)) {
  fv <- feat_vars[k]
  d  <- delta_feats[,fv]
  p  <- perm_p(d)
  cat(sprintf("  %-8s  %+.5f  %.5f  %.3f%s\n",
              feat_labels[k], mean(d), se_fn(d), p,
              ifelse(p<.05,"  *","")))
}

cat("\n── Translation stage reference ──────────────────────────────────────────\n")
cat("c_nmt    Δllh=+0.00297  p=0.001 *\n")
cat("H_e      Δllh=-0.00025  p=1.000\n")
cat("f_e      Δllh=-0.00006  p=0.573\n")
cat("f_eos    Δllh=-0.00031  p=1.000\n")
cat("f_recv   Δllh=-0.00031  p=0.919\n")
cat("f_cross  Δllh=-0.00014  p=0.705\n")
