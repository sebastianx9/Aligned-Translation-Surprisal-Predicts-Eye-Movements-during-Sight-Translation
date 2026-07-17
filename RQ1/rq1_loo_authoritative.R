# ── RQ1 AUTHORITATIVE 200-fold LOO-CV (translate stage, all seven predictors) ──
#
# Replaces the previously assembled Table 1, whose SE column had been taken
# from the reading-stage run. Everything here is produced in ONE run on ONE
# sample, under ONE convention:
#
#   * Random effects IDENTICAL in baseline and target: (1|participant) + (1|sentence_id).
#     This is required for Delta llh to isolate the target predictor's fixed effect.
#     A by-participant random slope for the target predictor cannot be held constant
#     (the baseline does not contain that predictor), so intercept-only is the
#     maximal structure consistent with a nested comparison. Lim et al. (2024) do
#     the same: their baseline and target models share the same random effects.
#   * sigma: each model's own training residual SD (sigma(model)) -- Lim's convention.
#   * Prediction on held-out sentence: fixed effects + the participant random
#     intercept learned in training; the held-out sentence's intercept is NOT used.
#   * 200-fold leave-one-sentence-out; paired sign-flip permutation test, 1000 perms, seed 42.
#   * Holm correction across the seven tests reported alongside raw p.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(dplyr); library(lme4)})

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
# Lim et al. (2024) normalise each attentional feature by its "dummy" value
# (the value it would take under uniform attention). Use the normalised file.
attn <- read.csv(file.path(DATA_DIR, "attention_features_6_norm.csv"),  stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep="\t",
                   header=TRUE, stringsAsFactors=FALSE, quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))

sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1, .groups="drop")

predictors <- nmt %>% left_join(sl, by="sentence_id") %>%
  mutate(word_length=nchar(word), word_position=word_index/(sent_len-1),
         word_lower=tolower(trimws(word))) %>%
  left_join(freq, by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id, word_index, mono_surprisal=surprisal_sum),
            by=c("sentence_id","word_index")) %>%
  left_join(attn %>% select(sentence_id, word_index, H_e=attn_entropy, f_e=attn_context,
                            f_eos=attn_eos, f_recv=attn_recv, f_cross=attn_cross),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position, log10_freq,
         nmt_surprisal=surprisal_soft, mono_surprisal, H_e, f_e, f_eos, f_recv, f_cross)

df <- fix %>% filter(stage == "translate") %>%
  left_join(predictors, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq),
         !is.na(H_e), !is.na(f_e), !is.na(f_eos), !is.na(f_recv), !is.na(f_cross)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L), by=c("sentence_id","word_index"))

z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
                     c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq),
                     c_He=z(H_e), c_fe=z(f_e), c_feos=z(f_eos),
                     c_frecv=z(f_recv), c_fcross=z(f_cross))

cat(sprintf("Translate stage: N obs = %d | N sentences = %d | N participants = %d\n\n",
            nrow(df), n_distinct(df$sentence_id), n_distinct(df$participant)))

ctrl <- lmerControl(optimizer="bobyqa")

held_out_llh <- function(model, test) {
  mu <- predict(model, newdata=test, re.form=~(1|participant), allow.new.levels=TRUE)
  mean(dnorm(test$log_tfd, mean=mu, sd=sigma(model), log=TRUE))   # each model's own sigma
}
se_fn  <- function(x) sd(x)/sqrt(length(x))
perm_p <- function(x, n_perm=1000) {
  obs <- mean(x); set.seed(42)
  mean(replicate(n_perm, mean(x*sample(c(-1,1), length(x), replace=TRUE))) >= obs)
}

base_f <- log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity + (1|participant) + (1|sentence_id)

pred_vars   <- c("c_nmt","c_mono","c_He","c_fe","c_feos","c_frecv","c_fcross")
pred_labels <- c("c_nmt","c_mono","H_e","f_e","f_eos","f_recv","f_cross")

sids <- sort(unique(df$sentence_id)); n_s <- length(sids)
deltas <- matrix(NA_real_, nrow=n_s, ncol=length(pred_vars),
                 dimnames=list(NULL, pred_vars))

t0 <- proc.time()
for (i in seq_along(sids)) {
  train <- df %>% filter(sentence_id != sids[i])
  test  <- df %>% filter(sentence_id == sids[i])
  m_base <- suppressMessages(lmer(base_f, data=train, REML=FALSE, control=ctrl))
  llh_b  <- held_out_llh(m_base, test)
  for (v in pred_vars) {
    f_t <- update(base_f, as.formula(paste(". ~ . +", v)))
    m_t <- suppressMessages(lmer(f_t, data=train, REML=FALSE, control=ctrl))
    deltas[i, v] <- held_out_llh(m_t, test) - llh_b
  }
  if (i %% 25 == 0) cat(sprintf("  [%3d/%d]  %.0fs\n", i, n_s, (proc.time()-t0)["elapsed"]))
}

dllh <- colMeans(deltas)
ses  <- apply(deltas, 2, se_fn)
praw <- apply(deltas, 2, perm_p)
pholm <- p.adjust(praw, method="holm")

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("RQ1 AUTHORITATIVE — translate stage, 200-fold LOO-CV, 1000 perms\n")
cat("RE identical in baseline & target: (1|participant)+(1|sentence_id)\n")
cat("sigma = each model's own; predict re.form = ~(1|participant)\n")
cat("══════════════════════════════════════════════════════════════════\n")
cat(sprintf("%-10s %-11s %-10s %-9s %-9s\n","predictor","dllh","SE","p_raw","p_holm"))
for (k in seq_along(pred_vars)) {
  cat(sprintf("%-10s %+.5f    %.5f    %.4f%s   %.4f%s\n",
              pred_labels[k], dllh[k], ses[k],
              praw[k], ifelse(praw[k] < .05, "*", " "),
              pholm[k], ifelse(pholm[k] < .05, "*", " ")))
}
saveRDS(list(deltas=deltas, dllh=dllh, se=ses, p_raw=praw, p_holm=pholm,
             n=nrow(df), labels=pred_labels),
        "/Users/sebastianx/Dissertation RQ1/rq1_loo_authoritative.rds")
cat("\nSaved rq1_loo_authoritative.rds\n")
