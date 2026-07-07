# ── RQ1 robustness check: maximal random effects + multiple-comparison correction ──
#
# Two questions:
#   (1) Does c_nmt's LOO-CV predictive gain survive when the predictor model
#       gets a by-participant random slope for the predictor (matching the
#       "maximal" philosophy used for RQ2's joint model), instead of the
#       intercept-only baseline RE used in the main analysis?
#   (2) Do the 7 predictors (c_nmt, c_mono, H_e, f_e, f_eos, f_recv, f_cross)
#       survive Holm correction for multiple comparisons?
#
# RE structure: M0 (baseline, no predictor of interest) keeps intercept-only
# RE, since there is no predictor slope to estimate. M1 (baseline+predictor)
# attempts (1+predictor|participant)+(1|sentence_id); falls back to
# intercept-only on convergence failure/singularity (same rule as
# analysis.R::fit_with_slope). Fallback frequency is logged per predictor.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(lme4); library(dplyr)})

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
attn <- read.csv(file.path(DATA_DIR, "attention_features_6_norm.csv"),  stringsAsFactors=FALSE)
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
  left_join(attn %>% select(sentence_id, word_index,
                             H_e=attn_entropy, f_e=attn_context,
                             f_eos=attn_eos, f_recv=attn_recv, f_cross=attn_cross),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal=surprisal_soft, mono_surprisal,
         log10_freq, H_e, f_e, f_eos, f_recv, f_cross)

df <- fix %>%
  left_join(predictors, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity))

df_trans <- df %>%
  filter(stage=="translate",
         !is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq),
         !is.na(H_e), !is.na(f_e), !is.na(f_eos), !is.na(f_recv), !is.na(f_cross)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df_trans <- df_trans %>%
  mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
         c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq),
         c_He=z(H_e), c_fe=z(f_e), c_feos=z(f_eos), c_frecv=z(f_recv), c_fcross=z(f_cross))

cat(sprintf("N obs = %d  |  N sentences = %d\n\n",
            nrow(df_trans), n_distinct(df_trans$sentence_id)))

ctrl   <- lmerControl(optimizer="bobyqa")
base_f <- log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity +
            (1|participant) + (1|sentence_id)

held_out_llh <- function(model, test_data) {
  mu <- predict(model, newdata=test_data,
                re.form=~(1|participant), allow.new.levels=TRUE)
  mean(dnorm(test_data$log_tfd, mean=mu, sd=sigma(model), log=TRUE))
}
se_fn  <- function(x) sd(x)/sqrt(length(x))
perm_p <- function(x, n_perm=1000, seed=42) {
  obs <- mean(x); set.seed(seed)
  mean(replicate(n_perm, mean(x*sample(c(-1,1),length(x),replace=TRUE))) >= obs)
}

# Predictor model with maximal RE attempt: (1+predictor|participant)+(1|sentence_id),
# falling back to intercept-only on convergence failure / singularity.
fit_predictor_model <- function(pred_var, train) {
  full_slope_f <- as.formula(paste(
    "log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity +", pred_var,
    "+ (1 +", pred_var, "|participant) + (1|sentence_id)"))
  full_int_f <- as.formula(paste(
    "log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity +", pred_var,
    "+ (1|participant) + (1|sentence_id)"))
  m <- tryCatch(suppressMessages(lmer(full_slope_f, data=train, REML=FALSE, control=ctrl)),
                error = function(e) NULL)
  fell_back <- FALSE
  if (is.null(m) || isSingular(m)) {
    m <- suppressMessages(lmer(full_int_f, data=train, REML=FALSE, control=ctrl))
    fell_back <- TRUE
  }
  list(model=m, fell_back=fell_back)
}

pred_vars   <- c("c_nmt","c_mono","c_He","c_fe","c_feos","c_frecv","c_fcross")
pred_labels <- c("c_nmt","c_mono","H_e","f_e","f_eos","f_recv","f_cross")

sids <- unique(df_trans$sentence_id)
n_s  <- length(sids)

delta_mat  <- matrix(0, n_s, length(pred_vars), dimnames=list(NULL, pred_vars))
fallback_count <- setNames(rep(0, length(pred_vars)), pred_vars)

cat("── 200-fold LOO-CV, maximal RE for predictor models ─────────────────────\n")
t0 <- proc.time()
for (i in seq_along(sids)) {
  s     <- sids[i]
  train <- df_trans %>% filter(sentence_id != s)
  test  <- df_trans %>% filter(sentence_id == s)

  m0   <- suppressMessages(lmer(base_f, data=train, REML=FALSE, control=ctrl))
  llh0 <- held_out_llh(m0, test)

  for (pv in pred_vars) {
    fit <- fit_predictor_model(pv, train)
    if (fit$fell_back) fallback_count[pv] <- fallback_count[pv] + 1
    delta_mat[i, pv] <- held_out_llh(fit$model, test) - llh0
  }
  if (i %% 25 == 0 || i == 1)
    cat(sprintf("  [%3d/200]  %.0fs elapsed\n", i, (proc.time()-t0)["elapsed"]))
}
cat(sprintf("\nDone: %.1f min\n", (proc.time()-t0)["elapsed"]/60))

cat("\n── Fallback-to-intercept-only frequency (out of 200 folds) ──────────────\n")
for (k in seq_along(pred_vars))
  cat(sprintf("  %-8s  %d/200 folds fell back\n", pred_labels[k], fallback_count[pred_vars[k]]))

# ── Results + Holm correction ────────────────────────────────────────────────

raw_p <- sapply(pred_vars, function(pv) perm_p(delta_mat[,pv]))
holm_p <- p.adjust(raw_p, method="holm")

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("RQ1 — Maximal-RE LOO-CV with Holm-corrected p-values (200-fold, 1000 perms)\n")
cat("══════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("%-10s %-10s %-10s %-10s %-12s\n","predictor","dllh","SE","p_raw","p_holm"))
for (k in seq_along(pred_vars)) {
  pv <- pred_vars[k]
  d  <- delta_mat[,pv]
  cat(sprintf("%-10s %+.5f   %.5f   %.4f%s     %.4f%s\n",
              pred_labels[k], mean(d), se_fn(d), raw_p[k],
              ifelse(raw_p[k]<.05,"*"," "),
              holm_p[k], ifelse(holm_p[k]<.05,"*"," ")))
}
cat("\nDone.\n")
