# ── RQ2: LOO-CV with Lim et al. normalised attention features ─────────────────
#
# Features from attention_features_6_norm.csv (divided by dummy attention):
#   H_e  / log(N_ctx)       — entropy relative to uniform
#   f_e  / (N_ctx/N)        — context flow relative to uniform
#   f_eos/ (1/N)            — eos flow relative to uniform  (= raw × N)
#   f_recv/(N_ctx/N)        — received flow relative to uniform
#   f_cross/(N_tgt/N)       — cross-attn flow relative to uniform
#
# Same 2×5 comparison as rq2_diagnostic.R:
#   vs baseline alone  |  vs baseline + c_nmt
# ──────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"

library(lme4)
library(dplyr)

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
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

z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
df_trans <- df_trans %>%
  mutate(c_nmt    = z(nmt_surprisal), c_wlen=z(word_length),
         c_wpos   = z(word_position),  c_freq=z(log10_freq),
         c_He     = z(H_e),   c_fe    = z(f_e),
         c_feos   = z(f_eos), c_frecv = z(f_recv), c_fcross=z(f_cross))

cat(sprintf("N obs = %d  |  N sentences = %d\n\n",
            nrow(df_trans), n_distinct(df_trans$sentence_id)))

# Feature distributions after normalisation
cat("── Normalised feature distributions ────────────────────────────────────\n")
for (f in c("H_e","f_e","f_eos","f_recv","f_cross")) {
  x <- df_trans[[f]]
  cat(sprintf("  %-8s  mean=%.3f  sd=%.3f  min=%.3f  max=%.3f\n",
              f, mean(x), sd(x), min(x), max(x)))
}

# Pearson r with log_tfd
cat("\n── Pearson r with log_tfd ───────────────────────────────────────────────\n")
for (f in c("c_nmt","c_He","c_fe","c_feos","c_frecv","c_fcross")) {
  r <- cor(df_trans[[f]], df_trans$log_tfd, use="complete.obs")
  cat(sprintf("  %-9s  r = %+.4f\n", f, r))
}

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
sids <- unique(df_trans$sentence_id)
n_s  <- length(sids)

delta_vs_base <- matrix(0, n_s, length(feat_vars), dimnames=list(NULL,feat_vars))
delta_vs_nmt  <- matrix(0, n_s, length(feat_vars), dimnames=list(NULL,feat_vars))

cat("\n── LOO-CV (200-fold, individual features × 2 baselines) ────────────────\n")
t0 <- proc.time()
for (i in seq_along(sids)) {
  s     <- sids[i]
  train <- df_trans %>% filter(sentence_id != s)
  test  <- df_trans %>% filter(sentence_id == s)

  m0  <- suppressMessages(lmer(base_f, data=train, REML=FALSE, control=ctrl))
  m_n <- suppressMessages(lmer(nmt_f,  data=train, REML=FALSE, control=ctrl))
  llh0 <- held_out_llh(m0,  test)
  llhn <- held_out_llh(m_n, test)

  for (fv in feat_vars) {
    f1 <- update(base_f, as.formula(paste(". ~ . +", fv)))
    f2 <- update(nmt_f,  as.formula(paste(". ~ . +", fv)))
    m1 <- suppressMessages(lmer(f1, data=train, REML=FALSE, control=ctrl))
    m2 <- suppressMessages(lmer(f2, data=train, REML=FALSE, control=ctrl))
    delta_vs_base[i, fv] <- held_out_llh(m1, test) - llh0
    delta_vs_nmt[i,  fv] <- held_out_llh(m2, test) - llhn
  }
  if (i %% 50 == 0 || i == 1)
    cat(sprintf("  [%3d/200]  %.0fs\n", i, (proc.time()-t0)["elapsed"]))
}
cat(sprintf("\nDone: %.1f min\n", (proc.time()-t0)["elapsed"]/60))

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("RQ2 — Normalised features — 200-fold LOO-CV (1000 perms, seed=42)\n")
cat(sprintf("%-10s  %-38s  %-38s\n","",
            "── vs baseline ────────────────────────",
            "── vs baseline+c_nmt ──────────────────"))
cat(sprintf("%-10s  %-8s %-8s %-6s    %-8s %-8s %-6s\n",
            "feature","Δllh","SE","p","Δllh","SE","p"))
cat(rep("─",82),"\n",sep="")

for (k in seq_along(feat_vars)) {
  fv <- feat_vars[k]
  d1 <- delta_vs_base[,fv]; d2 <- delta_vs_nmt[,fv]
  p1 <- perm_p(d1); p2 <- perm_p(d2)
  cat(sprintf("%-10s  %+.5f  %.5f  %.3f%s   %+.5f  %.5f  %.3f%s\n",
              feat_labels[k],
              mean(d1), se_fn(d1), p1, ifelse(p1<.05,"*"," "),
              mean(d2), se_fn(d2), p2, ifelse(p2<.05,"*"," ")))
}
cat(rep("─",82),"\n",sep="")
cat(sprintf("%-10s  %+.5f  %.5f  %.3f%s\n","c_nmt(ref)",0.00297,0.00118,0.001,"*"))
cat("\n── Raw (unnormalised) reference ─────────────────────────────────────────\n")
cat("vs baseline:   H_e p=0.659  f_e p=0.416  f_eos p=0.432  f_recv p=0.769  f_cross p=0.593\n")
cat("vs base+c_nmt: H_e p=0.752  f_e p=0.290  f_eos p=0.281  f_recv p=0.685  f_cross p=0.520\n")
