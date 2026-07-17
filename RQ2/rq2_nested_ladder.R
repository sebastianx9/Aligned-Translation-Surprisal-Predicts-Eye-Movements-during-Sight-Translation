# ── RQ2: nested out-of-sample ladder (terrain map, frequentist) ───────────────
#
# Answers RQ2 inside the same Delta llh framework as RQ1/RQ3, avoiding the
# cross-sample non-commensurability that makes a naive translate-vs-read
# Delta llh comparison invalid: every comparison below is a NESTED pair
# scored on the SAME held-out observations.
#
#   M1 = controls + cond + c_mono + cond:c_mono
#   M2 = M1 + c_nmt                  -> Delta llh(M2|M1) : c_nmt beyond c_mono
#   M3 = M2 + cond:c_nmt             -> Delta llh(M3|M2) : c_nmt's effect is task-modulated
#   N1 = controls + cond + c_nmt + cond:c_nmt
#                                     -> Delta llh(M3|N1) : does c_mono add beyond c_nmt?
#   (note N2 == M3, so four fits per fold suffice)
#
# Convention (identical to rq1_loo_authoritative.R; see Dissertation_PROVENANCE.md):
#   * RE identical across every compared pair: (1|participant) + (1|sentence_id)
#   * sigma = each model's own training residual SD
#   * predict re.form = ~(1|participant); held-out sentence's intercept unused
#   * 200-fold leave-one-sentence-out; sign-flip permutation test, 1000 perms, seed 42
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(dplyr); library(lme4)})

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep="\t",
                   header=TRUE, stringsAsFactors=FALSE, quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))

sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1, .groups="drop")
pred <- nmt %>% left_join(sl, by="sentence_id") %>%
  mutate(word_length=nchar(word), word_position=word_index/(sent_len-1),
         word_lower=tolower(trimws(word))) %>%
  left_join(freq, by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id, word_index, mono_surprisal=surprisal_sum),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal=surprisal_soft, mono_surprisal, log10_freq)

df <- fix %>% filter(stage %in% c("translate","read")) %>%
  left_join(pred, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity),
         cond = as.integer(stage == "translate")) %>%       # 1 = sight translation, 0 = reading aloud
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L), by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)   # pooled across both stages
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
                     c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))

cat(sprintf("Pooled: N = %d (translate %d, read %d) | sentences = %d | participants = %d\n\n",
            nrow(df), sum(df$cond==1), sum(df$cond==0),
            n_distinct(df$sentence_id), n_distinct(df$participant)))

ctrl <- lmerControl(optimizer="bobyqa")
RE   <- "(1|participant) + (1|sentence_id)"
CTRL <- "c_wlen + c_wpos + c_freq + ambiguity"

f_M1 <- as.formula(paste("log_tfd ~", CTRL, "+ cond + c_mono + cond:c_mono +", RE))
f_M2 <- as.formula(paste("log_tfd ~", CTRL, "+ cond + c_mono + cond:c_mono + c_nmt +", RE))
f_M3 <- as.formula(paste("log_tfd ~", CTRL, "+ cond + c_mono + cond:c_mono + c_nmt + cond:c_nmt +", RE))
f_N1 <- as.formula(paste("log_tfd ~", CTRL, "+ cond + c_nmt + cond:c_nmt +", RE))

held_out_llh <- function(model, test) {
  mu <- predict(model, newdata=test, re.form=~(1|participant), allow.new.levels=TRUE)
  mean(dnorm(test$log_tfd, mean=mu, sd=sigma(model), log=TRUE))
}
se_fn  <- function(x) sd(x)/sqrt(length(x))
perm_p <- function(x, n_perm=1000) {
  obs <- mean(x); set.seed(42)
  mean(replicate(n_perm, mean(x*sample(c(-1,1), length(x), replace=TRUE))) >= obs)
}

sids <- sort(unique(df$sentence_id)); n_s <- length(sids)
d_beyond <- d_modul <- d_reverse <- numeric(n_s)
sing <- 0
t0 <- proc.time()

for (i in seq_along(sids)) {
  train <- df %>% filter(sentence_id != sids[i])
  test  <- df %>% filter(sentence_id == sids[i])
  m1 <- suppressMessages(lmer(f_M1, train, REML=FALSE, control=ctrl))
  m2 <- suppressMessages(lmer(f_M2, train, REML=FALSE, control=ctrl))
  m3 <- suppressMessages(lmer(f_M3, train, REML=FALSE, control=ctrl))
  n1 <- suppressMessages(lmer(f_N1, train, REML=FALSE, control=ctrl))
  sing <- sing + sum(sapply(list(m1,m2,m3,n1), isSingular))

  l1 <- held_out_llh(m1, test); l2 <- held_out_llh(m2, test)
  l3 <- held_out_llh(m3, test); ln <- held_out_llh(n1, test)

  d_beyond[i]  <- l2 - l1     # c_nmt beyond c_mono
  d_modul[i]   <- l3 - l2     # c_nmt's task modulation
  d_reverse[i] <- l3 - ln     # c_mono (+interaction) beyond c_nmt (+interaction)

  if (i %% 25 == 0) cat(sprintf("  [%3d/%d]  %.0fs\n", i, n_s, (proc.time()-t0)["elapsed"]))
}

report <- function(nm, d, question) {
  cat(sprintf("\n%s\n  %s\n  dllh = %+.5f   SE = %.5f   p_perm = %.4f%s\n",
              nm, question, mean(d), se_fn(d), perm_p(d),
              ifelse(perm_p(d) < .05, "  *", "")))
}
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("RQ2 NESTED LADDER — pooled, 200-fold LOO-CV, 1000 perms\n")
cat(sprintf("RE (identical across every pair): %s\n", RE))
cat(sprintf("singular fits: %d / %d\n", sing, 4*n_s))
cat("══════════════════════════════════════════════════════════════════\n")
report("Delta llh(M2 | M1)", d_beyond,  "Does c_nmt add predictive power beyond c_mono?")
report("Delta llh(M3 | M2)", d_modul,   "Is c_nmt's effect task-modulated?")
report("Delta llh(M3 | N1)", d_reverse, "Does c_mono add anything beyond c_nmt?")

saveRDS(list(beyond=d_beyond, modul=d_modul, reverse=d_reverse, n=nrow(df)),
        "/Users/sebastianx/Dissertation RQ2/rq2_nested_ladder.rds")
cat("\nSaved rq2_nested_ladder.rds\n")
