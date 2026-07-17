# ── Reading-aloud stage LOO-CV for c_mono (fills the gap in rq2_read_stage.R) ──
# Same pipeline / baseline / RE structure / seed as rq2_read_stage.R so the
# resulting Delta llh is directly comparable to the other reading-stage values.

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(dplyr); library(lme4)})

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
attn <- read.csv(file.path(DATA_DIR, "attention_features_6.csv"),       stringsAsFactors=FALSE)
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
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal=surprisal_soft, mono_surprisal, log10_freq,
         H_e, f_e, f_eos, f_recv, f_cross)

df <- fix %>% left_join(predictors, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity))

# Same reading-stage filter as rq2_read_stage.R, plus mono non-NA
df_read <- df %>%
  filter(stage=="read",
         !is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq),
         !is.na(H_e), !is.na(f_e), !is.na(f_eos), !is.na(f_recv), !is.na(f_cross)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L), by=c("sentence_id","word_index"))

z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
df_read <- df_read %>% mutate(c_mono = z(mono_surprisal), c_wlen=z(word_length),
                              c_wpos = z(word_position),  c_freq=z(log10_freq))

cat(sprintf("Reading-aloud: N obs = %d  |  N sentences = %d\n",
            nrow(df_read), n_distinct(df_read$sentence_id)))

ctrl   <- lmerControl(optimizer="bobyqa")
base_f <- log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity + (1|participant) + (1|sentence_id)
mono_f <- update(base_f, .~.+c_mono)

held_out_llh <- function(model, test_data) {
  mu <- predict(model, newdata=test_data, re.form=~(1|participant), allow.new.levels=TRUE)
  mean(dnorm(test_data$log_tfd, mean=mu, sd=sigma(model), log=TRUE))
}
se_fn  <- function(x) sd(x)/sqrt(length(x))
perm_p <- function(x, n_perm=1000) {
  obs <- mean(x); set.seed(42)
  mean(replicate(n_perm, mean(x*sample(c(-1,1),length(x),replace=TRUE))) >= obs)
}

sids <- unique(df_read$sentence_id); n_s <- length(sids)
d_mono <- numeric(n_s)
for (i in seq_along(sids)) {
  s <- sids[i]
  train <- df_read %>% filter(sentence_id != s)
  test  <- df_read %>% filter(sentence_id == s)
  m_base <- suppressMessages(lmer(base_f, data=train, REML=FALSE, control=ctrl))
  m_mono <- suppressMessages(lmer(mono_f, data=train, REML=FALSE, control=ctrl))
  d_mono[i] <- held_out_llh(m_mono, test) - held_out_llh(m_base, test)
  if (i %% 50 == 0) cat(sprintf("  [%3d/%d]\n", i, n_s))
}

cat("\n══════ Reading-aloud stage, c_mono vs controls-only baseline ══════\n")
cat(sprintf("c_mono   Delta llh = %+.5f  SE = %.5f  p_perm = %.3f\n",
            mean(d_mono), se_fn(d_mono), perm_p(d_mono)))
