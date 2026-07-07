# ── Diagnostic: single-stage fixed-effect slopes (answers Q1/Q2/Q3) ──────────
#
# Q1: c_mono alone, TRANSLATE data      -> is its slope significant?
# Q2: c_nmt  alone, READ data           -> is its slope significant?
# Q3: c_mono + c_nmt jointly, per stage -> what happens to each slope?
#
# Uses lmerTest for p-values (in-sample coefficient significance), which is a
# DIFFERENT question from RQ1's out-of-sample LOO delta-llh.
# Pooled z-scoring so slopes are comparable to the joint brms model.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(lmerTest); library(dplyr)})

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
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
                    c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))

ctrl  <- lmerControl(optimizer="bobyqa")
re    <- "(1|participant) + (1|sentence_id)"
base  <- paste("log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity")

fit <- function(data, add) {
  f <- as.formula(paste(base, "+", add, "+", re))
  suppressMessages(lmer(f, data=data, REML=FALSE, control=ctrl))
}
row <- function(m, term) {
  co <- summary(m)$coefficients
  sprintf("  %-8s beta=%+.4f  SE=%.4f  t=%+.2f  p=%.4f%s",
          term, co[term,"Estimate"], co[term,"Std. Error"],
          co[term,"t value"], co[term,"Pr(>|t|)"],
          ifelse(co[term,"Pr(>|t|)"]<.05,"  *",""))
}

d_tr <- df %>% filter(stage=="translate")
d_rd <- df %>% filter(stage=="read")
cat(sprintf("N translate=%d  N read=%d\n\n", nrow(d_tr), nrow(d_rd)))

cat("── Q1/Q2: SINGLE predictor, per stage ──────────────────────────────\n")
cat("TRANSLATE, c_mono alone:\n"); cat(row(fit(d_tr,"c_mono"),"c_mono"),"\n")
cat("TRANSLATE, c_nmt  alone:\n"); cat(row(fit(d_tr,"c_nmt"), "c_nmt"), "\n")
cat("READ,      c_nmt  alone:\n"); cat(row(fit(d_rd,"c_nmt"), "c_nmt"), "\n")
cat("READ,      c_mono alone:\n"); cat(row(fit(d_rd,"c_mono"),"c_mono"),"\n")

cat("\n── Q3: BOTH predictors jointly, per stage ──────────────────────────\n")
m_tr <- fit(d_tr, "c_nmt + c_mono")
cat("TRANSLATE, c_nmt + c_mono:\n")
cat(row(m_tr,"c_nmt"),"\n"); cat(row(m_tr,"c_mono"),"\n")
m_rd <- fit(d_rd, "c_nmt + c_mono")
cat("READ,      c_nmt + c_mono:\n")
cat(row(m_rd,"c_nmt"),"\n"); cat(row(m_rd,"c_mono"),"\n")
