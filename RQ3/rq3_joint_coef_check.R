# ── RQ3: in-sample joint-model coefficient check (matches RQ2's non-redundancy
# ── operationalisation: does the coefficient survive, not does predictive
# ── power improve). Full-data lmerTest fit, no LOO-CV.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(lmerTest); library(dplyr)})

eye  <- read.csv(file.path(DATA_DIR, "eye_measures_word.csv"),          stringsAsFactors=FALSE)
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

df <- eye %>%
  filter(stage=="translate") %>%
  left_join(predictors, by=c("sentence_id","word_index")) %>%
  mutate(ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>%
  mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
         c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))

df_gd  <- df %>% filter(gd_ms  > 0) %>% mutate(log_y = log(gd_ms))
df_rrt <- df %>% filter(regress_in == 1, rrt_ms > 0) %>% mutate(log_y = log(rrt_ms))

cat(sprintf("N GD = %d  |  N RRT = %d\n\n", nrow(df_gd), nrow(df_rrt)))

ctrl <- lmerControl(optimizer="bobyqa")
f_solo_nmt  <- log_y ~ c_wlen + c_wpos + c_freq + ambiguity + c_nmt  + (1|participant) + (1|sentence_id)
f_solo_mono <- log_y ~ c_wlen + c_wpos + c_freq + ambiguity + c_mono + (1|participant) + (1|sentence_id)
f_joint     <- log_y ~ c_wlen + c_wpos + c_freq + ambiguity + c_nmt + c_mono + (1|participant) + (1|sentence_id)

row <- function(m, term) {
  co <- summary(m)$coefficients
  sprintf("  %-8s beta=%+.4f  SE=%.4f  t=%+.2f  p=%.4f%s",
          term, co[term,"Estimate"], co[term,"Std. Error"],
          co[term,"t value"], co[term,"Pr(>|t|)"],
          ifelse(co[term,"Pr(>|t|)"]<.05,"  *",""))
}

cat("══════════════════════════════════════════════════════════════════════════\n")
cat("RQ3 joint-model coefficient check (in-sample, matches RQ2's non-redundancy test)\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")

cat("── GD ──\n")
m_solo <- lmer(f_solo_nmt, data=df_gd, REML=FALSE, control=ctrl)
cat("c_nmt alone:\n"); cat(row(m_solo, "c_nmt"), "\n")
m_joint_gd <- lmer(f_joint, data=df_gd, REML=FALSE, control=ctrl)
cat("c_nmt + c_mono jointly:\n")
cat(row(m_joint_gd, "c_nmt"), "\n")
cat(row(m_joint_gd, "c_mono"), "\n\n")

cat("── RRT ──\n")
m_solo_r <- lmer(f_solo_nmt, data=df_rrt, REML=FALSE, control=ctrl)
cat("c_nmt alone:\n"); cat(row(m_solo_r, "c_nmt"), "\n")
m_joint_rrt <- lmer(f_joint, data=df_rrt, REML=FALSE, control=ctrl)
cat("c_nmt + c_mono jointly:\n")
cat(row(m_joint_rrt, "c_nmt"), "\n")
cat(row(m_joint_rrt, "c_mono"), "\n")

cat("\nDone.\n")
