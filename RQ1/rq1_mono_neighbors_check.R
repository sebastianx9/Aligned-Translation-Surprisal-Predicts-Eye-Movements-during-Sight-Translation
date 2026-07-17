# Robustness check: does current-word c_mono survive neighbour-surprisal controls?
# Motivated by the c_mono alignment bug: the shifted (successor) values predicted
# reading fixations numerically better than the corrected values, hinting that
# neighbouring words' difficulty may matter. Here spillover (i-1) and successor
# (i+1) monolingual surprisal enter alongside the current word's, both stages.
# lmerTest (Satterthwaite) for speed; interior words only (both neighbours exist).
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(dplyr); library(lmerTest)})

fix  <- read.csv(file.path(DATA_DIR,"fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR,"monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="") %>%
  select(word=Word,log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))

m <- mono %>% select(sentence_id, word_index, mono_curr = surprisal_sum) %>%
  group_by(sentence_id) %>% arrange(word_index, .by_group=TRUE) %>%
  mutate(mono_prev = lag(mono_curr), mono_next = lead(mono_curr)) %>% ungroup()

sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1,.groups="drop")
pred <- nmt %>% left_join(sl,by="sentence_id") %>%
  mutate(word_length=nchar(word), word_position=word_index/(sent_len-1),
         word_lower=tolower(trimws(word))) %>%
  left_join(freq,by=c("word_lower"="word")) %>%
  left_join(m, by=c("sentence_id","word_index")) %>%
  select(sentence_id,word_index,word_length,word_position,log10_freq,
         nmt_surprisal=surprisal_soft,mono_curr,mono_prev,mono_next)

base <- fix %>% filter(stage %in% c("translate","read")) %>%
  left_join(pred,by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal),!is.na(mono_curr),!is.na(mono_prev),
         !is.na(mono_next),!is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))

z <- function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

for (st in c("read","translate")) {
  df <- base %>% filter(stage==st) %>%
    mutate(c_mono=z(mono_curr), c_prev=z(mono_prev), c_next=z(mono_next),
           c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))
  cat(sprintf("\n══════════ %s stage (interior words, n=%d) ══════════\n", st, nrow(df)))
  m1 <- lmer(log_tfd ~ c_wlen+c_wpos+c_freq+ambiguity + c_mono +
             (1|participant)+(1|sentence_id), data=df, REML=FALSE)
  m2 <- lmer(log_tfd ~ c_wlen+c_wpos+c_freq+ambiguity + c_prev + c_mono + c_next +
             (1|participant)+(1|sentence_id), data=df, REML=FALSE)
  cat("-- current-word only --\n")
  print(round(summary(m1)$coefficients["c_mono",,drop=FALSE],4))
  cat("-- with neighbours --\n")
  print(round(summary(m2)$coefficients[c("c_prev","c_mono","c_next"),],4))
  cat(sprintf("LRT m2 vs m1: Chisq=%.2f, p=%.4f\n",
      anova(m1,m2)$Chisq[2], anova(m1,m2)$`Pr(>Chisq)`[2]))
}
cat("\nDONE\n")
