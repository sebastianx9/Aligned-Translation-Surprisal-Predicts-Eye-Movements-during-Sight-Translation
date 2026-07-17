# ── RQ1: fixed-effect estimates for the full sight-translation model ──────────
# Source of Table 2 in main.tex (the in-sample coefficient table).
# Primary RE structure: (1|participant) + (1|sentence_id), matching the LOO-CV.
# lmerTest / Satterthwaite. Translate stage, stoplight excluded, n = 5,149.

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(dplyr); library(lmerTest)})

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

df <- fix %>% filter(stage=="translate") %>%
  left_join(pred, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L), by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_wlen=z(word_length),
                     c_wpos=z(word_position), c_freq=z(log10_freq))

cat(sprintf("n = %d | sentences = %d | participants = %d\n\n",
            nrow(df), n_distinct(df$sentence_id), n_distinct(df$participant)))

m <- lmer(log_tfd ~ c_nmt + c_wlen + c_wpos + c_freq + ambiguity +
            (1|participant) + (1|sentence_id), data=df, REML=TRUE)

cat("── Table 2: fixed-effect estimates (primary RE, Satterthwaite) ──\n")
print(round(summary(m)$coefficients, 4))
cat(sprintf("\nc_nmt on exponentiated scale: %.1f%% longer TFD per SD\n",
            100*(exp(fixef(m)["c_nmt"]) - 1)))
