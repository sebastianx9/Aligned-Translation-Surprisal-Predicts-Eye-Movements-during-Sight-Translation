# Pre-modelling screening of the c_nmt distribution (artefact-exclusion figure).
# Exact main-pipeline sample BEFORE the stoplight exclusion (translate stage,
# NA-filtered on mono/freq, n = 5,152 observations): z-scored at observation
# level, one point per unique word position, in corpus order. Reproduces
# stoplight = +11.81 SD. Copy-scan positions (tab:copy) as open circles.
OUT <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
RQ1 <- "/Users/sebastianx/Dissertation RQ1"
suppressMessages({library(dplyr); library(ggplot2)})

fix  <- read.csv(file.path(DATA_DIR,"fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR,"monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="") %>%
  select(word=Word,log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))
pred <- nmt %>% mutate(word_lower=tolower(trimws(word))) %>%
  left_join(freq,by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id,word_index,mono_surprisal=surprisal_sum),
            by=c("sentence_id","word_index"))
df <- fix %>% filter(stage=="translate") %>%
  left_join(pred %>% select(sentence_id,word_index,surprisal_soft,mono_surprisal,log10_freq),
            by=c("sentence_id","word_index")) %>%
  filter(!is.na(surprisal_soft),!is.na(mono_surprisal),!is.na(log10_freq))
df$c_nmt <- (df$surprisal_soft - mean(df$surprisal_soft)) / sd(df$surprisal_soft)

wp <- df %>% distinct(sentence_id, word_index, c_nmt) %>%
  arrange(sentence_id, word_index) %>% mutate(idx = row_number())
sl_val <- wp %>% filter(sentence_id=="S003", word_index==3)
nx <- max(wp$c_nmt[!(wp$sentence_id=="S003" & wp$word_index==3)])
cat(sprintf("obs=%d positions=%d | stoplight=%+.2f SD | next largest=%+.2f SD\n",
            nrow(df), nrow(wp), sl_val$c_nmt, nx))

copy <- read.csv(file.path(RQ1,"copy_failures_full_scan.csv"), stringsAsFactors=FALSE) %>%
  left_join(wp, by=c("sentence_id","word_index")) %>%
  filter(!is.na(c_nmt)) %>%
  mutate(status = ifelse(sentence_id=="S003" & word_index==3, "excluded", "retained"))
sl <- filter(copy, status=="excluded")

p <- ggplot(wp, aes(idx, c_nmt)) +
  geom_point(colour="#0072B2", alpha=0.35, size=0.9) +
  geom_point(data=copy, colour="#0072B2", size=1.8, shape=1, stroke=0.8) +
  annotate("text", x=sl$idx+60, y=sl$c_nmt, colour="grey30", size=3.6, fontface="italic",
           label=sprintf("stoplight (%+.2f SD)", sl$c_nmt), hjust=0) +
  labs(x="Source-word positions (corpus order, S001 to S200)",
       y=expression(italic(c)[nmt]~"(SD units, before exclusion)")) +
  theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(),
        axis.line.y=element_line(colour="black", linewidth=0.3))
ggsave(file.path(OUT,"rq1_cnmt_screening.pdf"), p, width=6.8, height=3.4, device="pdf")
cat("saved rq1_cnmt_screening.pdf\n")
