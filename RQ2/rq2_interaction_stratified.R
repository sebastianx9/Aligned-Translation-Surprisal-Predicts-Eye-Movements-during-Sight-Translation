# HISTORICAL exploratory script for the superseded four-model ladder. It is
# not used by the authoritative CSF pipeline and expects deprecated M2/M3
# caches. Use RQ2/rq2_interaction_kfold.R for the reported comparison.
# Stratified evaluation of the condition:c_nmt interaction gain.
# The pooled M3-M2 held-out gain dilutes a translation-related stage interaction over
# the mostly-reading sample (see "Attenuation of interaction gains" in Methods).
# No refitting: models are fit on the pooled data; we only re-aggregate the
# already-computed pointwise elpd differences over the translation held-out rows.
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"; OUT <- "/Users/sebastianx/Dissertation RQ2"
suppressMessages({library(dplyr)})

kf <- readRDS(file.path(OUT, "rq2_kfold_elpd.rds"))   # pointwise$M2,$M3 ; sentence_id (df order)
ptw <- kf$pointwise; sid_saved <- kf$sentence_id; N <- kf$N

# ── re-derive df (identical prep to rq2_kfold_elpd.R) to recover cond in row order ──
fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep="\t", header=TRUE,
                   stringsAsFactors=FALSE, quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))
sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1, .groups="drop")
pred <- nmt %>% left_join(sl, by="sentence_id") %>%
  mutate(word_length=nchar(word), word_position=word_index/(sent_len-1), word_lower=tolower(trimws(word))) %>%
  left_join(freq, by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id, word_index, mono_surprisal=surprisal_sum), by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position, nmt_surprisal=surprisal_soft, mono_surprisal, log10_freq)
df <- fix %>% filter(stage %in% c("translate","read")) %>%
  left_join(pred, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity), cond=as.integer(stage=="translate")) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L), by=c("sentence_id","word_index"))

stopifnot(nrow(df) == N, all(df$sentence_id == sid_saved))   # row-order alignment check
cat(sprintf("Alignment OK: N=%d (translate %d, read %d)\n", N, sum(df$cond==1), sum(df$cond==0)))

d_i <- ptw$M3 - ptw$M2   # condition:c_nmt interaction, pointwise held-out elpd gain
sid <- df$sentence_id; cond <- df$cond

signflip <- function(d_s, n=10000){o<-sum(d_s); set.seed(42)
  mean(replicate(n, sum(d_s*sample(c(-1,1),length(d_s),replace=TRUE)))>=o)}

eval_on <- function(mask, tag){
  di <- d_i[mask]; ss <- sid[mask]; ds <- tapply(di, ss, sum); S <- length(ds)
  data.frame(subset=tag, n=sum(mask), n_sent=S, elpd_diff=sum(di),
             se_cluster=sd(ds)*sqrt(S), z=sum(di)/(sd(ds)*sqrt(S)), p=signflip(ds))
}
res <- rbind(
  eval_on(rep(TRUE,N),  "pooled (all obs)"),
  eval_on(cond==1,      "translation only"),
  eval_on(cond==0,      "reading only")
)
res$p_note <- ifelse(res$p<.05,"*","")
cat("\n=== condition:c_nmt interaction gain (M3 - M2), by evaluation subset ===\n")
print(format(res, digits=3, nsmall=3))
saveRDS(res, file.path(OUT, "rq2_interaction_stratified.rds"))
cat("\nSaved rq2_interaction_stratified.rds\nDONE\n")
