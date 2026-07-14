# RQ2 joint model, FULLY MAXIMAL symmetric random effects (item 5b).
# RE: (1 + condition + c_nmt + c_mono + condition:c_nmt + condition:c_mono | participant)
#   + (1 + condition + c_nmt + c_mono + condition:c_nmt + condition:c_mono | sentence_id)
# Big model: high adapt_delta. Reports fixef, stage-specific slopes, diagnostics,
# and the old-vs-new interaction coefficients for comparison. Falls back is manual.
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"; OUT <- "/Users/sebastianx/Dissertation RQ2"
suppressMessages({library(brms); library(dplyr)}); options(mc.cores=4)
fix  <- read.csv(file.path(DATA_DIR,"fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR,"monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="") %>%
  select(word=Word,log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))
sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1,.groups="drop")
pred <- nmt %>% left_join(sl,by="sentence_id") %>%
  mutate(word_length=nchar(word),word_position=word_index/(sent_len-1),word_lower=tolower(trimws(word))) %>%
  left_join(freq,by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id,word_index,mono_surprisal=surprisal_sum),by=c("sentence_id","word_index")) %>%
  select(sentence_id,word_index,word_length,word_position,nmt_surprisal=surprisal_soft,mono_surprisal,log10_freq)
df <- fix %>% filter(stage %in% c("translate","read")) %>%
  left_join(pred,by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms),ambiguity=factor(ambiguity),condition=as.integer(stage=="translate")) %>%
  filter(!is.na(nmt_surprisal),!is.na(mono_surprisal),!is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df<-df%>%mutate(c_nmt=z(nmt_surprisal),c_mono=z(mono_surprisal),c_wlen=z(word_length),c_wpos=z(word_position),c_freq=z(log10_freq))
cat(sprintf("Pooled N=%d (translate %d, read %d)\n", nrow(df), sum(df$condition==1), sum(df$condition==0)))
priors <- c(prior(normal(0,1),class=b),prior(normal(6,1),class=Intercept),
            prior(exponential(1),class=sd),prior(exponential(1),class=sigma),prior(lkj(2),class=cor))
re <- "(1 + condition + c_nmt + c_mono + condition:c_nmt + condition:c_mono | participant) + (1 + condition + c_nmt + c_mono + condition:c_nmt + condition:c_mono | sentence_id)"
form <- as.formula(paste("log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity + condition*c_nmt + condition*c_mono +", re))
t0<-proc.time()
m<-brm(form,data=df,prior=priors,control=list(adapt_delta=0.99,max_treedepth=15),
       chains=4,iter=4000,warmup=2000,save_pars=save_pars(all=TRUE),silent=2,refresh=0)
saveRDS(m,file.path(DATA_DIR,"brm_cache","rq2_joint_maximal.rds"))
cat(sprintf("\nfit time %.0f min\n",(proc.time()-t0)["elapsed"]/60))
cat("\n=== fixef (MAXIMAL RE) ===\n"); print(round(fixef(m),4))
h<-hypothesis(m,c(nmt_read="c_nmt=0", nmt_trans="c_nmt + condition:c_nmt = 0",
                  mono_read="c_mono=0", mono_trans="c_mono + condition:c_mono = 0"))
cat("\n=== stage-specific slopes ===\n"); print(h)
nd<-sum(subset(nuts_params(m),Parameter=="divergent__")$Value)
cat(sprintf("\nRhat max: %.4f | divergences: %d\n", round(max(rhat(m),na.rm=TRUE),4), nd))
cat("\n=== by-sentence & by-participant RE SDs (degeneracy check) ===\n"); print(VarCorr(m))
cat("DONE\n")
