# RQ1 in-sample coefficients (Table 2), brms version — translate stage full model.
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(brms); library(dplyr)})
options(mc.cores=4)
fix  <- read.csv(file.path(DATA_DIR,"fixation_durations_word.csv"),stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))
sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1,.groups="drop")
pred <- nmt %>% left_join(sl,by="sentence_id") %>%
  mutate(word_length=nchar(word),word_position=word_index/(sent_len-1),word_lower=tolower(trimws(word))) %>%
  left_join(freq,by=c("word_lower"="word")) %>%
  select(sentence_id,word_index,word_length,word_position,nmt_surprisal=surprisal_soft,log10_freq)
df <- fix %>% filter(stage=="translate") %>% left_join(pred,by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms),ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal),!is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df<-df%>%mutate(c_nmt=z(nmt_surprisal),c_wlen=z(word_length),c_wpos=z(word_position),c_freq=z(log10_freq))
cat(sprintf("translate n=%d\n",nrow(df)))
pri<-c(prior(normal(0,1),class=b),prior(normal(6,1),class=Intercept),prior(exponential(1),class=sd),prior(exponential(1),class=sigma))
m<-brm(log_tfd~c_wlen+c_wpos+c_freq+ambiguity+c_nmt+(1|participant)+(1|sentence_id),
       data=df,prior=pri,chains=4,iter=2000,warmup=1000,silent=2,refresh=0)
saveRDS(m,file.path(DATA_DIR,"brm_cache","rq1_coef_brms.rds"))
print(round(fixef(m),4))
cat("\nRhat max:",round(max(rhat(m),na.rm=TRUE),4)," done\n")
