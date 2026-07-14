# RQ3 brms kfold elpd — c_nmt / c_mono on GD and RRT (translate stage).
# Bayesian replacement for the lmer two-stage LOO-CV. Same convention as
# rq1_kfold_elpd.R: 10-fold sentence-grouped elpd, sentence-clustered SE,
# sentence-level sign-flip permutation. Fixes the two issues in the old lmer
# script (shared sigma(m_base); it also used the un-normalised attention file,
# irrelevant here since only c_nmt/c_mono are reported).
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"; OUT <- "/Users/sebastianx/Dissertation RQ3"
suppressMessages({library(brms); library(dplyr)}); options(mc.cores = 4)

em   <- read.csv(file.path(DATA_DIR,"eye_measures_word.csv"),          stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR,"monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="") %>%
  select(word=Word,log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))
sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1,.groups="drop")
pred <- nmt %>% left_join(sl,by="sentence_id") %>%
  mutate(word_length=nchar(word),word_position=word_index/(sent_len-1),word_lower=tolower(trimws(word))) %>%
  left_join(freq,by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id,word_index,mono_surprisal=surprisal_sum),by=c("sentence_id","word_index")) %>%
  select(sentence_id,word_index,word_length,word_position,log10_freq,nmt_surprisal=surprisal_soft,mono_surprisal)
base_df <- em %>% filter(stage=="translate") %>% left_join(pred,by=c("sentence_id","word_index")) %>%
  mutate(ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal),!is.na(mono_surprisal),!is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
base_df <- base_df %>% mutate(c_nmt=z(nmt_surprisal),c_mono=z(mono_surprisal),
                              c_wlen=z(word_length),c_wpos=z(word_position),c_freq=z(log10_freq))
df_gd  <- base_df %>% filter(!is.na(gd_ms),  gd_ms  > 0) %>% mutate(log_gd  = log(gd_ms))
df_rrt <- base_df %>% filter(regress_in==1, rrt_ms > 0)  %>% mutate(log_rrt = log(rrt_ms))
cat(sprintf("GD n=%d  RRT n=%d\n", nrow(df_gd), nrow(df_rrt)))

pri<-c(prior(normal(0,1),class=b),prior(normal(6,1),class=Intercept),prior(exponential(1),class=sd),prior(exponential(1),class=sigma))
CACHE<-file.path(DATA_DIR,"brm_cache"); CTRL<-"c_wlen+c_wpos+c_freq+ambiguity"; RE<-"(1|participant)+(1|sentence_id)"
sfp<-function(ds,n=10000){o<-sum(ds);set.seed(42);mean(replicate(n,sum(ds*sample(c(-1,1),length(ds),replace=TRUE)))>=o)}

run_outcome <- function(df, y, tag) {
  df <- df[!is.na(df[[y]]), ]; N<-nrow(df); sid<-df$sentence_id; S<-n_distinct(sid)
  set.seed(42); fv<-loo::kfold_split_grouped(K=10,x=sid)
  f<-function(rhs) as.formula(paste(y,"~",CTRL,rhs,"+",RE))
  fitkf<-function(nm,form){p<-file.path(CACHE,sprintf("rq3kf_%s_%s.rds",tag,nm)); if(file.exists(p))return(readRDS(p))
    m<-brm(form,data=df,prior=pri,control=list(adapt_delta=0.95,max_treedepth=12),chains=4,iter=2000,warmup=1000,silent=2,refresh=0)
    kf<-kfold(m,folds=fv,chains=4,iter=2000,warmup=1000,silent=2,refresh=0);saveRDS(kf,p);kf}
  pb<-fitkf("base",f(""))$pointwise[,"elpd_kfold"]
  res<-data.frame(outcome=tag,predictor=c("c_nmt","c_mono"),elpd_diff=NA,se_cluster=NA,p=NA)
  for(k in 1:2){v<-c("c_nmt","c_mono")[k]; pt<-fitkf(v,f(paste("+",v)))$pointwise[,"elpd_kfold"]
    di<-pt-pb; ds<-tapply(di,sid,sum); res[k,3:5]<-c(sum(di),sd(ds)*sqrt(S),sfp(ds))}
  res
}
r_gd  <- run_outcome(df_gd,  "log_gd",  "GD")
r_rrt <- run_outcome(df_rrt, "log_rrt", "RRT")
res <- rbind(r_gd, r_rrt)
saveRDS(res, file.path(OUT,"rq3_kfold_elpd.rds"))
cat("=== RQ3 brms kfold elpd (GD, RRT) ===\n"); print(format(res, digits=3))
cat("DONE\n")
