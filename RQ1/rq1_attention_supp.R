# Item 1: do attention features add predictive power OVER a baseline that already
# contains c_nmt? Baseline = controls + c_nmt; targets = baseline + each attention.
# brms 10-fold sentence-grouped elpd, sentence-clustered SE, sign-flip perm. Translate stage.
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"; OUT <- "/Users/sebastianx/Dissertation RQ1"
suppressMessages({library(brms); library(dplyr)}); options(mc.cores=4)
fix  <- read.csv(file.path(DATA_DIR,"fixation_durations_word.csv"),stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),stringsAsFactors=FALSE)
attn <- read.csv(file.path(DATA_DIR,"attention_features_6_norm.csv"),stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="") %>%
  select(word=Word,log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))
sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1,.groups="drop")
pred <- nmt %>% left_join(sl,by="sentence_id") %>%
  mutate(word_length=nchar(word),word_position=word_index/(sent_len-1),word_lower=tolower(trimws(word))) %>%
  left_join(freq,by=c("word_lower"="word")) %>%
  left_join(attn %>% select(sentence_id,word_index,H_e=attn_entropy,f_e=attn_context,
                            f_eos=attn_eos,f_recv=attn_recv,f_cross=attn_cross),by=c("sentence_id","word_index")) %>%
  select(sentence_id,word_index,word_length,word_position,log10_freq,nmt_surprisal=surprisal_soft,H_e,f_e,f_eos,f_recv,f_cross)
df <- fix %>% filter(stage=="translate") %>% left_join(pred,by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms),ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal),!is.na(log10_freq),!is.na(H_e),!is.na(f_e),!is.na(f_eos),!is.na(f_recv),!is.na(f_cross)) %>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df<-df%>%mutate(c_nmt=z(nmt_surprisal),c_wlen=z(word_length),c_wpos=z(word_position),c_freq=z(log10_freq),
                c_He=z(H_e),c_fe=z(f_e),c_feos=z(f_eos),c_frecv=z(f_recv),c_fcross=z(f_cross))
N<-nrow(df); sid<-df$sentence_id; S<-n_distinct(sid)
set.seed(42); fold_vec<-loo::kfold_split_grouped(K=10,x=df$sentence_id)
pri<-c(prior(normal(0,1),class=b),prior(normal(6,1),class=Intercept),prior(exponential(1),class=sd),prior(exponential(1),class=sigma))
CACHE<-file.path(DATA_DIR,"brm_cache"); BASE<-"c_wlen+c_wpos+c_freq+ambiguity+c_nmt"; RE<-"(1|participant)+(1|sentence_id)"
f<-function(rhs)as.formula(paste("log_tfd ~",BASE,rhs,"+",RE))
fitkf<-function(nm,form){p<-file.path(CACHE,sprintf("rq1supp_%s.rds",nm)); if(file.exists(p))return(readRDS(p))
  m<-brm(form,data=df,prior=pri,control=list(adapt_delta=0.95,max_treedepth=12),chains=4,iter=2000,warmup=1000,silent=2,refresh=0)
  kf<-kfold(m,folds=fold_vec,chains=4,iter=2000,warmup=1000,silent=2,refresh=0); saveRDS(kf,p); kf}
pb<-fitkf("base_cnmt",f(""))$pointwise[,"elpd_kfold"]
sfp<-function(ds,n=10000){o<-sum(ds);set.seed(42);mean(replicate(n,sum(ds*sample(c(-1,1),length(ds),replace=TRUE)))>=o)}
vars<-c("c_He","c_fe","c_feos","c_frecv","c_fcross"); labs<-c("H_e","f_e","f_eos","f_recv","f_cross")
res<-data.frame(feature=labs,elpd_diff=NA,se_cluster=NA,p=NA,per_word=NA)
for(k in seq_along(vars)){pt<-fitkf(vars[k],f(paste("+",vars[k])))$pointwise[,"elpd_kfold"]
  di<-pt-pb; ds<-tapply(di,sid,sum); res[k,-1]<-c(sum(di),sd(ds)*sqrt(S),sfp(ds),sum(di)/N)}
res$p_holm<-p.adjust(res$p,method="holm")
saveRDS(res,file.path(OUT,"rq1_attention_supp.rds"))
cat("=== Attention supplementary power over controls+c_nmt (translate) ===\n"); print(format(res,digits=3))
cat("DONE\n")
