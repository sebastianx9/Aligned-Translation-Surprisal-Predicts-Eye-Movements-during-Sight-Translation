# Quick exploratory check: do the five attention features predict GD / RRT?
# lmer, 10-fold sentence-grouped held-out Delta-llh over a controls baseline
# (per-model sigma; predict with participant RE, marginalise sentence RE;
# sentence-level sign-flip permutation; Holm across the 5 features per outcome).
# Normalised attention file only (per methods decision).
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(dplyr); library(lme4)})
em   <- read.csv(file.path(DATA_DIR,"eye_measures_word.csv"),         stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),   stringsAsFactors=FALSE)
attn <- read.csv(file.path(DATA_DIR,"attention_features_6_norm.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="") %>%
  select(word=Word,log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))
sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1,.groups="drop")
pred <- nmt %>% left_join(sl,by="sentence_id") %>%
  mutate(word_length=nchar(word),word_position=word_index/(sent_len-1),word_lower=tolower(trimws(word))) %>%
  left_join(freq,by=c("word_lower"="word")) %>%
  left_join(attn %>% select(sentence_id,word_index,H_e=attn_entropy,f_e=attn_context,
                            f_eos=attn_eos,f_recv=attn_recv,f_cross=attn_cross),by=c("sentence_id","word_index")) %>%
  select(sentence_id,word_index,word_length,word_position,log10_freq,H_e,f_e,f_eos,f_recv,f_cross)
base <- em %>% filter(stage=="translate") %>% left_join(pred,by=c("sentence_id","word_index")) %>%
  mutate(ambiguity=factor(ambiguity)) %>%
  filter(!is.na(log10_freq),!is.na(H_e),!is.na(f_e),!is.na(f_eos),!is.na(f_recv),!is.na(f_cross)) %>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
base <- base %>% mutate(c_wlen=z(word_length),c_wpos=z(word_position),c_freq=z(log10_freq),
                        c_He=z(H_e),c_fe=z(f_e),c_feos=z(f_eos),c_frecv=z(f_recv),c_fcross=z(f_cross))
feats <- c(H_e="c_He",f_e="c_fe",f_eos="c_feos",f_recv="c_frecv",f_cross="c_fcross")
CTRL <- "c_wlen+c_wpos+c_freq+ambiguity"; RE <- "(1|participant)+(1|sentence_id)"
ctrl <- lmerControl(optimizer="bobyqa")

heldout_llh <- function(df, rhs){
  set.seed(42); fold <- loo::kfold_split_grouped(K=10, x=df$sentence_id)
  ll <- numeric(nrow(df))
  for(k in 1:10){ tr<-df[fold!=k,]; te<-df[fold==k,]
    m<-lmer(as.formula(paste("y ~",CTRL,rhs,"+",RE)),data=tr,REML=FALSE,control=ctrl)
    mu<-predict(m,newdata=te,allow.new.levels=TRUE,re.form=~(1|participant))
    ll[fold==k]<-dnorm(te$y,mu,sigma(m),log=TRUE) }
  ll }
signflip <- function(ds,n=5000){o<-sum(ds);set.seed(42);mean(replicate(n,sum(ds*sample(c(-1,1),length(ds),TRUE)))>=o)}

run <- function(df,label){
  df<-df[!is.na(df$y),]; sid<-df$sentence_id; S<-length(unique(sid))
  b<-heldout_llh(df,"")
  res<-data.frame(outcome=label,feature=names(feats),dllh=NA,se=NA,p=NA)
  for(i in seq_along(feats)){ t<-heldout_llh(df,paste("+",feats[i]))
    di<-t-b; ds<-tapply(di,sid,sum); res[i,3:5]<-c(mean(di),sd(ds)*sqrt(S)/nrow(df),signflip(ds)) }
  res$p_holm<-p.adjust(res$p,method="holm"); res }
gd  <- run(base %>% filter(gd_ms>0)  %>% mutate(y=log(gd_ms)),  "GD")
rrt <- run(base %>% filter(reread_occurrence==1, rrt_ms>0) %>% mutate(y=log(rrt_ms)), "RRT")
cat("\n=== Attention features -> GD / RRT (lmer 10-fold held-out Delta-llh) ===\n")
print(format(rbind(gd,rrt),digits=3)); cat("DONE\n")
