# Item 7: symmetric converse. N1=controls+cond+c_nmt+cond:c_nmt (cached);
# N2=N1+c_mono; N3=N2+cond:c_mono=M3 (cached). Gives c_mono main (N2|N1) and
# c_mono interaction (M3|N2), mirroring c_nmt's M2|M1 and M3|M2.
DATA_DIR<-"/Users/sebastianx/Dissertation_Data"; OUT<-"/Users/sebastianx/Dissertation RQ2"
suppressMessages({library(brms);library(dplyr)}); options(mc.cores=4)
fix<-read.csv(file.path(DATA_DIR,"fixation_durations_word.csv"),stringsAsFactors=FALSE)
nmt<-read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),stringsAsFactors=FALSE)
mono<-read.csv(file.path(DATA_DIR,"monolingual_surprisal_word.csv"),stringsAsFactors=FALSE)
freq<-read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="")%>%
  select(word=Word,log10_freq=Lg10WF)%>%mutate(word=tolower(trimws(word)))
sl<-nmt%>%group_by(sentence_id)%>%summarise(sent_len=max(word_index)+1,.groups="drop")
pred<-nmt%>%left_join(sl,by="sentence_id")%>%
  mutate(word_length=nchar(word),word_position=word_index/(sent_len-1),word_lower=tolower(trimws(word)))%>%
  left_join(freq,by=c("word_lower"="word"))%>%
  left_join(mono%>%select(sentence_id,word_index,mono_surprisal=surprisal_sum),by=c("sentence_id","word_index"))%>%
  select(sentence_id,word_index,word_length,word_position,nmt_surprisal=surprisal_soft,mono_surprisal,log10_freq)
df<-fix%>%filter(stage%in%c("translate","read"))%>%left_join(pred,by=c("sentence_id","word_index"))%>%
  mutate(log_tfd=log(total_fixation_duration_ms),ambiguity=factor(ambiguity),cond=as.integer(stage=="translate"))%>%
  filter(!is.na(nmt_surprisal),!is.na(mono_surprisal),!is.na(log10_freq))%>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df<-df%>%mutate(c_nmt=z(nmt_surprisal),c_mono=z(mono_surprisal),c_wlen=z(word_length),c_wpos=z(word_position),c_freq=z(log10_freq))
N<-nrow(df);sid<-df$sentence_id;S<-n_distinct(sid)
set.seed(42);fold_vec<-loo::kfold_split_grouped(K=10,x=df$sentence_id)
pri<-c(prior(normal(0,1),class=b),prior(normal(6,1),class=Intercept),prior(exponential(1),class=sd),prior(exponential(1),class=sigma))
CACHE<-file.path(DATA_DIR,"brm_cache")
RE<-"(1|participant)+(1|sentence_id)"; CTRL<-"c_wlen+c_wpos+c_freq+ambiguity"
fN2<-as.formula(paste("log_tfd ~",CTRL,"+ cond + c_nmt + cond:c_nmt + c_mono +",RE))
p<-file.path(CACHE,"kfold_N2.rds")
if(file.exists(p)){kfN2<-readRDS(p)}else{
  m<-brm(fN2,data=df,prior=pri,control=list(adapt_delta=0.95,max_treedepth=12),chains=4,iter=2000,warmup=1000,silent=2,refresh=0)
  kfN2<-kfold(m,folds=fold_vec,chains=4,iter=2000,warmup=1000,silent=2,refresh=0);saveRDS(kfN2,p)}
ptN2<-kfN2$pointwise[,"elpd_kfold"]
ptN1<-readRDS(file.path(CACHE,"kfold_N1.rds"))$pointwise[,"elpd_kfold"]
ptM3<-readRDS(file.path(CACHE,"kfold_M3.rds"))$pointwise[,"elpd_kfold"]
sfp<-function(ds,n=10000){o<-sum(ds);set.seed(42);mean(replicate(n,sum(ds*sample(c(-1,1),length(ds),replace=TRUE)))>=o)}
rep<-function(nm,di){ds<-tapply(di,sid,sum);cat(sprintf("%-14s elpd_diff=%+.2f se=%.2f p=%.3f per-word=%+.5f\n",nm,sum(di),sd(ds)*sqrt(S),sfp(ds),sum(di)/N))}
cat("=== c_mono mirror ladder ===\n")
rep("N2|N1 c_mono main", ptN2-ptN1)
rep("M3|N2 c_mono int",  ptM3-ptN2)
saveRDS(list(N2N1=ptN2-ptN1, M3N2=ptM3-ptN2, sid=sid, N=N, S=S), file.path(OUT,"rq2_mirror_ladder.rds"))
cat("DONE\n")
