# Stage decomposition of BOTH interaction gains (condition:c_nmt, condition:c_mono).
# condition:c_nmt gain = M3 - M2   (M2, M3 already cached)
# condition:c_mono gain = M3 - M3_nomonoint, where M3_nomonoint drops cond:c_mono
#   (M3 minus the c_mono stage interaction) -- fit here.
# No stratified FITTING: all models are fit jointly on the pooled data; we only
# partition each model pair's already-computed pointwise held-out density by stage.
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"; OUT <- "/Users/sebastianx/Dissertation RQ2"
suppressMessages({library(brms); library(dplyr)}); options(mc.cores=4)

kf <- readRDS(file.path(OUT,"rq2_kfold_elpd.rds")); ptw <- kf$pointwise; N <- kf$N; sid_saved <- kf$sentence_id

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
  mutate(log_tfd=log(total_fixation_duration_ms),ambiguity=factor(ambiguity),cond=as.integer(stage=="translate")) %>%
  filter(!is.na(nmt_surprisal),!is.na(mono_surprisal),!is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal),c_mono=z(mono_surprisal),c_wlen=z(word_length),c_wpos=z(word_position),c_freq=z(log10_freq))
stopifnot(nrow(df)==N, all(df$sentence_id==sid_saved))
sid <- df$sentence_id; cond <- df$cond

set.seed(42); fold_vec <- loo::kfold_split_grouped(K=10, x=df$sentence_id)   # EXACT reproduction of rq2_kfold_elpd.R folds
priors <- c(prior(normal(0,1),class=b),prior(normal(6,1),class=Intercept),prior(exponential(1),class=sd),prior(exponential(1),class=sigma))
CACHE <- file.path(DATA_DIR,"brm_cache")
form_nomonoint <- log_tfd ~ c_wlen+c_wpos+c_freq+ambiguity + cond + c_mono + c_nmt + cond:c_nmt + (1|participant)+(1|sentence_id)
p_no <- file.path(CACHE,"kfold_M3_nomonoint.rds")
if (file.exists(p_no)) { cat("[M3_nomonoint] cached\n"); kf_no <- readRDS(p_no) } else {
  cat("[M3_nomonoint] fitting + 10-fold ...\n"); t0<-proc.time()
  m <- brm(form_nomonoint,data=df,prior=priors,control=list(adapt_delta=0.95,max_treedepth=12),chains=4,iter=2000,warmup=1000,silent=2,refresh=0)
  kf_no <- kfold(m,folds=fold_vec,chains=4,iter=2000,warmup=1000,silent=2,refresh=0)
  saveRDS(kf_no,p_no); cat(sprintf("[M3_nomonoint] %.0f min\n",(proc.time()-t0)["elapsed"]/60)) }
ptw_no <- kf_no$pointwise[,"elpd_kfold"]

signflip <- function(d_s,n=10000){o<-sum(d_s);set.seed(42);mean(replicate(n,sum(d_s*sample(c(-1,1),length(d_s),replace=TRUE)))>=o)}
decomp <- function(d_i,arm){
  do.call(rbind,lapply(list(c("pooled",NA),c("translation",1),c("reading",0)),function(s){
    mask <- if(is.na(s[2])) rep(TRUE,N) else cond==as.integer(s[2])
    di<-d_i[mask];ds<-tapply(di,sid[mask],sum);S<-length(ds)
    data.frame(interaction=arm,subset=s[1],n=sum(mask),elpd_diff=sum(di),se=sd(ds)*sqrt(S),z=sum(di)/(sd(ds)*sqrt(S)),p=signflip(ds))}))}
res <- rbind(decomp(ptw$M3-ptw$M2,"condition:c_nmt"), decomp(ptw$M3-ptw_no,"condition:c_mono"))
res$sig <- ifelse(res$p<.05,"*","")
cat("\n=== Interaction gains decomposed by stage ===\n"); print(format(res,digits=3,nsmall=3))
saveRDS(res, file.path(OUT,"rq2_interaction_decomp.rds")); cat("\nSaved rq2_interaction_decomp.rds\nDONE\n")
