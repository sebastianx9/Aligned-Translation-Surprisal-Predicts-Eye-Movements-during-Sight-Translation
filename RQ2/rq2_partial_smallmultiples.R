# RQ2 data-side figure: partial residuals of log(TFD) against ALL continuous
# predictors (both surprisals + the three controls), one panel each, from the
# MAXIMAL joint model. Points + lines are split by stage; solid = linear fit,
# dashed = loess. The interaction shows directly: the two surprisals' stage
# lines DIVERGE (stage-dependent slope) while the controls' COINCIDE (no
# stage interaction), and c_nmt's slope is seen against the control slopes.
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(brms); library(dplyr); library(ggplot2)})
m <- readRDS(file.path(DATA_DIR,"brm_cache","rq2_joint_maximal.rds"))

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

fe <- fixef(m); b0<-fe["Intercept","Estimate"]
full <- fitted(m)[,"Estimate"]; resid <- df$log_tfd - full
Stage <- factor(ifelse(df$condition==1,"Translation","Reading"), levels=c("Reading","Translation"))

# per-predictor stage-specific slope (surprisals interact with stage; controls do not)
slopes <- list(
  "c[nmt]"  = fe["c_nmt","Estimate"]  + df$condition*fe["condition:c_nmt","Estimate"],
  "c[mono]" = fe["c_mono","Estimate"] + df$condition*fe["condition:c_mono","Estimate"],
  "c[wlen]" = rep(fe["c_wlen","Estimate"], nrow(df)),
  "c[wpos]" = rep(fe["c_wpos","Estimate"], nrow(df)),
  "c[freq]" = rep(fe["c_freq","Estimate"], nrow(df)))
xvars <- list("c[nmt]"=df$c_nmt,"c[mono]"=df$c_mono,"c[wlen]"=df$c_wlen,"c[wpos]"=df$c_wpos,"c[freq]"=df$c_freq)
lv <- c("c[nmt]","c[mono]","c[wlen]","c[wpos]","c[freq]")
long <- bind_rows(lapply(lv, function(v)
  data.frame(predictor=v, stage=Stage, x=xvars[[v]], y=b0 + slopes[[v]]*xvars[[v]] + resid))) %>%
  mutate(predictor=factor(predictor, levels=lv))
scol <- c("Reading"="#D55E00","Translation"="#0072B2")

# solid line = MODEL-implied fit (stage-invariant slope for controls -> the two
# stage lines coincide; stage-specific for the surprisals -> they diverge)
mslope <- function(v, cond) switch(v,
  "c[nmt]" =fe["c_nmt","Estimate"] +cond*fe["condition:c_nmt","Estimate"],
  "c[mono]"=fe["c_mono","Estimate"]+cond*fe["condition:c_mono","Estimate"],
  "c[wlen]"=fe["c_wlen","Estimate"], "c[wpos]"=fe["c_wpos","Estimate"], "c[freq]"=fe["c_freq","Estimate"])
xr <- long %>% group_by(predictor) %>% summarise(xmin=min(x),xmax=max(x),.groups="drop")
mline <- bind_rows(lapply(lv, function(v) bind_rows(lapply(c("Reading","Translation"), function(st){
  r<-xr[xr$predictor==v,]; xs<-seq(r$xmin,r$xmax,length.out=60)
  data.frame(predictor=factor(v,levels=lv), stage=factor(st,levels=c("Reading","Translation")),
             x=xs, y=b0 + mslope(v, as.integer(st=="Translation"))*xs) }))))

p <- ggplot(long, aes(x,y,colour=stage)) +
  geom_point(alpha=0.03, size=0.35) +
  geom_smooth(method="loess", se=FALSE, linewidth=0.7, linetype="dashed") +
  geom_line(data=mline, linewidth=0.9) +
  facet_wrap(~predictor, scales="free_x", labeller=label_parsed, ncol=3) +
  scale_colour_manual(values=scol, name=NULL) +
  labs(x="z-scored predictor (SD units)", y="Partial residual of log(TFD)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold"),
        legend.position="top", axis.line=element_line(colour="black",linewidth=0.3))
ggsave(file.path(OUT,"rq2_interaction_smallmultiples.pdf"), p, width=8.4, height=5.6, device="pdf")
cat("saved rq2_interaction_smallmultiples.pdf (5 predictors, stage-split, lm+loess)\n")
