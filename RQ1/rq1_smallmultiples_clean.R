# Clean small-multiples marginal-effect plot (main figure): predicted fixation
# duration vs each continuous predictor, controls at mean, random effects
# marginalised; 95% credible band + rug for coverage, no partial-residual cloud.
OUT <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
suppressMessages({library(brms); library(dplyr); library(ggplot2)})
m <- readRDS(file.path(DATA_DIR,"brm_cache","rq1_coef_maximal.rds"))
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
BLUE<-"#0072B2"
ambref <- factor(levels(df$ambiguity)[1],levels=levels(df$ambiguity))
vars <- c(c_nmt="c[nmt]~(NMT~surprisal)", c_wlen="c[wlen]~(word~length)",
          c_wpos="c[wpos]~(position)",   c_freq="c[freq]~(frequency)")
curve <- bind_rows(lapply(names(vars), function(v){
  g <- data.frame(c_nmt=0,c_wlen=0,c_wpos=0,c_freq=0,ambiguity=ambref)
  xs <- seq(quantile(df[[v]],.005),quantile(df[[v]],.995),length.out=120)
  g <- g[rep(1,length(xs)),]; g[[v]] <- xs
  lp <- posterior_linpred(m,newdata=g,re_formula=NA)
  data.frame(panel=factor(vars[v],levels=vars), x=xs,
             fit=exp(colMeans(lp)), lo=exp(apply(lp,2,quantile,.025)), hi=exp(apply(lp,2,quantile,.975)))
}))
rug <- bind_rows(lapply(names(vars), function(v)
  data.frame(panel=factor(vars[v],levels=vars), x=df[[v]])))
p <- ggplot(curve, aes(x,fit)) +
  geom_ribbon(aes(ymin=lo,ymax=hi), alpha=0.22, fill=BLUE) +
  geom_line(colour=BLUE, linewidth=0.9) +
  geom_rug(data=rug, aes(x=x), inherit.aes=FALSE, alpha=0.012, length=unit(0.03,"npc")) +
  facet_wrap(~panel, scales="free_x", labeller=label_parsed) +
  labs(x="z-scored predictor", y="Predicted fixation duration (ms)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold"),
        axis.line=element_line(colour="black",linewidth=0.3))
ggsave(file.path(OUT,"rq1_effect_smallmultiples.pdf"), p, width=8.0, height=5.6, device="pdf")
cat("saved rq1_effect_smallmultiples.pdf\n")
