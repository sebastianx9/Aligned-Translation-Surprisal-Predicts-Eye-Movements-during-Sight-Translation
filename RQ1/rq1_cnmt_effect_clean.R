# Clean marginal-effect plot for c_nmt (ggpredict style): predicted fixation
# duration vs c_nmt with 95% credible band, controls at their mean and random
# effects marginalised. No partial-residual cloud; a rug shows data coverage.
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
grid <- data.frame(c_nmt=seq(quantile(df$c_nmt,.005),quantile(df$c_nmt,.995),length.out=140),
                   c_wlen=0,c_wpos=0,c_freq=0,ambiguity=factor(levels(df$ambiguity)[1],levels=levels(df$ambiguity)))
lp <- posterior_linpred(m,newdata=grid,re_formula=NA)
grid$fit<-exp(colMeans(lp)); grid$lo<-exp(apply(lp,2,quantile,.025)); grid$hi<-exp(apply(lp,2,quantile,.975))
p <- ggplot(grid, aes(c_nmt,fit)) +
  geom_ribbon(aes(ymin=lo,ymax=hi), alpha=0.22, fill=BLUE) +
  geom_line(colour=BLUE, linewidth=1.0) +
  geom_rug(data=df, aes(x=c_nmt), inherit.aes=FALSE, alpha=0.015, length=unit(0.03,"npc")) +
  labs(x=expression(c[nmt]*"  (z-scored aligned NMT surprisal)"),
       y="Predicted fixation duration (ms)") +
  theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(), axis.line=element_line(colour="black",linewidth=0.3))
ggsave(file.path(OUT,"rq1_cnmt_effect_clean.pdf"), p, width=6.5, height=4.2, device="pdf")
cat(sprintf("saved rq1_cnmt_effect_clean.pdf | ms range %.0f--%.0f\n", min(grid$fit), max(grid$fit)))
