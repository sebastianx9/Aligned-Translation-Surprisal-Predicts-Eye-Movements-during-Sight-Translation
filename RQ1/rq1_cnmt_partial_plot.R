# RQ1 partial-residual (effect) plot for c_nmt, from the MAXIMAL-RE coefficient model.
# Partial residual = Intercept + b_cnmt*c_nmt + (log_tfd - full_fitted); the line is
# the fixed-effect c_nmt slope with controls at their mean and random effects at 0.
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

fe <- fixef(m); b0 <- fe["Intercept","Estimate"]; b1 <- fe["c_nmt","Estimate"]
full <- fitted(m)[,"Estimate"]                      # posterior-mean expectation incl. random effects
df$presid <- b0 + b1*df$c_nmt + (df$log_tfd - full) # component + residual

grid <- data.frame(c_nmt=seq(min(df$c_nmt),max(df$c_nmt),length.out=120),
                   c_wlen=0,c_wpos=0,c_freq=0,
                   ambiguity=factor(levels(df$ambiguity)[1],levels=levels(df$ambiguity)))
lp <- posterior_linpred(m,newdata=grid,re_formula=NA)
grid$fit<-colMeans(lp); grid$lo<-apply(lp,2,quantile,.025); grid$hi<-apply(lp,2,quantile,.975)

qs <- quantile(df$c_nmt, seq(0,1,.1)); df$bin <- cut(df$c_nmt,qs,include.lowest=TRUE)
binned <- df %>% group_by(bin) %>% summarise(x=mean(c_nmt),y=mean(presid),se=sd(presid)/sqrt(n()),.groups="drop")

p <- ggplot() +
  geom_point(data=df, aes(c_nmt,presid), alpha=0.05, size=0.5, colour="#0072B2") +
  geom_ribbon(data=grid, aes(c_nmt,ymin=lo,ymax=hi), alpha=0.22, fill="#0072B2") +
  geom_line(data=grid, aes(c_nmt,fit), colour="#0072B2", linewidth=0.9) +
  labs(x=expression(c[nmt]*"  (z-scored aligned NMT surprisal)"),
       y="Partial residual  (log fixation duration, ms)") +
  theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(), axis.line=element_line(colour="black",linewidth=0.3))
ggsave(file.path(OUT,"rq1_cnmt_partial.pdf"), p, width=6.5, height=4.4, device="pdf")
cat(sprintf("saved rq1_cnmt_partial.pdf | slope b1=%.4f | n=%d | decile means:\n", b1, nrow(df)))
print(binned %>% mutate(across(where(is.numeric),~round(.,3))))
