# Two variants of the RQ1 partial-residual display, from the maximal-RE model:
#  (1) c_nmt effect plot with a loess overlay (empirical shape vs linear fit)
#  (2) small multiples: partial-residual effect plots for all continuous predictors
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

fe <- fixef(m); b0 <- fe["Intercept","Estimate"]
full <- fitted(m)[,"Estimate"]; resid <- df$log_tfd - full
BLUE <- "#0072B2"

## (1) c_nmt with loess overlay ------------------------------------------------
df$presid <- b0 + fe["c_nmt","Estimate"]*df$c_nmt + resid
grid <- data.frame(c_nmt=seq(min(df$c_nmt),max(df$c_nmt),length.out=120),
                   c_wlen=0,c_wpos=0,c_freq=0,ambiguity=factor(levels(df$ambiguity)[1],levels=levels(df$ambiguity)))
lp <- posterior_linpred(m,newdata=grid,re_formula=NA)
grid$fit<-colMeans(lp); grid$lo<-apply(lp,2,quantile,.025); grid$hi<-apply(lp,2,quantile,.975)
p1 <- ggplot() +
  geom_point(data=df, aes(c_nmt,presid), alpha=0.05, size=0.5, colour=BLUE) +
  geom_ribbon(data=grid, aes(c_nmt,ymin=lo,ymax=hi), alpha=0.22, fill=BLUE) +
  geom_line(data=grid, aes(c_nmt,fit), colour=BLUE, linewidth=0.9) +
  geom_smooth(data=df, aes(c_nmt,presid), method="loess", se=FALSE,
              colour="grey20", linetype="dashed", linewidth=0.8) +
  labs(x=expression(c[nmt]*"  (z-scored aligned NMT surprisal)"),
       y="Partial residual  (log fixation duration, ms)") +
  theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(), axis.line=element_line(colour="black",linewidth=0.3))
ggsave(file.path(OUT,"rq1_cnmt_partial_smooth.pdf"), p1, width=6.5, height=4.4, device="pdf")

## (2) small multiples for all continuous predictors ---------------------------
vars <- c(c_nmt="c[nmt]~(NMT~surprisal)", c_wlen="c[wlen]~(word~length)",
          c_wpos="c[wpos]~(position)", c_freq="c[freq]~(frequency)")
bk <- sapply(names(vars), function(v) fe[v,"Estimate"])
long <- bind_rows(lapply(names(vars), function(v)
  data.frame(panel=factor(vars[v],levels=vars), x=df[[v]], presid=b0 + bk[v]*df[[v]] + resid)))
p2 <- ggplot(long, aes(x, presid)) +
  geom_point(alpha=0.04, size=0.4, colour=BLUE) +
  geom_smooth(method="lm", se=TRUE, colour=BLUE, fill=BLUE, linewidth=0.9, alpha=0.2) +
  geom_smooth(method="loess", se=FALSE, colour="grey20", linetype="dashed", linewidth=0.7) +
  facet_wrap(~panel, scales="free_x", labeller=label_parsed) +
  labs(x="z-scored predictor", y="Partial residual  (log fixation duration, ms)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold"),
        axis.line=element_line(colour="black",linewidth=0.3))
ggsave(file.path(OUT,"rq1_partial_smallmultiples.pdf"), p2, width=8.0, height=6.0, device="pdf")
cat("saved rq1_cnmt_partial_smooth.pdf and rq1_partial_smallmultiples.pdf\n")
cat(sprintf("slopes: c_nmt=%.3f c_wlen=%.3f c_wpos=%.3f c_freq=%.3f\n", bk["c_nmt"],bk["c_wlen"],bk["c_wpos"],bk["c_freq"]))
