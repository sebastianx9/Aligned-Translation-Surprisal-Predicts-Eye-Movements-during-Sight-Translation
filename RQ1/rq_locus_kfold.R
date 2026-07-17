# ── Locus family: predictive (kfold elpd) tests of neighbour-surprisal effects ──
# Option B: separate hypothesis family ("where does the mono signal fall"),
# distinct from RQ1's 7-candidate family. All on interior words (both
# neighbours present), with baselines refit on the same subsets.
#
#   A) reading TFD  : base | +s_prev | +s_curr | +s_next   (locus of mono, oral reading)
#   B) GD  (translate): base | +s_next                      (preview in first pass)
#   C) RRT (translate): base | +s_next                      (successor in re-reading)
#   D) translate TFD : base+window(3 mono) | +c_nmt         (c_nmt beyond the EVS window)
#
# Same conventions as rq1/rq3 kfold scripts: sentence-grouped 10-fold (seed 42),
# identical RE (1|participant)+(1|sentence_id), sentence-clustered SE,
# 10,000-permutation sign-flip test.
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT      <- "/Users/sebastianx/Dissertation RQ1"
suppressMessages({library(brms); library(dplyr)}); options(mc.cores=4)

fix  <- read.csv(file.path(DATA_DIR,"fixation_durations_word.csv"),    stringsAsFactors=FALSE)
em   <- read.csv(file.path(DATA_DIR,"eye_measures_word.csv"),          stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR,"nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR,"monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR,"subtlex_us.csv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="") %>%
  select(word=Word,log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))

m <- mono %>% select(sentence_id, word_index, mono_curr=surprisal_sum) %>%
  group_by(sentence_id) %>% arrange(word_index,.by_group=TRUE) %>%
  mutate(mono_prev=lag(mono_curr), mono_next=lead(mono_curr)) %>% ungroup()
sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1,.groups="drop")
pred <- nmt %>% left_join(sl,by="sentence_id") %>%
  mutate(word_length=nchar(word), word_position=word_index/(sent_len-1),
         word_lower=tolower(trimws(word))) %>%
  left_join(freq,by=c("word_lower"="word")) %>%
  left_join(m,by=c("sentence_id","word_index")) %>%
  select(sentence_id,word_index,word_length,word_position,log10_freq,
         nmt_surprisal=surprisal_soft,mono_curr,mono_prev,mono_next)

prep <- function(df, ycol) {
  df <- df %>% left_join(pred,by=c("sentence_id","word_index")) %>%
    mutate(ambiguity=factor(ambiguity)) %>%
    filter(!is.na(nmt_surprisal),!is.na(mono_curr),!is.na(mono_prev),
           !is.na(mono_next),!is.na(log10_freq)) %>%
    anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
  z <- function(x)(x-mean(x))/sd(x)
  df %>% mutate(y=log(.data[[ycol]]),
                c_nmt=z(nmt_surprisal), c_mono=z(mono_curr),
                c_prev=z(mono_prev), c_next=z(mono_next),
                c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))
}

pri <- c(prior(normal(0,1),class=b), prior(normal(6,1),class=Intercept),
         prior(exponential(1),class=sd), prior(exponential(1),class=sigma))
CACHE <- file.path(DATA_DIR,"brm_cache")
CTRL <- "c_wlen + c_wpos + c_freq + ambiguity"
RE   <- "(1|participant) + (1|sentence_id)"

sfp <- function(ds, n=10000){o<-sum(ds); set.seed(42)
  mean(replicate(n, sum(ds*sample(c(-1,1),length(ds),replace=TRUE)))>=o)}

run_family <- function(df, tag, rhs_list) {
  N<-nrow(df); sid<-df$sentence_id; S<-n_distinct(sid)
  set.seed(42); fv<-loo::kfold_split_grouped(K=10, x=sid)
  cat(sprintf("\n══════ %s: N=%d, S=%d ══════\n", tag, N, S))
  fitkf <- function(nm, rhs){
    p <- file.path(CACHE, sprintf("locus_%s_%s.rds", tag, nm))
    if (file.exists(p)) { cat(sprintf("[%s/%s] cached\n",tag,nm)); return(readRDS(p)) }
    cat(sprintf("[%s/%s] fitting + 10-fold ...\n", tag, nm)); t0<-proc.time()
    mm <- brm(as.formula(paste("y ~",CTRL,rhs,"+",RE)), data=df, prior=pri,
              control=list(adapt_delta=0.95,max_treedepth=12),
              chains=4, iter=2000, warmup=1000, silent=2, refresh=0)
    kf <- kfold(mm, folds=fv, chains=4, iter=2000, warmup=1000, silent=2, refresh=0)
    saveRDS(kf,p); cat(sprintf("[%s/%s] done %.0f min\n",tag,nm,(proc.time()-t0)["elapsed"]/60)); kf
  }
  pb <- fitkf("base", rhs_list$base)$pointwise[,"elpd_kfold"]
  res <- list()
  for (nm in setdiff(names(rhs_list),"base")) {
    pt <- fitkf(nm, rhs_list[[nm]])$pointwise[,"elpd_kfold"]
    di <- pt-pb; ds <- tapply(di,sid,sum)
    res[[nm]] <- c(elpd=sum(di), se=sd(ds)*sqrt(S), p=sfp(ds))
    cat(sprintf("  %s vs base: elpd=%+7.3f  se=%6.3f  p=%.4f%s\n",
        nm, sum(di), sd(ds)*sqrt(S), res[[nm]]["p"],
        ifelse(res[[nm]]["p"]<.05,"  *","")))
  }
  res
}

d_read <- prep(fix %>% filter(stage=="read", total_fixation_duration_ms>0), "total_fixation_duration_ms")
d_trans<- prep(fix %>% filter(stage=="translate", total_fixation_duration_ms>0), "total_fixation_duration_ms")
d_gd   <- prep(em %>% filter(stage=="translate", !is.na(gd_ms), gd_ms>0), "gd_ms")
d_rrt  <- prep(em %>% filter(stage=="translate", regress_in==1, rrt_ms>0), "rrt_ms")

results <- list(
  reading = run_family(d_read, "readTFD",
    list(base="", s_prev="+ c_prev", s_curr="+ c_mono", s_next="+ c_next")),
  gd      = run_family(d_gd, "GD",   list(base="", s_next="+ c_next")),
  rrt     = run_family(d_rrt, "RRT", list(base="", s_next="+ c_next")),
  cnmt_beyond_window = run_family(d_trans, "transWIN",
    list(base="+ c_prev + c_mono + c_next", c_nmt="+ c_prev + c_mono + c_next + c_nmt"))
)
saveRDS(results, file.path(OUT,"rq_locus_kfold.rds"))
cat("\nSaved rq_locus_kfold.rds\nALL LOCUS FAMILIES DONE\n")
