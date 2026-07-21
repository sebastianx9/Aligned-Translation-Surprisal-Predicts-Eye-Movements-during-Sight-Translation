#!/usr/bin/env Rscript

# RQ3 brms kfold elpd — c_nmt / c_mono on FFD, GD, go-past, and conditional RRT.
# Bayesian replacement for the lmer two-stage LOO-CV. Same convention as
# rq1_kfold_elpd.R: 10-fold sentence-grouped elpd, sentence-clustered SE,
# sentence-level sign-flip permutation. Fixes the two issues in the old lmer
# script (shared sigma(m_base); it also used the un-normalised attention file,
# irrelevant here since only c_nmt/c_mono are reported).
suppressMessages({library(brms); library(dplyr)}); options(mc.cores = 4)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork=TRUE)
source(file.path(repo_root, "R", "analysis_design.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}
DATA_DIR <- normalizePath(
  get_arg("--data-dir", Sys.getenv("DISSERTATION_DATA_DIR", ".")),
  mustWork = TRUE
)
OUT <- get_arg("--output-dir", Sys.getenv("DISSERTATION_OUTPUT_DIR", DATA_DIR))
exclude_contrastive <- parse_bool(
  get_arg("--exclude-contrastive", "false"), "--exclude-contrastive"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

eye_path <- file.path(DATA_DIR, "eye_measures_word.csv")
nmt_path <- file.path(DATA_DIR, "nmt_surprisal_soft_word.csv")
mono_path <- file.path(DATA_DIR, "monolingual_surprisal_word.csv")
freq_path <- file.path(DATA_DIR, "subtlex_us.csv")
input_hashes <- analysis_input_hashes(c(
  eye_measures=eye_path, nmt_surprisal=nmt_path,
  monolingual_surprisal=mono_path, frequency=freq_path,
  analysis_design=file.path(repo_root, "R", "analysis_design.R")
))

em   <- read.csv(eye_path,  stringsAsFactors=FALSE)
nmt  <- read.csv(nmt_path,  stringsAsFactors=FALSE)
mono <- read.csv(mono_path, stringsAsFactors=FALSE)
freq <- read.table(freq_path, sep="\t", header=TRUE,
                   stringsAsFactors=FALSE, quote="") %>%
  transmute(word_lower=lexical_form(Word),log10_freq=Lg10WF)
required_eye_columns <- c(
  "participant", "sentence_id", "stage", "word_index",
  "ffd_ms", "gd_ms", "go_past_ms", "rrt_ms", "reread_occurrence",
  "first_encounter_status"
)
missing_eye_columns <- setdiff(required_eye_columns, names(em))
if (length(missing_eye_columns)) {
  stop(
    "eye_measures_word.csv must be regenerated with the current extractor; ",
    "missing columns: ", paste(missing_eye_columns, collapse=", ")
  )
}
sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1,.groups="drop")
pred <- nmt %>% left_join(sl,by="sentence_id") %>%
  mutate(word_lower=lexical_form(word),word_length=nchar(word_lower, type="chars"),
         word_position=word_index/(sent_len-1)) %>%
  left_join(freq,by="word_lower") %>%
  left_join(mono %>% select(sentence_id,word_index,mono_surprisal=surprisal_sum),by=c("sentence_id","word_index")) %>%
  select(sentence_id,word_index,word_length,word_position,log10_freq,nmt_surprisal=surprisal_soft,mono_surprisal)
base_df <- em %>% filter(stage=="translate") %>% left_join(pred,by=c("sentence_id","word_index")) %>%
  mutate(ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal),!is.na(mono_surprisal),!is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003",word_index=3L),by=c("sentence_id","word_index"))
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
base_df <- base_df %>% mutate(c_nmt=z(nmt_surprisal),c_mono=z(mono_surprisal),
                              c_wlen=z(word_length),c_wpos=z(word_position),c_freq=z(log10_freq))
if (any(
  !is.na(base_df$go_past_ms) &
    (!is.finite(base_df$go_past_ms) | base_df$go_past_ms <= 0)
)) {
  stop("Observed go_past_ms values must be finite and positive.")
}
if (any(
  !is.na(base_df$go_past_ms) & !is.na(base_df$gd_ms) &
    base_df$go_past_ms + 0.01 < base_df$gd_ms
)) {
  stop("Observed go_past_ms cannot be shorter than gd_ms.")
}
df_ffd <- base_df %>% filter(!is.na(ffd_ms), ffd_ms > 0) %>% mutate(log_ffd = log(ffd_ms))
df_gd  <- base_df %>% filter(!is.na(gd_ms),  gd_ms  > 0) %>% mutate(log_gd  = log(gd_ms))
df_go_past <- base_df %>%
  filter(!is.na(go_past_ms), go_past_ms > 0) %>%
  mutate(log_go_past = log(go_past_ms))
df_rrt <- base_df %>%
  filter(reread_occurrence == 1, rrt_ms > 0) %>%
  mutate(log_rrt = log(rrt_ms))
cat(sprintf("FFD n=%d  GD n=%d  go-past n=%d  conditional RRT n=%d\n",
            nrow(df_ffd), nrow(df_gd), nrow(df_go_past), nrow(df_rrt)))

pri<-c(prior(normal(0,1),class=b),prior(normal(6,1),class=Intercept),prior(exponential(1),class=sd),prior(exponential(1),class=sigma))
CACHE<-file.path(DATA_DIR,"brm_cache"); dir.create(CACHE, showWarnings=FALSE)
CTRL<-"c_wlen+c_wpos+c_freq+ambiguity"; RE<-"(1|participant)+(1|sentence_id)"
sfp<-function(ds,n=10000){o<-sum(ds);set.seed(42)
  permuted<-replicate(n,sum(ds*sample(c(-1,1),length(ds),replace=TRUE)))
  (1+sum(permuted>=o))/(n+1)}

run_outcome <- function(df, y, tag) {
  df <- df[!is.na(df[[y]]), ]
  fv_full <- make_sentence_folds(df$sentence_id, K=10, seed=42)
  keep <- contrastive_keep(df$sentence_id, exclude_contrastive)
  df <- df[keep, , drop=FALSE]
  fv <- fv_full[keep]
  N<-nrow(df); sid<-df$sentence_id
  counts <- design_counts(sid)
  J <- unname(counts[["n_sentence_ids"]])
  G <- unname(counts[["n_inference_clusters"]])
  assert_contrastive_fold_binding(fv, sid)
  f<-function(rhs) as.formula(paste(y,"~",CTRL,rhs,"+",RE))
  fitkf<-function(nm,form){p<-file.path(
      CACHE,
      variant_filename(sprintf("rq3kf_v3_%s_%s.rds",tag,nm),
                       exclude_contrastive)
    ); if(file.exists(p)){
      cached <- readRDS(p)
      assert_analysis_input_hashes(cached, input_hashes, p)
      stopifnot(
        length(cached$pointwise[,"elpd_kfold"]) == N,
        identical(as.integer(attr(cached,"folds")), as.integer(fv))
      )
      return(cached)
    }
    m<-brm(form,data=df,prior=pri,control=list(adapt_delta=0.95,max_treedepth=12),chains=4,iter=2000,warmup=1000,seed=42,silent=2,refresh=0)
    kf<-kfold(m,folds=fv,chains=4,iter=2000,warmup=1000,seed=42,silent=2,refresh=0)
    kf<-set_analysis_input_hashes(kf,input_hashes)
    saveRDS(kf,p);kf}
  pb<-fitkf("base",f(""))$pointwise[,"elpd_kfold"]
  res<-data.frame(outcome=tag,predictor=c("c_nmt","c_mono"),elpd_diff=NA,
                  se_cluster=NA,p=NA,n_observations=N,n_sentence_ids=J,
                  n_inference_clusters=G)
  for(k in 1:2){v<-c("c_nmt","c_mono")[k]; pt<-fitkf(v,f(paste("+",v)))$pointwise[,"elpd_kfold"]
    di<-pt-pb; ds<-cluster_delta_sums(di,sid); res[k,3:5]<-c(sum(di),sd(ds)*sqrt(G),sfp(ds))}
  res
}
r_ffd <- run_outcome(df_ffd, "log_ffd", "FFD")
r_gd  <- run_outcome(df_gd,  "log_gd",  "GD")
r_go_past <- run_outcome(df_go_past, "log_go_past", "Go-past")
r_rrt <- run_outcome(df_rrt, "log_rrt", "RRT")
res <- rbind(r_ffd, r_gd, r_go_past, r_rrt)
res$exclude_contrastive <- exclude_contrastive
res$p_holm_secondary <- NA_real_
for (current_predictor in unique(res$predictor)) {
  secondary_rows <- res$predictor == current_predictor &
    res$outcome %in% c("FFD", "GD", "Go-past")
  res$p_holm_secondary[secondary_rows] <- p.adjust(
    res$p[secondary_rows], method="holm"
  )
}
attr(res, "outcome_definitions") <- c(
  FFD="first fixation duration on a word's first encounter",
  GD="gaze duration during a word's first encounter",
  `Go-past`=paste(
    "go-past fixation time from first landing through the fixation before",
    "the first subsequent rightward crossing; structurally undefined without",
    "a crossing or when first encountered by regression-in"
  ),
  RRT="re-reading duration conditional on at least one post-first-encounter fixation"
)
res <- set_analysis_input_hashes(res, input_hashes)
saveRDS(res, file.path(
  OUT, variant_filename("rq3_kfold_elpd.rds", exclude_contrastive)
))
cat("=== RQ3 brms kfold elpd (FFD, GD, go-past, conditional RRT) ===\n")
print(format(res, digits=3))
cat("DONE\n")
