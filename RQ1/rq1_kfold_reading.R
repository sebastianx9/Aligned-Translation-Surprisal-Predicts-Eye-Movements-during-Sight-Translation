# ── RQ1 brms kfold elpd — 7 predictors vs baseline, READING stage ─────────────
# Reading-stage companion to rq1_kfold_elpd.R (translate). Identical convention:
# sentence-grouped 10-fold, per-model pointwise elpd, sentence-clustered SE,
# sentence-level sign-flip permutation, within-READING z-scoring.
# WITHIN-STAGE ONLY — not a cross-stage test; specificity stays the RQ2 coeff.
# lmer 200-fold reading values preserved in rq1_loo_reading_authoritative.rds.
# ─────────────────────────────────────────────────────────────────────────────
DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT      <- "/Users/sebastianx/Dissertation RQ1"
suppressMessages({library(brms); library(dplyr)})
options(mc.cores = 4)

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
attn <- read.csv(file.path(DATA_DIR, "attention_features_6_norm.csv"),  stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep="\t",
                   header=TRUE, stringsAsFactors=FALSE, quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))

sl <- nmt %>% group_by(sentence_id) %>% summarise(sent_len=max(word_index)+1, .groups="drop")
pred <- nmt %>% left_join(sl, by="sentence_id") %>%
  mutate(word_length=nchar(word), word_position=word_index/(sent_len-1),
         word_lower=tolower(trimws(word))) %>%
  left_join(freq, by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id, word_index, mono_surprisal=surprisal_sum),
            by=c("sentence_id","word_index")) %>%
  left_join(attn %>% select(sentence_id, word_index, H_e=attn_entropy, f_e=attn_context,
                            f_eos=attn_eos, f_recv=attn_recv, f_cross=attn_cross),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position, log10_freq,
         nmt_surprisal=surprisal_soft, mono_surprisal, H_e, f_e, f_eos, f_recv, f_cross)

df <- fix %>% filter(stage == "read") %>%
  left_join(pred, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity)) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq),
         !is.na(H_e), !is.na(f_e), !is.na(f_eos), !is.na(f_recv), !is.na(f_cross)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L), by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
                    c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq),
                    c_He=z(H_e), c_fe=z(f_e), c_feos=z(f_eos),
                    c_frecv=z(f_recv), c_fcross=z(f_cross))
N <- nrow(df); sid <- df$sentence_id; S <- n_distinct(sid)
cat(sprintf("Reading: N=%d  sentences=%d  participants=%d\n\n", N, S, n_distinct(df$participant)))

set.seed(42)
fold_vec <- loo::kfold_split_grouped(K=10, x=df$sentence_id)
stopifnot(all(tapply(fold_vec, df$sentence_id, function(f) length(unique(f)))==1))

priors <- c(prior(normal(0,1), class=b), prior(normal(6,1), class=Intercept),
            prior(exponential(1), class=sd), prior(exponential(1), class=sigma))
CACHE <- file.path(DATA_DIR, "brm_cache"); dir.create(CACHE, showWarnings=FALSE)
CTRL <- "c_wlen + c_wpos + c_freq + ambiguity"
RE   <- "(1|participant) + (1|sentence_id)"
f <- function(rhs) as.formula(paste("log_tfd ~", CTRL, rhs, "+", RE))

pred_vars   <- c("c_nmt","c_mono","c_He","c_fe","c_feos","c_frecv","c_fcross")
pred_labels <- c("c_nmt","c_mono","H_e","f_e","f_eos","f_recv","f_cross")

fit_kfold <- function(name, formula) {
  path <- file.path(CACHE, sprintf("rq1kfR_%s.rds", name))
  if (file.exists(path)) { cat(sprintf("[%s] cached\n", name)); return(readRDS(path)) }
  cat(sprintf("[%s] fitting + 10-fold ...\n", name)); t0 <- proc.time()
  m  <- brm(formula, data=df, prior=priors, control=list(adapt_delta=0.95, max_treedepth=12),
            chains=4, iter=2000, warmup=1000, silent=2, refresh=0)
  kf <- kfold(m, folds=fold_vec, chains=4, iter=2000, warmup=1000, silent=2, refresh=0)
  saveRDS(kf, path); cat(sprintf("[%s] %.0f min\n", name, (proc.time()-t0)["elapsed"]/60)); kf
}

ptw_base <- fit_kfold("base", f(""))$pointwise[, "elpd_kfold"]
sign_flip_p <- function(d_s, n=10000){ obs<-sum(d_s); set.seed(42)
  mean(replicate(n, sum(d_s*sample(c(-1,1),length(d_s),replace=TRUE))) >= obs) }

res <- data.frame(predictor=pred_labels, elpd_diff=NA, se_cluster=NA, se_naive=NA,
                  z=NA, p=NA, per_word=NA)
for (k in seq_along(pred_vars)) {
  ptw <- fit_kfold(pred_vars[k], f(paste("+", pred_vars[k])))$pointwise[, "elpd_kfold"]
  d_i <- ptw - ptw_base; d_s <- tapply(d_i, sid, sum)
  res[k, -1] <- c(sum(d_i), sd(d_s)*sqrt(S), sd(d_i)*sqrt(N),
                  sum(d_i)/(sd(d_s)*sqrt(S)), sign_flip_p(d_s), sum(d_i)/N)
}
res$p_holm <- p.adjust(res$p, method="holm")

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("RQ1 READING — brms sentence-grouped 10-fold elpd, sentence-clustered SE\n")
cat(sprintf("N=%d obs, S=%d sentences (WITHIN-STAGE ONLY)\n", N, S))
cat("══════════════════════════════════════════════════════════════════\n")
print(format(res, digits=3, nsmall=3))
saveRDS(list(res=res, pointwise_base=ptw_base, sid=sid, N=N, S=S),
        file.path(OUT, "rq1_kfold_reading.rds"))
cat("\nSaved rq1_kfold_reading.rds\nDONE\n")
