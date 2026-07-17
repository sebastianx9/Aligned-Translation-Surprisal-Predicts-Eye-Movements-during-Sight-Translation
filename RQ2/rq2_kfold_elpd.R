# ── RQ2 nested ladder, Bayesian elpd version (PROTOTYPE) ─────────────────────
#
# Same three nested comparisons as rq2_nested_ladder.R, but scored with
# brms + kfold() instead of lmer + plug-in Gaussian density.
#
#   M1 = controls + cond + c_mono + cond:c_mono
#   M2 = M1 + c_nmt                  -> c_nmt beyond c_mono
#   M3 = M2 + cond:c_nmt             -> c_nmt's effect is task-modulated
#   N1 = controls + cond + c_nmt + cond:c_nmt   (N2 == M3)
#
# Three deliberate choices, each stated in Methods:
#
# 1. RE structure identical across every compared pair: (1|participant) +
#    (1|sentence_id). A random slope present in the target but not the
#    baseline would be absorbed into elpd_diff, so the gain could no longer
#    be attributed to the predictor's fixed effect.
#
# 2. Folds are grouped by sentence (K=10, matching Lim et al. 2024), and the
#    SAME fold assignment is passed to all four models, so every comparison
#    is scored on identical held-out rows.
#
# 3. Uncertainty is CLUSTERED AT THE SENTENCE LEVEL, not taken from
#    loo_compare()'s default se_diff. The default computes sd(d_i)*sqrt(N)
#    over pointwise differences, which assumes the N word-level observations
#    are independent. They are not: words are nested in sentences and in
#    participants. The script reports both so the inflation is visible.
#
#    Draft Methods wording (PENDING experiment confirmation; NOT yet in main.tex):
#    "We report elpd differences from sentence-grouped 10-fold cross-validation.
#    Standard errors are clustered at the sentence level rather than taken from
#    the default pointwise computation, which assumes independent observations."
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT      <- "/Users/sebastianx/Dissertation RQ2"
suppressMessages({library(brms); library(dplyr)})

options(mc.cores = 4)

# ── Data prep (identical to rq2_nested_ladder.R / rq2_brm_joint.R) ───────────
fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
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
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal=surprisal_soft, mono_surprisal, log10_freq)

df <- fix %>% filter(stage %in% c("translate","read")) %>%
  left_join(pred, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity),
         cond = as.integer(stage == "translate")) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L), by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
                    c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))

N <- nrow(df)
cat(sprintf("Pooled: N = %d (translate %d, read %d) | sentences = %d | participants = %d\n\n",
            N, sum(df$cond==1), sum(df$cond==0),
            n_distinct(df$sentence_id), n_distinct(df$participant)))

# ── Shared, sentence-grouped fold assignment ─────────────────────────────────
set.seed(42)
fold_vec <- loo::kfold_split_grouped(K = 10, x = df$sentence_id)
stopifnot(length(fold_vec) == N)
# sanity: no sentence straddles two folds
stopifnot(all(tapply(fold_vec, df$sentence_id, function(f) length(unique(f))) == 1))
cat("Fold sizes (observations):\n"); print(table(fold_vec))
cat("\n")

# ── Models ───────────────────────────────────────────────────────────────────
RE   <- "(1|participant) + (1|sentence_id)"
CTRL <- "c_wlen + c_wpos + c_freq + ambiguity"
f <- function(rhs) as.formula(paste("log_tfd ~", CTRL, "+", rhs, "+", RE))

forms <- list(
  M1 = f("cond + c_mono + cond:c_mono"),
  M2 = f("cond + c_mono + cond:c_mono + c_nmt"),
  M3 = f("cond + c_mono + cond:c_mono + c_nmt + cond:c_nmt"),
  N1 = f("cond + c_nmt + cond:c_nmt")
)

priors <- c(prior(normal(0, 1), class=b),
            prior(normal(6, 1), class=Intercept),
            prior(exponential(1), class=sd),
            prior(exponential(1), class=sigma))

CACHE <- file.path(DATA_DIR, "brm_cache")
dir.create(CACHE, showWarnings=FALSE)

ptw <- list()   # pointwise elpd_kfold, one vector of length N per model
kfs <- list()   # full kfold objects, kept for the official loo_compare() cross-check
for (nm in names(forms)) {
  path <- file.path(CACHE, sprintf("kfold_%s.rds", nm))
  if (file.exists(path)) {
    cat(sprintf("[%s] loading cached kfold\n", nm))
    kf <- readRDS(path)
  } else {
    cat(sprintf("[%s] fitting + 10-fold refit ...\n", nm))
    t0 <- proc.time()
    m <- brm(forms[[nm]], data=df, prior=priors,
             control=list(adapt_delta=0.95, max_treedepth=12),
             chains=4, iter=2000, warmup=1000, silent=2, refresh=0)
    kf <- kfold(m, folds = fold_vec, chains = 4, iter = 2000, warmup = 1000,
                silent = 2, refresh = 0)
    saveRDS(kf, path)
    cat(sprintf("[%s] done in %.0f min\n", nm, (proc.time()-t0)["elapsed"]/60))
  }
  ptw[[nm]] <- kf$pointwise[, "elpd_kfold"]
  kfs[[nm]] <- kf
}

# ── The three rungs ──────────────────────────────────────────────────────────
sid <- df$sentence_id
S   <- n_distinct(sid)

sign_flip_p <- function(d_s, n_perm = 10000) {
  obs <- sum(d_s); set.seed(42)
  mean(replicate(n_perm, sum(d_s * sample(c(-1,1), length(d_s), replace=TRUE))) >= obs)
}

report <- function(nm, question, d_i) {
  d_s <- tapply(d_i, sid, sum)                 # one number per sentence (a cluster)
  elpd_diff  <- sum(d_i)
  se_cluster <- sd(d_s) * sqrt(S)              # sentence-clustered
  se_naive   <- sd(d_i) * sqrt(N)              # what loo_compare() would print
  p          <- sign_flip_p(d_s)
  cat(sprintf("\n%s\n  %s\n", nm, question))
  cat(sprintf("  elpd_diff = %+8.2f\n", elpd_diff))
  cat(sprintf("    SE (sentence-clustered) = %6.2f   -> z = %+5.2f   p_signflip = %.4f%s\n",
              se_cluster, elpd_diff/se_cluster, p, ifelse(p < .05, "  *", "")))
  cat(sprintf("    SE (loo_compare default) = %6.2f   -> z = %+5.2f   [ANTI-CONSERVATIVE]\n",
              se_naive, elpd_diff/se_naive))
  cat(sprintf("    inflation factor SE_cluster / SE_naive = %.2fx\n", se_cluster/se_naive))
  cat(sprintf("  per-word Delta llh = elpd_diff / N = %+.5f\n", elpd_diff/N))
  invisible(list(d_i=d_i, d_s=d_s, elpd_diff=elpd_diff,
                 se_cluster=se_cluster, se_naive=se_naive, p=p))
}

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("RQ2 NESTED LADDER — brms + sentence-grouped 10-fold, elpd\n")
cat(sprintf("RE (identical across every pair): %s\n", RE))
cat(sprintf("N = %d observations, S = %d sentences\n", N, S))
cat("══════════════════════════════════════════════════════════════════\n")

res <- list(
  beyond  = report("elpd_diff(M2 - M1)", "Does c_nmt add predictive power beyond c_mono?",
                   ptw$M2 - ptw$M1),
  modul   = report("elpd_diff(M3 - M2)", "Is c_nmt's effect task-modulated?",
                   ptw$M3 - ptw$M2),
  reverse = report("elpd_diff(M3 - N1)", "Does c_mono add anything beyond c_nmt?",
                   ptw$M3 - ptw$N1)
)

# ── Cross-check: the official brms/loo loo_compare() output ──────────────────
# NOT used for inference. Printed to (a) confirm our hand-computed elpd_diff
# matches brms/loo exactly, and (b) expose the default se_diff — the
# observation-level SE we replace with the sentence-clustered se_cluster.
# se_naive in report() above should match loo_compare()'s se_diff here.
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("Official loo_compare() cross-check (se_diff = the default we replace)\n")
cat("══════════════════════════════════════════════════════════════════\n")
official <- function(nm, kf_base, kf_target) {
  lc <- loo_compare(list(base = kf_base, target = kf_target))
  cat(sprintf("\n%s\n", nm))
  print(round(lc[, c("elpd_diff", "se_diff")], 3))
  cat(sprintf("  loo_compare se_diff = %.3f   (cf. se_naive above)\n",
              max(lc[, "se_diff"])))   # best model's row is 0; the other is the diff SE
}
official("M2 vs M1", kfs$M1, kfs$M2)
official("M3 vs M2", kfs$M2, kfs$M3)
official("M3 vs N1", kfs$N1, kfs$M3)

saveRDS(list(pointwise=ptw, sentence_id=sid, res=res, N=N, S=S),
        file.path(OUT, "rq2_kfold_elpd.rds"))
cat("\nSaved rq2_kfold_elpd.rds\n")
