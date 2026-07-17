# ── M0-M4 Unified Framework: brms + 10-fold kfold, sentence-grouped ──────────
#
# Five nested comparisons, all vs M0 baseline:
#
#   M0 = controls + condition (intercept-only baseline)
#   M1 = M0 + c_nmt + c_mono (both surprisals, no interactions)
#   M2 = M0 + c_nmt + c_mono + condition:c_nmt (single interaction 1)
#   M3 = M0 + c_nmt + c_mono + condition:c_mono (single interaction 2)
#   M4 = M0 + c_nmt + c_mono + condition:c_nmt + condition:c_mono (both interactions)
#
# All comparisons yield pooled elpd_diff with sentence-clustered SE and
# permutation-based sign-flip p-values to account for non-independence.
#
# Design choices:
#   1. RE structure (1|participant) + (1|sentence_id) identical across all 5 models
#   2. Fold assignment sentence-grouped (K=10) and shared across all models
#   3. SE clustered at sentence level, not pointwise (observation-level assumption violated)
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT      <- "/Users/sebastianx/Dissertation RQ2"
suppressMessages({library(brms); library(dplyr)})

options(mc.cores = 4)

# ── Data prep (identical to rq2_joint_maximal.R) ──────────────────────────────
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
         condition=as.integer(stage=="translate")) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L), by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
                    c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))

N <- nrow(df)
cat(sprintf("Pooled: N = %d (translate %d, read %d) | sentences = %d | participants = %d\n\n",
            N, sum(df$condition==1), sum(df$condition==0),
            n_distinct(df$sentence_id), n_distinct(df$participant)))

# ── Shared, sentence-grouped fold assignment ──────────────────────────────────
set.seed(42)
fold_vec <- loo::kfold_split_grouped(K = 10, x = df$sentence_id)
stopifnot(length(fold_vec) == N)
stopifnot(all(tapply(fold_vec, df$sentence_id, function(f) length(unique(f))) == 1))
cat("Fold sizes (observations):\n"); print(table(fold_vec))
cat("\n")

# ── Models ────────────────────────────────────────────────────────────────────
RE   <- "(1|participant) + (1|sentence_id)"
CTRL <- "c_wlen + c_wpos + c_freq + ambiguity"
f <- function(rhs) as.formula(paste("log_tfd ~ condition +", CTRL, "+", rhs, "+", RE))

forms <- list(
  M0 = as.formula(paste("log_tfd ~ condition +", CTRL, "+", RE)),
  M1 = f("c_nmt + c_mono"),
  M2 = f("c_nmt + c_mono + condition:c_nmt"),
  M3 = f("c_nmt + c_mono + condition:c_mono"),
  M4 = f("c_nmt + c_mono + condition:c_nmt + condition:c_mono")
)

priors <- c(prior(normal(0, 1), class=b),
            prior(normal(6, 1), class=Intercept),
            prior(exponential(1), class=sd),
            prior(exponential(1), class=sigma))

CACHE <- file.path(DATA_DIR, "brm_cache")
dir.create(CACHE, showWarnings=FALSE)

ptw <- list()
kfs <- list()
for (nm in names(forms)) {
  path <- file.path(CACHE, sprintf("kfold_m_unified_%s.rds", nm))
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

# ── All contrasts vs M0 ───────────────────────────────────────────────────────
sid <- df$sentence_id
S   <- n_distinct(sid)

sign_flip_p <- function(d_s, n_perm = 1000) {
  obs <- sum(d_s); set.seed(42)
  mean(replicate(n_perm, sum(d_s * sample(c(-1,1), length(d_s), replace=TRUE))) >= obs)
}

report <- function(nm, question, d_i) {
  d_s <- tapply(d_i, sid, sum)
  elpd_diff  <- sum(d_i)
  se_cluster <- sd(d_s) * sqrt(S)
  se_naive   <- sd(d_i) * sqrt(N)
  p          <- sign_flip_p(d_s)
  cat(sprintf("\n%s\n  %s\n", nm, question))
  cat(sprintf("  elpd_diff = %+8.2f\n", elpd_diff))
  cat(sprintf("    SE (sentence-clustered) = %6.2f   -> z = %+5.2f   p_signflip = %.4f%s\n",
              se_cluster, elpd_diff/se_cluster, p, ifelse(p < .05, "  *", "")))
  cat(sprintf("    SE (loo_compare default) = %6.2f   -> z = %+5.2f   [ANTI-CONSERVATIVE]\n",
              se_naive, elpd_diff/se_naive))
  cat(sprintf("    inflation factor SE_cluster / SE_naive = %.2fx\n", se_cluster/se_naive))
  cat(sprintf("  per-word Delta elpd = elpd_diff / N = %+.5f\n", elpd_diff/N))
  invisible(list(d_i=d_i, d_s=d_s, elpd_diff=elpd_diff,
                 se_cluster=se_cluster, se_naive=se_naive, p=p))
}

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("RQ2 M0-M4 UNIFIED FRAMEWORK — brms + sentence-grouped 10-fold, elpd\n")
cat(sprintf("RE (identical across all 5 models): %s\n", RE))
cat(sprintf("N = %d observations, S = %d sentences\n", N, S))
cat("══════════════════════════════════════════════════════════════════\n")

results <- list(
  m1_vs_m0 = report("M1 vs M0", "Both surprisals beyond controls?",
                    ptw$M1 - ptw$M0),
  m2_vs_m0 = report("M2 vs M0", "M0 + surprisals + c_nmt interaction?",
                    ptw$M2 - ptw$M0),
  m3_vs_m0 = report("M3 vs M0", "M0 + surprisals + c_mono interaction?",
                    ptw$M3 - ptw$M0),
  m4_vs_m0 = report("M4 vs M0", "M0 + surprisals + both interactions?",
                    ptw$M4 - ptw$M0),
  m2_vs_m1 = report("M2 vs M1", "Marginal gain: c_nmt interaction alone?",
                    ptw$M2 - ptw$M1),
  m3_vs_m1 = report("M3 vs M1", "Marginal gain: c_mono interaction alone?",
                    ptw$M3 - ptw$M1),
  m4_vs_m1 = report("M4 vs M1", "Marginal gain: both interactions together?",
                    ptw$M4 - ptw$M1)
)

# ── Cross-check: official loo_compare() ───────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("Official loo_compare() cross-check (se_diff = default, not used for inference)\n")
cat("══════════════════════════════════════════════════════════════════\n")
for (nm in c("M1", "M2", "M3", "M4")) {
  lc <- loo_compare(list(M0 = kfs$M0, Target = kfs[[nm]]))
  cat(sprintf("\n%s vs M0:\n", nm))
  print(round(lc[, c("elpd_diff", "se_diff")], 3))
}

saveRDS(list(pointwise=ptw, sentence_id=sid, results=results, N=N, S=S, forms=forms),
        file.path(OUT, "rq2_kfold_m0_m4_unified.rds"))
cat("\n\nSaved rq2_kfold_m0_m4_unified.rds\n")
