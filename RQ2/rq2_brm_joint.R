# ── RQ2: Joint Bayesian interaction model (double dissociation) ──────────────
#
# Single joint model with BOTH predictors and BOTH stage interactions:
#   M_joint: baseline + condition*c_nmt + condition*c_mono
#
# Evidence = conditional effects (stage-specific slopes), each predictor
#            estimated while controlling for the other.
#   c_nmt  read      = beta_cnmt
#   c_nmt  translate = beta_cnmt  + beta_cnmt:cond
#   c_mono read      = beta_cmono
#   c_mono translate = beta_cmono + beta_cmono:cond
#
# No loo_compare: the LOO elpd contrast is structurally diluted (each
# predictor carries signal in only one stage) and adds nothing beyond the
# conditional-effect CIs. Dropped deliberately.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"

library(brms)
library(dplyr)

# ── Data prep (identical to rq2_brm_loo.R) ───────────────────────────────────

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep="\t",
                   header=TRUE, stringsAsFactors=FALSE, quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>%
  mutate(word=tolower(trimws(word)))

sent_lengths <- nmt %>%
  group_by(sentence_id) %>%
  summarise(sent_len = max(word_index) + 1, .groups="drop")

predictors <- nmt %>%
  left_join(sent_lengths, by="sentence_id") %>%
  mutate(word_length   = nchar(word),
         word_position = word_index / (sent_len - 1),
         word_lower    = tolower(trimws(word))) %>%
  left_join(freq, by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id, word_index, mono_surprisal=surprisal_sum),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal=surprisal_soft, mono_surprisal, log10_freq)

df <- fix %>%
  filter(stage %in% c("translate","read")) %>%
  left_join(predictors, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd   = log(total_fixation_duration_ms),
         ambiguity = factor(ambiguity),
         condition = as.integer(stage == "translate")) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
df <- df %>%
  mutate(c_nmt  = z(nmt_surprisal),
         c_mono = z(mono_surprisal),
         c_wlen = z(word_length),
         c_wpos = z(word_position),
         c_freq = z(log10_freq))

cat(sprintf("Pooled N = %d  (translate: %d  read: %d)\n\n",
            nrow(df), sum(df$condition==1), sum(df$condition==0)))

# ── Model fitting ─────────────────────────────────────────────────────────────

priors <- c(
  prior(normal(0, 1), class=b),
  prior(normal(6, 1), class=Intercept),
  prior(exponential(1), class=sd),
  prior(exponential(1), class=sigma),
  prior(lkj(2), class=cor)
)

ctrl <- list(adapt_delta=0.95, max_treedepth=12)

CACHE_DIR <- file.path(DATA_DIR, "brm_cache")
dir.create(CACHE_DIR, showWarnings=FALSE)

cache <- function(name, fit_fn) {
  path <- file.path(CACHE_DIR, paste0(name, ".rds"))
  if (file.exists(path)) { cat(sprintf("Loading %s from cache...\n", name)); readRDS(path) }
  else { m <- fit_fn(); saveRDS(m, path); m }
}

# Maximal random effects: both predictors get random slopes by participant.
# If this fails to converge, fall back to (1+condition|participant).
M_joint <- cache("M_joint", function() {
  cat("Fitting M_joint (condition * c_nmt + condition * c_mono)...\n")
  brm(log_tfd ~ c_wlen + c_wpos + c_freq + ambiguity +
        condition * c_nmt + condition * c_mono +
        (1 + condition + c_nmt + c_mono | participant) +
        (1 + condition | sentence_id),
      data=df, prior=priors, control=ctrl,
      chains=4, iter=4000, warmup=2000,
      save_pars=save_pars(all=TRUE), silent=2, refresh=0)
})

# ── Diagnostics ───────────────────────────────────────────────────────────────

cat("\n── Convergence diagnostics ──────────────────────────────────────────────\n")
sm <- summary(M_joint)
cat("Max Rhat (fixed):", round(max(sm$fixed[,"Rhat"]), 4), "\n")
cat("Min Bulk_ESS (fixed):", round(min(sm$fixed[,"Bulk_ESS"]), 0), "\n")
cat("Divergences:", sum(subset(nuts_params(M_joint), Parameter=="divergent__")$Value), "\n")

# ── Interaction coefficients ──────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("Interaction coefficients (joint model)\n")
cat("══════════════════════════════════════════════════════════════════════════\n")
cat("\ncondition:c_nmt\n");  print(fixef(M_joint)["condition:c_nmt",])
cat("\ncondition:c_mono\n"); print(fixef(M_joint)["condition:c_mono",])

# ── Conditional effects: the double dissociation ──────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("Conditional effects (stage-specific slopes, each controlling for the other)\n")
cat("══════════════════════════════════════════════════════════════════════════\n")

cat("\nc_nmt  in READ       (= beta_cnmt):\n")
print(fixef(M_joint)["c_nmt",])

cat("\nc_nmt  in TRANSLATE  (= beta_cnmt + beta_cnmt:cond):\n")
print(hypothesis(M_joint, "c_nmt + condition:c_nmt = 0")$hypothesis[,c("Estimate","CI.Lower","CI.Upper","Star")])

cat("\nc_mono in READ       (= beta_cmono):\n")
print(fixef(M_joint)["c_mono",])

cat("\nc_mono in TRANSLATE  (= beta_cmono + beta_cmono:cond):\n")
print(hypothesis(M_joint, "c_mono + condition:c_mono = 0")$hypothesis[,c("Estimate","CI.Lower","CI.Upper","Star")])

cat("\nDone.\n")
