#!/usr/bin/env Rscript

# Translation-stage nested comparison:
#   controls + c_mono + c_nmt  versus  controls + c_mono.
# Uses the same complete-case sample, grouped folds, priors, and cache naming as
# RQ1/rq1_kfold_elpd.R, so the c_mono model can be reused after RQ1 finishes.

suppressMessages({library(brms); library(dplyr)})
options(mc.cores = 4, warn = 1)

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
OUT <- get_arg(
  "--output-dir", Sys.getenv("DISSERTATION_OUTPUT_DIR", DATA_DIR)
)
exclude_contrastive <- parse_bool(
  get_arg("--exclude-contrastive", "false"), "--exclude-contrastive"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fix_path <- file.path(DATA_DIR, "fixation_durations_word.csv")
nmt_path <- file.path(DATA_DIR, "nmt_surprisal_soft_word.csv")
mono_path <- file.path(DATA_DIR, "monolingual_surprisal_word.csv")
attn_path <- file.path(DATA_DIR, "attention_features_6_norm.csv")
freq_path <- file.path(DATA_DIR, "subtlex_us.csv")
input_hashes <- analysis_input_hashes(c(
  fixation=fix_path, nmt_surprisal=nmt_path,
  monolingual_surprisal=mono_path, attention_features=attn_path,
  frequency=freq_path,
  analysis_design=file.path(repo_root, "R", "analysis_design.R")
))

fix <- read.csv(fix_path, stringsAsFactors = FALSE)
nmt <- read.csv(nmt_path, stringsAsFactors = FALSE)
mono <- read.csv(mono_path, stringsAsFactors = FALSE)
attn <- read.csv(attn_path, stringsAsFactors = FALSE)
freq <- read.table(
  freq_path, sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, quote = ""
) %>%
  transmute(word_lower = lexical_form(Word), log10_freq = Lg10WF)

sentence_lengths <- nmt %>%
  group_by(sentence_id) %>%
  summarise(sentence_length = max(word_index) + 1, .groups = "drop")
predictors <- nmt %>%
  left_join(sentence_lengths, by = "sentence_id") %>%
  mutate(
    word_lower = lexical_form(word),
    word_length = nchar(word_lower, type = "chars"),
    word_position = word_index / (sentence_length - 1)
  ) %>%
  left_join(freq, by = "word_lower") %>%
  left_join(
    mono %>% select(
      sentence_id, word_index, mono_surprisal = surprisal_sum
    ),
    by = c("sentence_id", "word_index")
  ) %>%
  left_join(
    attn %>% select(
      sentence_id, word_index, H_e = attn_entropy,
      f_e = attn_context, f_eos = attn_eos,
      f_recv = attn_recv, f_cross = attn_cross
    ),
    by = c("sentence_id", "word_index")
  ) %>%
  select(
    sentence_id, word_index, word_length, word_position, log10_freq,
    nmt_surprisal = surprisal_soft, mono_surprisal,
    H_e, f_e, f_eos, f_recv, f_cross
  )

df <- fix %>%
  filter(stage == "translate") %>%
  left_join(predictors, by = c("sentence_id", "word_index")) %>%
  mutate(
    log_tfd = log(total_fixation_duration_ms),
    ambiguity = factor(ambiguity)
  ) %>%
  filter(
    !is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq),
    !is.na(H_e), !is.na(f_e), !is.na(f_eos), !is.na(f_recv),
    !is.na(f_cross)
  ) %>%
  anti_join(tibble(sentence_id = "S003", word_index = 3L),
            by = c("sentence_id", "word_index"))

z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
df <- df %>% mutate(
  c_nmt = z(nmt_surprisal), c_mono = z(mono_surprisal),
  c_wlen = z(word_length), c_wpos = z(word_position),
  c_freq = z(log10_freq)
)
fold_vec_full <- make_sentence_folds(df$sentence_id, K = 10, seed = 42)
keep <- contrastive_keep(df$sentence_id, exclude_contrastive)
df <- df[keep, , drop = FALSE]
fold_vec <- fold_vec_full[keep]
N <- nrow(df)
sid <- df$sentence_id
counts <- design_counts(sid)
J <- unname(counts[["n_sentence_ids"]])
G <- unname(counts[["n_inference_clusters"]])
cat(sprintf(
  "Translate: N=%d | sentence IDs=%d | inference clusters=%d | participants=%d | contrastive pair=%s\n",
  N, J, G, n_distinct(df$participant),
  ifelse(exclude_contrastive, "excluded", "included and paired")
))

stopifnot(
  N > 0L,
  all(tapply(fold_vec, sid, function(x) length(unique(x))) == 1L)
)
assert_contrastive_fold_binding(fold_vec, sid)

priors <- c(
  prior(normal(0, 1), class = b),
  prior(normal(6, 1), class = Intercept),
  prior(exponential(1), class = sd),
  prior(exponential(1), class = sigma)
)
CACHE <- file.path(DATA_DIR, "brm_cache")
dir.create(CACHE, showWarnings = FALSE)
CTRL <- "c_wlen + c_wpos + c_freq + ambiguity"
RE <- "(1 | participant) + (1 | sentence_id)"

fit_kfold <- function(name, formula) {
  cache_path <- file.path(
    CACHE,
    variant_filename(sprintf("rq1kf_v4_%s.rds", name),
                     exclude_contrastive)
  )
  if (file.exists(cache_path)) {
    cat(sprintf("[%s] loading cached kfold\n", name))
    result <- readRDS(cache_path)
    assert_analysis_input_hashes(result, input_hashes, cache_path)
    stopifnot(isTRUE(all.equal(
      as.numeric(attr(result, "folds")), as.numeric(fold_vec)
    )))
    return(result)
  }
  cat(sprintf("[%s] fitting + 10-fold refit\n", name))
  model <- brm(
    formula, data = df, prior = priors,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    chains = 4, iter = 4000, warmup = 2000, seed = 42,
    silent = 2, refresh = 0
  )
  result <- kfold(
    model, folds = fold_vec, chains = 4, iter = 4000, warmup = 2000,
    seed = 42, silent = 2, refresh = 0
  )
  result <- set_analysis_input_hashes(result, input_hashes)
  saveRDS(result, cache_path)
  result
}

mono_formula <- as.formula(paste(
  "log_tfd ~", CTRL, "+ c_mono +", RE
))
both_formula <- as.formula(paste(
  "log_tfd ~", CTRL, "+ c_mono + c_nmt +", RE
))
ptw_mono <- fit_kfold("c_mono", mono_formula)$pointwise[, "elpd_kfold"]
ptw_both <- fit_kfold("c_mono_nmt", both_formula)$pointwise[, "elpd_kfold"]

pointwise_delta <- ptw_both - ptw_mono
sentence_delta <- cluster_delta_sums(pointwise_delta, sid)
delta_elpd <- sum(pointwise_delta)
clustered_se <- sd(sentence_delta) * sqrt(G)
set.seed(42)
permuted <- replicate(
  10000L,
  sum(sentence_delta * sample(c(-1, 1), G, replace = TRUE))
)
p_signflip <- (1 + sum(permuted >= delta_elpd)) / 10001

result <- list(
  elpd_diff = delta_elpd,
  se_cluster = clustered_se,
  p = p_signflip,
  per_word = delta_elpd / N,
  pointwise_mono = ptw_mono,
  pointwise_both = ptw_both,
  sid = sid,
  N = N,
  J = J,
  G = G,
  S = G,
  exclude_contrastive = exclude_contrastive,
  folds = fold_vec,
  input_hashes = input_hashes
)
print(data.frame(
  contrast = "c_nmt beyond controls + c_mono",
  delta_elpd = delta_elpd,
  clustered_se = clustered_se,
  p_signflip_one_sided = p_signflip,
  n_observations = N,
  n_sentence_ids = J,
  n_inference_clusters = G
))
saveRDS(
  result,
  file.path(OUT, variant_filename("rq2_beyond_kfold.rds",
                                  exclude_contrastive))
)
