#!/usr/bin/env Rscript

# RQ1 sensitivity checks:
#   (1) c_nmt beyond effective soft-alignment mass;
#   (2) c_nmt with the S003/stoplight observations restored.

suppressPackageStartupMessages({library(brms); library(dplyr)})
options(mc.cores = 4)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}
data_dir <- normalizePath(
  get_arg("--data-dir", Sys.getenv("DISSERTATION_DATA_DIR", ".")),
  mustWork = TRUE
)
output_dir <- get_arg(
  "--output-dir", Sys.getenv("DISSERTATION_OUTPUT_DIR", data_dir)
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
path <- function(...) file.path(data_dir, ...)

fix <- read.csv(path("fixation_durations_word.csv"), stringsAsFactors = FALSE)
nmt <- read.csv(path("nmt_surprisal_soft_word.csv"), stringsAsFactors = FALSE)
mono <- read.csv(path("monolingual_surprisal_word.csv"), stringsAsFactors = FALSE)
mass <- read.csv(path("nmt_alignment_mass_word.csv"), stringsAsFactors = FALSE)
attn <- read.csv(path("attention_features_6_norm.csv"), stringsAsFactors = FALSE)
freq <- read.table(path("subtlex_us.csv"), sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, quote = "") %>%
  transmute(word_lower = tolower(trimws(Word)), log10_freq = Lg10WF)

# The mass extractor must reproduce the primary c_nmt values exactly at the
# stored six-decimal precision.
mass_check <- nmt %>%
  select(sentence_id, word_index, surprisal_soft) %>%
  inner_join(
    mass %>% select(sentence_id, word_index, surprisal_soft_recomputed),
    by = c("sentence_id", "word_index")
  )
stopifnot(
  nrow(mass_check) == nrow(nmt),
  max(abs(mass_check$surprisal_soft -
            mass_check$surprisal_soft_recomputed)) == 0
)

sentence_lengths <- nmt %>%
  group_by(sentence_id) %>%
  summarise(sentence_length = max(word_index) + 1, .groups = "drop")

predictors <- nmt %>%
  left_join(
    mass %>% select(sentence_id, word_index, alignment_mass,
                    surprisal_per_mass),
    by = c("sentence_id", "word_index")
  ) %>%
  left_join(sentence_lengths, by = "sentence_id") %>%
  mutate(
    word_length = nchar(word),
    word_position = word_index / (sentence_length - 1),
    word_lower = tolower(trimws(word))
  ) %>%
  left_join(freq, by = "word_lower") %>%
  left_join(
    mono %>% select(sentence_id, word_index,
                    mono_surprisal = surprisal_sum),
    by = c("sentence_id", "word_index")
  ) %>%
  left_join(
    attn %>% select(sentence_id, word_index,
                    H_e = attn_entropy, f_e = attn_context,
                    f_eos = attn_eos, f_recv = attn_recv,
                    f_cross = attn_cross),
    by = c("sentence_id", "word_index")
  ) %>%
  select(
    sentence_id, word_index, word_length, word_position, log10_freq,
    nmt_surprisal = surprisal_soft, mono_surprisal, alignment_mass,
    surprisal_per_mass, H_e, f_e, f_eos, f_recv, f_cross
  )

required_predictors <- c(
  "nmt_surprisal", "mono_surprisal", "alignment_mass",
  "surprisal_per_mass", "log10_freq", "H_e", "f_e", "f_eos",
  "f_recv", "f_cross"
)
translate_all <- fix %>%
  filter(stage == "translate") %>%
  left_join(predictors, by = c("sentence_id", "word_index")) %>%
  mutate(log_tfd = log(total_fixation_duration_ms),
         ambiguity = factor(ambiguity)) %>%
  filter(if_all(all_of(required_predictors), ~ !is.na(.x)))
translate_clean <- translate_all %>%
  anti_join(tibble(sentence_id = "S003", word_index = 3L),
            by = c("sentence_id", "word_index"))

z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
scale_data <- function(data) {
  data %>% mutate(
    c_nmt = z(nmt_surprisal), c_mass = z(alignment_mass),
    c_wlen = z(word_length), c_wpos = z(word_position),
    c_freq = z(log10_freq), c_mono = z(mono_surprisal)
  )
}
translate_clean <- scale_data(translate_clean)
translate_all <- scale_data(translate_all)
stopifnot(
  nrow(translate_clean) > 0L,
  nrow(translate_all) >= nrow(translate_clean),
  nrow(translate_all) - nrow(translate_clean) ==
    sum(translate_all$sentence_id == "S003" &
          translate_all$word_index == 3L)
)

priors <- c(
  prior(normal(0, 1), class = b),
  prior(normal(6, 1), class = Intercept),
  prior(exponential(1), class = sd),
  prior(exponential(1), class = sigma)
)
cache_dir <- path("brm_cache")
dir.create(cache_dir, showWarnings = FALSE)
control_terms <- "c_wlen + c_wpos + c_freq + ambiguity"
random_intercepts <- "(1 | participant) + (1 | sentence_id)"
make_formula <- function(extra) {
  as.formula(paste("log_tfd ~", control_terms, extra, "+",
                   random_intercepts))
}

make_folds <- function(data) {
  set.seed(42)
  folds <- loo::kfold_split_grouped(K = 10, x = data$sentence_id)
  stopifnot(all(tapply(folds, data$sentence_id,
                       function(x) length(unique(x))) == 1L))
  folds
}

fit_kfold <- function(name, formula, data, folds) {
  cache_path <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(cache_path)) return(readRDS(cache_path))
  model <- brm(
    formula, data = data, prior = priors,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    chains = 4, iter = 2000, warmup = 1000, seed = 42,
    silent = 2, refresh = 0
  )
  result <- kfold(
    model, folds = folds, chains = 4, iter = 2000, warmup = 1000,
    seed = 42, silent = 2, refresh = 0
  )
  saveRDS(result, cache_path)
  result
}

sign_flip <- function(sentence_delta, n_perm = 1000L) {
  observed <- sum(sentence_delta)
  set.seed(42)
  permuted <- replicate(
    n_perm,
    sum(sentence_delta * sample(c(-1, 1), length(sentence_delta), TRUE))
  )
  (1 + sum(permuted >= observed)) / (n_perm + 1)
}

compare_kfold <- function(target, baseline, sentence_id, contrast) {
  delta <- target$pointwise[, "elpd_kfold"] -
    baseline$pointwise[, "elpd_kfold"]
  sentence_delta <- tapply(delta, sentence_id, sum)
  tibble(
    contrast = contrast,
    delta_elpd = sum(delta),
    clustered_se = sqrt(length(sentence_delta)) * sd(sentence_delta),
    p_signflip_one_sided = sign_flip(sentence_delta),
    n_observations = length(delta),
    n_sentences = length(sentence_delta)
  )
}

fold_clean <- make_folds(translate_clean)
primary_cache <- path("brm_cache", "rq1kf_c_nmt.rds")
if (file.exists(primary_cache)) {
  stopifnot(isTRUE(all.equal(
    as.numeric(fold_clean),
    as.numeric(attr(readRDS(primary_cache), "folds"))
  )))
}

kf_mass <- fit_kfold(
  "rq1rob_mass_base", make_formula("+ c_mass"),
  translate_clean, fold_clean
)
kf_mass_nmt <- fit_kfold(
  "rq1rob_mass_cnmt", make_formula("+ c_mass + c_nmt"),
  translate_clean, fold_clean
)
result_mass <- compare_kfold(
  kf_mass_nmt, kf_mass, translate_clean$sentence_id,
  "c_nmt gain beyond alignment mass"
)

fold_all <- make_folds(translate_all)
kf_stop_base <- fit_kfold(
  "rq1rob_stoplight_base", make_formula(""), translate_all, fold_all
)
kf_stop_nmt <- fit_kfold(
  "rq1rob_stoplight_cnmt", make_formula("+ c_nmt"),
  translate_all, fold_all
)
result_stoplight <- compare_kfold(
  kf_stop_nmt, kf_stop_base, translate_all$sentence_id,
  "c_nmt gain with stoplight retained"
)

results <- bind_rows(result_mass, result_stoplight)
write.csv(results,
          file.path(output_dir, "rq1_mass_stoplight_predictive_results.csv"),
          row.names = FALSE)
print(results)

# Full-data coefficient checks use the same maximal random slopes as the
# corresponding RQ1 coefficient model.
fit_model <- function(name, formula, data) {
  cache_path <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(cache_path)) return(readRDS(cache_path))
  model <- brm(
    formula, data = data, prior = c(priors, prior(lkj(2), class = cor)),
    control = list(adapt_delta = 0.99, max_treedepth = 14),
    chains = 4, iter = 4000, warmup = 2000, seed = 42,
    silent = 2, refresh = 0
  )
  saveRDS(model, cache_path)
  model
}

mass_formula <- as.formula(paste(
  "log_tfd ~", control_terms, "+ c_mass + c_nmt +",
  "(1 + c_nmt | participant) + (1 + c_nmt | sentence_id)"
))
stop_formula <- as.formula(paste(
  "log_tfd ~", control_terms, "+ c_nmt +",
  "(1 + c_nmt | participant) + (1 + c_nmt | sentence_id)"
))
coef_mass <- fit_model(
  "rq1rob_coef_mass_cnmt", mass_formula, translate_clean
)
coef_stop <- fit_model(
  "rq1rob_coef_stoplight_cnmt", stop_formula, translate_all
)

extract_term <- function(model, term, analysis) {
  estimates <- fixef(model)
  tibble(
    analysis = analysis, term = term,
    estimate = estimates[term, "Estimate"],
    posterior_se = estimates[term, "Est.Error"],
    ci_low = estimates[term, "Q2.5"],
    ci_high = estimates[term, "Q97.5"]
  )
}
coefficient_results <- bind_rows(
  extract_term(coef_mass, "c_nmt", "c_nmt controlling alignment mass"),
  extract_term(coef_mass, "c_mass", "alignment mass controlling c_nmt"),
  extract_term(coef_stop, "c_nmt", "c_nmt with stoplight retained")
)
write.csv(
  coefficient_results,
  file.path(output_dir, "rq1_mass_stoplight_coefficient_results.csv"),
  row.names = FALSE
)
print(coefficient_results)
