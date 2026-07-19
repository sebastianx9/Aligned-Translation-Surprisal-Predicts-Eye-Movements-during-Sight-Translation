#!/usr/bin/env Rscript

# Sequential posterior sensitivity check for restoring the seven pooled
# read/translate S003/stoplight observations to the RQ2 joint model. The added
# rows are scored with the primary model's clean-sample centring and scaling.
# Pareto k and importance ESS diagnose the posterior reweighting.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(posterior)
})

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
path <- function(...) file.path(data_dir, ...)

model <- readRDS(path("brm_cache", "rq2_joint_maximal_v2.rds"))
fix <- read.csv(path("fixation_durations_word.csv"), stringsAsFactors = FALSE)
nmt <- read.csv(path("nmt_surprisal_soft_word.csv"), stringsAsFactors = FALSE)
mono <- read.csv(path("monolingual_surprisal_word.csv"), stringsAsFactors = FALSE)
freq <- read.table(path("subtlex_us.csv"), sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, quote = "") %>%
  transmute(word_lower = tolower(trimws(Word)), log10_freq = Lg10WF)

sentence_lengths <- nmt %>%
  group_by(sentence_id) %>%
  summarise(sentence_length = max(word_index) + 1, .groups = "drop")
predictors <- nmt %>%
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
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal = surprisal_soft, mono_surprisal, log10_freq)

raw <- fix %>%
  filter(stage %in% c("translate", "read")) %>%
  left_join(predictors, by = c("sentence_id", "word_index")) %>%
  mutate(
    log_tfd = log(total_fixation_duration_ms),
    ambiguity = factor(ambiguity),
    condition = as.integer(stage == "translate")
  ) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal),
         !is.na(log10_freq))
clean <- raw %>%
  anti_join(tibble(sentence_id = "S003", word_index = 3L),
            by = c("sentence_id", "word_index"))
added <- raw %>% filter(sentence_id == "S003", word_index == 3L)

z_from <- function(x, reference) {
  (x - mean(reference)) / sd(reference)
}
added <- added %>% mutate(
  c_nmt = z_from(nmt_surprisal, clean$nmt_surprisal),
  c_mono = z_from(mono_surprisal, clean$mono_surprisal),
  c_wlen = z_from(word_length, clean$word_length),
  c_wpos = z_from(word_position, clean$word_position),
  c_freq = z_from(log10_freq, clean$log10_freq)
)
stopifnot(nrow(added) == 7L)

log_likelihood <- log_lik(
  model, newdata = added, re_formula = NULL, allow_new_levels = FALSE
)
psis_result <- loo::psis(rowSums(log_likelihood))
weights <- weights(psis_result, normalize = TRUE, log = FALSE)
importance_ess <- 1 / sum(weights^2)
pareto_k <- loo::pareto_k_values(psis_result)

draws <- as_draws_df(model)
interaction <- draws[["b_condition:c_nmt"]]
weighted_quantile <- function(x, w, probs) {
  order_index <- order(x)
  approx(
    cumsum(w[order_index]), x[order_index], xout = probs,
    method = "linear", ties = "ordered", rule = 2
  )$y
}
interval <- weighted_quantile(interaction, weights, c(.025, .5, .975))
result <- tibble(
  added_observations = nrow(added),
  pareto_k = pareto_k,
  importance_ess = importance_ess,
  posterior_draws = length(weights),
  interaction_mean = sum(weights * interaction),
  interaction_median = interval[2],
  ci_low = interval[1],
  ci_high = interval[3],
  posterior_pr_positive = sum(weights[interaction > 0])
)
write.csv(result, path("rq2_stoplight_importance_results.csv"),
          row.names = FALSE)
print(result)
