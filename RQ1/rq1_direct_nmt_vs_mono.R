#!/usr/bin/env Rscript

# Direct paired RQ1 comparison:
#   ELPD(M_nmt) - ELPD(M_mono)
#
# Uses pointwise held-out log predictive densities from the two models fitted
# on identical sentence-grouped folds. Differences are summed within sentence
# before the clustered SE and sign-flip test are computed.

suppressPackageStartupMessages(library(dplyr))

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
n_perm <- as.integer(get_arg("--n-perm", "1000"))
seed <- as.integer(get_arg("--seed", "42"))
stopifnot(n_perm > 0L, !is.na(seed))

path <- function(...) file.path(data_dir, ...)
mono_path <- path("monolingual_surprisal_word.csv")
nmt_path <- path("nmt_surprisal_soft_word.csv")
fix_path <- path("fixation_durations_word.csv")
freq_path <- path("subtlex_us.csv")
nmt_cache_path <- path("brm_cache", "rq1kf_c_nmt.rds")
mono_cache_path <- path("brm_cache", "rq1kf_c_mono.rds")

required <- c(mono_path, nmt_path, fix_path, freq_path,
              nmt_cache_path, mono_cache_path)
missing_files <- required[!file.exists(required)]
if (length(missing_files)) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

# Prevent reuse of the pre-alignment-fix monolingual model.
if (file.info(mono_cache_path)$mtime < file.info(mono_path)$mtime) {
  stop("The c_mono cache predates the current monolingual-surprisal CSV.")
}

kf_nmt <- readRDS(nmt_cache_path)
kf_mono <- readRDS(mono_cache_path)
for (item in list(kf_nmt, kf_mono)) {
  stopifnot(inherits(item, "kfold"), "elpd_kfold" %in% colnames(item$pointwise))
}

pw_nmt <- kf_nmt$pointwise[, "elpd_kfold"]
pw_mono <- kf_mono$pointwise[, "elpd_kfold"]
folds_nmt <- attr(kf_nmt, "folds")
folds_mono <- attr(kf_mono, "folds")
stopifnot(
  length(pw_nmt) == length(pw_mono),
  identical(folds_nmt, folds_mono),
  identical(attr(kf_nmt, "yhash"), attr(kf_mono, "yhash"))
)

fix <- read.csv(fix_path, stringsAsFactors = FALSE)
nmt <- read.csv(nmt_path, stringsAsFactors = FALSE)
mono <- read.csv(mono_path, stringsAsFactors = FALSE)
freq <- read.table(freq_path, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, quote = "") %>%
  transmute(word_lower = tolower(trimws(Word)), log10_freq = Lg10WF)
stopifnot(!anyDuplicated(mono[c("sentence_id", "word_index")]))

predictors <- nmt %>%
  mutate(word_lower = tolower(trimws(word))) %>%
  left_join(freq, by = "word_lower") %>%
  select(sentence_id, word_index, nmt_surprisal = surprisal_soft,
         log10_freq) %>%
  left_join(
    mono %>% select(sentence_id, word_index,
                    mono_surprisal = surprisal_sum),
    by = c("sentence_id", "word_index")
  )

rq1_data <- fix %>%
  left_join(predictors, by = c("sentence_id", "word_index")) %>%
  filter(
    stage == "translate",
    !is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq),
    !(sentence_id == "S003" & word_index == 3L &
        tolower(trimws(word)) == "stoplight")
  )

stopifnot(
  nrow(rq1_data) == length(pw_nmt),
  nrow(rq1_data) == 5149L,
  n_distinct(rq1_data$sentence_id) == 200L,
  all(tapply(folds_nmt, rq1_data$sentence_id,
             function(x) length(unique(x))) == 1L)
)

pointwise_delta <- pw_nmt - pw_mono
sentence_results <- tibble(
  sentence_id = rq1_data$sentence_id,
  fold = folds_nmt,
  delta_elpd = pointwise_delta
) %>%
  group_by(sentence_id) %>%
  summarise(
    fold = first(fold), n_observations = n(),
    delta_elpd = sum(delta_elpd), .groups = "drop"
  ) %>%
  arrange(sentence_id)

sentence_delta <- sentence_results$delta_elpd
n_sentences <- length(sentence_delta)
delta_elpd <- sum(pointwise_delta)
clustered_se <- sqrt(n_sentences) * sd(sentence_delta)

set.seed(seed)
permuted <- replicate(
  n_perm,
  sum(sentence_delta * sample(c(-1, 1), n_sentences, replace = TRUE))
)
p_one <- (1 + sum(permuted >= delta_elpd)) / (n_perm + 1)
p_two <- (1 + sum(abs(permuted) >= abs(delta_elpd))) / (n_perm + 1)

results <- tibble(
  contrast = "M_nmt - M_mono",
  delta_elpd = delta_elpd,
  sentence_clustered_se = clustered_se,
  ci_95_low = delta_elpd - 1.96 * clustered_se,
  ci_95_high = delta_elpd + 1.96 * clustered_se,
  p_signflip_one_sided = p_one,
  p_signflip_two_sided = p_two,
  n_observations = length(pointwise_delta),
  n_sentences = n_sentences,
  n_permutations = n_perm,
  seed = seed
)

write.csv(results, path("rq1_direct_nmt_vs_mono_results.csv"),
          row.names = FALSE)
write.csv(sentence_results,
          path("rq1_direct_nmt_vs_mono_sentence_deltas.csv"),
          row.names = FALSE)
print(results)
