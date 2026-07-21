#!/usr/bin/env Rscript

# Locus family: sentence-grouped predictive checks of neighbouring-word
# monolingual surprisal.  These checks are separate from RQ1's seven-predictor
# family and are restricted to interior words, for which both neighbours are
# defined:
#
#   A) reading TFD    : controls | +c_prev | +c_mono | +c_next
#   B) translate GD   : controls | +c_next
#   C) translate RRT  : controls | +c_next
#   D) translate TFD  : controls + three-word mono window | +c_nmt
#
# RRT is conditional on a regression-in.  The current eye-measure file also
# contains go-past time, but go-past is part of the RQ3 outcome profile rather
# than this neighbouring-word locus check.

suppressMessages({
  library(brms)
  library(dplyr)
})
options(mc.cores = 4)

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
repo_root <- normalizePath(
  file.path(dirname(script_file), ".."), mustWork = TRUE
)
design_path <- file.path(repo_root, "R", "analysis_design.R")
source(design_path)

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
  get_arg("--exclude-contrastive", "false"),
  "--exclude-contrastive"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fix_path <- file.path(DATA_DIR, "fixation_durations_word.csv")
eye_path <- file.path(DATA_DIR, "eye_measures_word.csv")
nmt_path <- file.path(DATA_DIR, "nmt_surprisal_soft_word.csv")
mono_path <- file.path(DATA_DIR, "monolingual_surprisal_word.csv")
freq_path <- file.path(DATA_DIR, "subtlex_us.csv")
input_hashes <- analysis_input_hashes(c(
  fixation = fix_path,
  eye_measures = eye_path,
  nmt_surprisal = nmt_path,
  monolingual_surprisal = mono_path,
  frequency = freq_path,
  analysis_design = design_path
))

fix <- read.csv(fix_path, stringsAsFactors = FALSE)
eye <- read.csv(eye_path, stringsAsFactors = FALSE)
nmt <- read.csv(nmt_path, stringsAsFactors = FALSE)
mono <- read.csv(mono_path, stringsAsFactors = FALSE)
freq <- read.table(
  freq_path, sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, quote = ""
) %>%
  transmute(
    word_lower = lexical_form(Word),
    log10_freq = Lg10WF
  )

required_fix_columns <- c(
  "participant", "sentence_id", "stage", "word_index", "ambiguity",
  "total_fixation_duration_ms"
)
missing_fix_columns <- setdiff(required_fix_columns, names(fix))
if (length(missing_fix_columns)) {
  stop(
    "fixation_durations_word.csv is missing columns: ",
    paste(missing_fix_columns, collapse = ", ")
  )
}
required_eye_columns <- c(
  "participant", "sentence_id", "stage", "word_index", "ambiguity",
  "gd_ms", "rrt_ms", "reread_occurrence"
)
missing_eye_columns <- setdiff(required_eye_columns, names(eye))
if (length(missing_eye_columns)) {
  stop(
    "eye_measures_word.csv is missing columns: ",
    paste(missing_eye_columns, collapse = ", ")
  )
}

# Form the neighbouring-word predictors at the sentence-template level before
# joining them to repeated participant observations.
mono_window <- mono %>%
  select(sentence_id, word_index, mono_curr = surprisal_sum) %>%
  group_by(sentence_id) %>%
  arrange(word_index, .by_group = TRUE) %>%
  mutate(
    mono_prev = lag(mono_curr),
    mono_next = lead(mono_curr)
  ) %>%
  ungroup()

sentence_lengths <- nmt %>%
  group_by(sentence_id) %>%
  summarise(sentence_length = max(word_index) + 1L, .groups = "drop")

predictors <- nmt %>%
  left_join(sentence_lengths, by = "sentence_id") %>%
  mutate(
    word_lower = lexical_form(word),
    word_length = nchar(word_lower, type = "chars"),
    word_position = word_index / (sentence_length - 1)
  ) %>%
  left_join(freq, by = "word_lower") %>%
  left_join(mono_window, by = c("sentence_id", "word_index")) %>%
  transmute(
    sentence_id, word_index, word_length, word_position, log10_freq,
    nmt_surprisal = surprisal_soft,
    mono_curr, mono_prev, mono_next
  )

stoplight <- tibble(sentence_id = "S003", word_index = 3L)
z_score <- function(x) {
  value_sd <- sd(x, na.rm = TRUE)
  if (!is.finite(value_sd) || value_sd == 0) {
    stop("Cannot z-score a predictor with zero or non-finite SD.")
  }
  (x - mean(x, na.rm = TRUE)) / value_sd
}

# Outcome-specific complete-case cleaning and scaling are performed on the
# full primary data.  Only afterwards is S031/S032 removed for the matched
# leave-pair-out sensitivity analysis, so both variants retain the primary
# centres, SDs, and fold labels.
prepare_outcome <- function(data, outcome_column) {
  data %>%
    left_join(predictors, by = c("sentence_id", "word_index")) %>%
    mutate(ambiguity = factor(ambiguity)) %>%
    filter(
      !is.na(.data[[outcome_column]]),
      is.finite(.data[[outcome_column]]),
      .data[[outcome_column]] > 0,
      !is.na(nmt_surprisal),
      !is.na(mono_curr),
      !is.na(mono_prev),
      !is.na(mono_next),
      !is.na(log10_freq),
      !is.na(word_length),
      !is.na(word_position),
      !is.na(ambiguity)
    ) %>%
    anti_join(stoplight, by = c("sentence_id", "word_index")) %>%
    mutate(
      y = log(.data[[outcome_column]]),
      c_nmt = z_score(nmt_surprisal),
      c_mono = z_score(mono_curr),
      c_prev = z_score(mono_prev),
      c_next = z_score(mono_next),
      c_wlen = z_score(word_length),
      c_wpos = z_score(word_position),
      c_freq = z_score(log10_freq)
    )
}

d_read_full <- prepare_outcome(
  fix %>% filter(stage == "read"),
  "total_fixation_duration_ms"
)
d_translate_full <- prepare_outcome(
  fix %>% filter(stage == "translate"),
  "total_fixation_duration_ms"
)
d_gd_full <- prepare_outcome(
  eye %>% filter(stage == "translate"),
  "gd_ms"
)
d_rrt_full <- prepare_outcome(
  eye %>% filter(stage == "translate", reread_occurrence == 1),
  "rrt_ms"
)

priors <- c(
  prior(normal(0, 1), class = b),
  prior(normal(6, 1), class = Intercept),
  prior(exponential(1), class = sd),
  prior(exponential(1), class = sigma)
)
CACHE <- file.path(DATA_DIR, "brm_cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
CTRL <- "c_wlen + c_wpos + c_freq + ambiguity"
RE <- "(1 | participant) + (1 | sentence_id)"
ANALYSIS_VERSION <- "v2"
N_SIGN_FLIPS <- 10000L
ANALYSIS_SEED <- 42L

sign_flip_p <- function(cluster_delta, n = N_SIGN_FLIPS) {
  observed <- sum(cluster_delta)
  set.seed(ANALYSIS_SEED)
  permuted <- replicate(
    n,
    sum(
      cluster_delta *
        sample(c(-1, 1), length(cluster_delta), replace = TRUE)
    )
  )
  (1 + sum(permuted >= observed)) / (n + 1)
}

make_formula <- function(extra_terms) {
  rhs <- c(CTRL, if (nzchar(extra_terms)) extra_terms, RE)
  as.formula(paste("y ~", paste(rhs, collapse = " + ")))
}

run_family <- function(data_full, family, outcome, model_terms) {
  fold_full <- make_sentence_folds(
    data_full$sentence_id, K = 10L, seed = ANALYSIS_SEED
  )
  keep <- contrastive_keep(
    data_full$sentence_id, exclude_contrastive
  )
  data <- data_full[keep, , drop = FALSE]
  folds <- fold_full[keep]
  sentence_id <- data$sentence_id
  counts <- design_counts(sentence_id)
  N <- nrow(data)
  J <- unname(counts[["n_sentence_ids"]])
  G <- unname(counts[["n_inference_clusters"]])

  if (!N || J < 10L || G < 10L) {
    stop(family, " has insufficient retained data for grouped 10-fold CV.")
  }
  assert_contrastive_fold_binding(folds, sentence_id)
  stopifnot(
    length(folds) == N,
    all(tapply(folds, sentence_id, function(x) {
      length(unique(x)) == 1L
    }))
  )

  cat(sprintf(
    "\n%s (%s): N=%d | sentence IDs=%d | inference clusters=%d | participants=%d | contrastive pair=%s\n",
    family, outcome, N, J, G, n_distinct(data$participant),
    ifelse(exclude_contrastive, "excluded", "included and paired")
  ))

  fit_kfold <- function(model_name, extra_terms) {
    cache_name <- variant_filename(
      sprintf(
        "rq_locus_kfold_%s_%s_%s.rds",
        ANALYSIS_VERSION, family, model_name
      ),
      exclude_contrastive
    )
    cache_path <- file.path(CACHE, cache_name)
    if (file.exists(cache_path)) {
      cat(sprintf("[%s/%s] loading cached kfold\n", family, model_name))
      cached <- readRDS(cache_path)
      assert_analysis_input_hashes(cached, input_hashes, cache_path)
      stopifnot(
        length(cached$pointwise[, "elpd_kfold"]) == N,
        identical(as.integer(attr(cached, "folds")), as.integer(folds))
      )
      return(cached)
    }

    cat(sprintf("[%s/%s] fitting + 10-fold refit\n", family, model_name))
    start_time <- proc.time()
    model <- brm(
      make_formula(extra_terms), data = data, prior = priors,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      chains = 4, iter = 2000, warmup = 1000,
      seed = ANALYSIS_SEED, silent = 2, refresh = 0
    )
    result <- kfold(
      model, folds = folds,
      chains = 4, iter = 2000, warmup = 1000,
      seed = ANALYSIS_SEED, silent = 2, refresh = 0
    )
    attr(result, "folds") <- as.integer(folds)
    result <- set_analysis_input_hashes(result, input_hashes)
    saveRDS(result, cache_path)
    cat(sprintf(
      "[%s/%s] completed in %.0f min\n",
      family, model_name,
      (proc.time() - start_time)[["elapsed"]] / 60
    ))
    result
  }

  if (!"base" %in% names(model_terms)) {
    stop(family, " model specification has no base model.")
  }
  base_pointwise <- fit_kfold(
    "base", model_terms[["base"]]
  )$pointwise[, "elpd_kfold"]
  candidate_names <- setdiff(names(model_terms), "base")
  candidate_pointwise <- list()
  rows <- vector("list", length(candidate_names))

  for (index in seq_along(candidate_names)) {
    model_name <- candidate_names[[index]]
    pointwise <- fit_kfold(
      model_name, model_terms[[model_name]]
    )$pointwise[, "elpd_kfold"]
    pointwise_delta <- pointwise - base_pointwise
    cluster_delta <- cluster_delta_sums(
      pointwise_delta, sentence_id
    )
    stopifnot(length(cluster_delta) == G)
    clustered_se <- sd(cluster_delta) * sqrt(G)
    naive_se <- sd(pointwise_delta) * sqrt(N)

    rows[[index]] <- data.frame(
      family = family,
      outcome = outcome,
      candidate = model_name,
      baseline_terms = model_terms[["base"]],
      candidate_terms = model_terms[[model_name]],
      elpd_diff = sum(pointwise_delta),
      se_cluster = clustered_se,
      se_naive = naive_se,
      z_cluster = sum(pointwise_delta) / clustered_se,
      p_signflip_one_sided = sign_flip_p(cluster_delta),
      per_observation = sum(pointwise_delta) / N,
      n_observations = N,
      n_sentence_ids = J,
      n_inference_clusters = G,
      n_participants = n_distinct(data$participant),
      exclude_contrastive = exclude_contrastive,
      stringsAsFactors = FALSE
    )
    candidate_pointwise[[model_name]] <- pointwise
  }

  table <- bind_rows(rows)
  cat("\n")
  print(format(table, digits = 3, nsmall = 3))
  list(
    table = table,
    pointwise_base = base_pointwise,
    pointwise_candidates = candidate_pointwise,
    sentence_id = sentence_id,
    folds = folds,
    N = N,
    J = J,
    G = G
  )
}

families <- list(
  reading = run_family(
    d_read_full, "reading", "reading TFD",
    c(
      base = "",
      s_prev = "c_prev",
      s_curr = "c_mono",
      s_next = "c_next"
    )
  ),
  gd = run_family(
    d_gd_full, "gd", "translation GD",
    c(base = "", s_next = "c_next")
  ),
  rrt = run_family(
    d_rrt_full, "rrt", "conditional translation RRT",
    c(base = "", s_next = "c_next")
  ),
  cnmt_beyond_window = run_family(
    d_translate_full, "cnmt_window", "translation TFD",
    c(
      base = "c_prev + c_mono + c_next",
      c_nmt = "c_prev + c_mono + c_next + c_nmt"
    )
  )
)

result_table <- bind_rows(lapply(families, `[[`, "table"))
# Retain nominal p-values and provide the Holm adjustment across this separate
# six-comparison locus family, leaving the reporting hierarchy explicit.
result_table$p_holm_locus_family <- p.adjust(
  result_table$p_signflip_one_sided, method = "holm"
)

result <- list(
  analysis_version = ANALYSIS_VERSION,
  table = result_table,
  families = families,
  metadata = list(
    seed = ANALYSIS_SEED,
    n_folds = 10L,
    n_sign_flips = N_SIGN_FLIPS,
    sign_flip_tail = "one-sided: candidate elpd gain greater than zero",
    sign_flip_correction = "(b + 1) / (B + 1)",
    inference_cluster = "sentence template; S031/S032 combined",
    random_effects = RE,
    controls = CTRL,
    scaling = paste(
      "outcome-specific primary complete-case sample before the optional",
      "S031/S032 exclusion"
    ),
    lexical_form = paste(
      "lower case after removing Unicode punctuation at token edges;",
      "used for both SUBTLEX matching and character length"
    ),
    stoplight_excluded = "S003 word_index 3",
    rrt_definition = "re-reading duration conditional on regression-in",
    exclude_contrastive = exclude_contrastive,
    input_hashes = input_hashes
  )
)
result <- set_analysis_input_hashes(result, input_hashes)
output_name <- variant_filename(
  sprintf("rq_locus_kfold_%s.rds", ANALYSIS_VERSION),
  exclude_contrastive
)
saveRDS(result, file.path(OUT, output_name))

cat("\nLocus-family results\n")
print(format(result_table, digits = 3, nsmall = 3))
cat("\nSaved ", output_name, "\nALL LOCUS FAMILIES DONE\n", sep = "")
