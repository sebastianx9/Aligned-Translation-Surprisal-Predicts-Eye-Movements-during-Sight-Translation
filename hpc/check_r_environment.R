#!/usr/bin/env Rscript

# Fail-fast CSF check: packages, corrected input data, and a tiny Stan fit.

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
skip_stan <- tolower(get_arg("--skip-stan", "false")) == "true"

required_packages <- c("brms", "dplyr", "loo", "posterior")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nRun hpc/csf3_install_packages.sbatch, then resubmit this check."
  )
}

required_files <- c(
  "fixation_durations_word.csv", "eye_measures_word.csv",
  "nmt_surprisal_soft_word.csv", "nmt_alignment_mass_word.csv",
  "monolingual_surprisal_word.csv", "attention_features_6_norm.csv",
  "subtlex_us.csv"
)
missing_files <- required_files[
  !file.exists(file.path(data_dir, required_files))
]
if (length(missing_files)) {
  stop("Missing analysis inputs: ", paste(missing_files, collapse = ", "))
}

fix <- read.csv(
  file.path(data_dir, "fixation_durations_word.csv"),
  stringsAsFactors = FALSE
)
eye <- read.csv(
  file.path(data_dir, "eye_measures_word.csv"),
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(fix) == nrow(eye),
  nrow(fix) == 13178L,
  length(unique(fix$participant)) == 42L,
  !any(fix$sentence_id %in% c("S031", "S032")),
  !any(eye$sentence_id %in% c("S031", "S032"))
)
cat(sprintf(
  "Corrected inputs: %d rows | %d participants | read=%d | translate=%d\n",
  nrow(fix), length(unique(fix$participant)),
  sum(fix$stage == "read"), sum(fix$stage == "translate")
))

suppressPackageStartupMessages(library(dplyr))
nmt <- read.csv(
  file.path(data_dir, "nmt_surprisal_soft_word.csv"),
  stringsAsFactors = FALSE
)
mono <- read.csv(
  file.path(data_dir, "monolingual_surprisal_word.csv"),
  stringsAsFactors = FALSE
)
attn <- read.csv(
  file.path(data_dir, "attention_features_6_norm.csv"),
  stringsAsFactors = FALSE
)
freq <- read.table(
  file.path(data_dir, "subtlex_us.csv"), sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, quote = ""
) %>%
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
  )

primary <- fix %>%
  filter(stage == "translate") %>%
  left_join(predictors, by = c("sentence_id", "word_index")) %>%
  filter(
    !is.na(surprisal_soft), !is.na(mono_surprisal), !is.na(log10_freq),
    !is.na(H_e), !is.na(f_e), !is.na(f_eos), !is.na(f_recv),
    !is.na(f_cross)
  ) %>%
  anti_join(tibble(sentence_id = "S003", word_index = 3L),
            by = c("sentence_id", "word_index"))
pooled <- fix %>%
  filter(stage %in% c("translate", "read")) %>%
  left_join(predictors, by = c("sentence_id", "word_index")) %>%
  filter(
    !is.na(surprisal_soft), !is.na(mono_surprisal), !is.na(log10_freq)
  ) %>%
  anti_join(tibble(sentence_id = "S003", word_index = 3L),
            by = c("sentence_id", "word_index"))
rrt <- eye %>%
  filter(stage == "translate") %>%
  left_join(predictors, by = c("sentence_id", "word_index")) %>%
  filter(
    !is.na(surprisal_soft), !is.na(mono_surprisal), !is.na(log10_freq),
    !(sentence_id == "S003" & word_index == 3L),
    regress_in == 1L, rrt_ms > 0
  )
stopifnot(
  nrow(primary) == 5483L,
  n_distinct(primary$sentence_id) == 198L,
  nrow(pooled) == 12018L,
  sum(pooled$stage == "translate") == 5483L,
  sum(pooled$stage == "read") == 6535L,
  nrow(rrt) == 2640L
)
cat(sprintf(
  paste0(
    "Analysis samples: RQ1=%d/%d sentences | ",
    "RQ2=%d (%d translate + %d read) | conditional RRT=%d\n"
  ),
  nrow(primary), n_distinct(primary$sentence_id), nrow(pooled),
  sum(pooled$stage == "translate"), sum(pooled$stage == "read"),
  nrow(rrt)
))

cat("\nPackage versions:\n")
for (package in required_packages) {
  cat(sprintf("  %-10s %s\n", package, as.character(packageVersion(package))))
}

if (!skip_stan) {
  # This deliberately small model validates the full brms -> Stan toolchain.
  set.seed(42)
  smoke_data <- data.frame(
    y = rnorm(40), x = rnorm(40),
    group = factor(rep(seq_len(10), each = 4))
  )
  smoke_model <- brms::brm(
    y ~ x + (1 | group), data = smoke_data,
    chains = 2, cores = 2, iter = 300, warmup = 150,
    seed = 42, refresh = 0, silent = 2
  )
  stopifnot(all(is.finite(brms::fixef(smoke_model))))
  cat("\nStan smoke fit completed successfully.\n\n")
} else {
  cat("\nStan smoke fit skipped by request.\n\n")
}
print(sessionInfo())
