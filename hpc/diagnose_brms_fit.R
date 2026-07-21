#!/usr/bin/env Rscript

# Summarise convergence diagnostics from a cached brms fit without refitting it.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

fit_path <- get_arg("--fit-path")
if (is.null(fit_path) || !nzchar(fit_path)) {
  stop("Supply --fit-path=/absolute/path/to/model.rds")
}
fit_path <- normalizePath(fit_path, mustWork = TRUE)

output_dir <- get_arg("--output-dir", dirname(fit_path))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

required <- c("brms", "posterior")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "))
}

cat("Reading: ", fit_path, "\n", sep = "")
fit <- readRDS(fit_path)
if (!inherits(fit, "brmsfit")) {
  stop("The cached object is not a brmsfit: ", fit_path)
}

draws <- posterior::as_draws_array(fit)
diagnostics <- posterior::summarise_draws(
  draws, "rhat", "ess_bulk", "ess_tail"
)
diagnostics <- as.data.frame(diagnostics, stringsAsFactors = FALSE)
diagnostics <- diagnostics[
  is.finite(diagnostics$rhat) |
    is.finite(diagnostics$ess_bulk) |
    is.finite(diagnostics$ess_tail),
  c("variable", "rhat", "ess_bulk", "ess_tail"),
  drop = FALSE
]

stem <- tools::file_path_sans_ext(basename(fit_path))
csv_path <- file.path(output_dir, paste0(stem, "_sampling_diagnostics.csv"))
write.csv(diagnostics, csv_path, row.names = FALSE)

n_draws <- posterior::ndraws(draws)
finite_rhat <- diagnostics$rhat[is.finite(diagnostics$rhat)]
finite_bulk <- diagnostics$ess_bulk[is.finite(diagnostics$ess_bulk)]
finite_tail <- diagnostics$ess_tail[is.finite(diagnostics$ess_tail)]

cat(sprintf(
  paste0(
    "Posterior draws: %d\nParameters diagnosed: %d\n",
    "Maximum R-hat: %.4f | >1.01: %d | >1.05: %d\n",
    "Minimum bulk ESS: %.1f | below 400: %d\n",
    "Minimum tail ESS: %.1f | below 400: %d\n"
  ),
  n_draws, nrow(diagnostics),
  max(finite_rhat), sum(finite_rhat > 1.01), sum(finite_rhat > 1.05),
  min(finite_bulk), sum(finite_bulk < 400),
  min(finite_tail), sum(finite_tail < 400)
))

print_rows <- function(title, data, order_by, decreasing = FALSE, n = 15L) {
  cat("\n", title, "\n", sep = "")
  keep <- is.finite(data[[order_by]])
  out <- data[keep, , drop = FALSE]
  out <- out[order(out[[order_by]], decreasing = decreasing), , drop = FALSE]
  print(utils::head(out, n), row.names = FALSE, digits = 4)
}

print_rows("Largest R-hat values:", diagnostics, "rhat", decreasing = TRUE)
print_rows("Lowest bulk ESS values:", diagnostics, "ess_bulk")
print_rows("Lowest tail ESS values:", diagnostics, "ess_tail")

key_parameters <- diagnostics[
  grepl("^(b_|Intercept$|sigma$|sd_|cor_)", diagnostics$variable),
  , drop = FALSE
]
print_rows(
  "Largest R-hat values among fixed effects and distribution/group parameters:",
  key_parameters, "rhat", decreasing = TRUE, n = 30L
)

cat("\nFull diagnostics written to: ", csv_path, "\n", sep = "")
