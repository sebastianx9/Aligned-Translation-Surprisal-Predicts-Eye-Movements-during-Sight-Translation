# ── RQ3: c_nmt outcome profile (vertical dot-and-error-bar) ─────────────────
# FFD, GD, go-past, and conditional RRT come from the regenerated primary
# RQ3 result;
# TFD comes from the corresponding primary RQ1 comparison. --data-dir contains
# both RDS files and --output-dir receives the PDF. Error bars are +-1.96
# sentence-clustered SE; fill denotes the nominal one-sided p-value.

suppressMessages({library(dplyr); library(ggplot2)})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value=TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}
default_data_dir <- Sys.getenv(
  "DISSERTATION_RESULTS_DIR",
  Sys.getenv("DISSERTATION_OUTPUT_DIR", ".")
)
data_dir <- normalizePath(
  get_arg("--data-dir", default_data_dir), mustWork=TRUE
)
output_dir <- get_arg(
  "--output-dir", Sys.getenv("DISSERTATION_FIGURE_DIR", data_dir)
)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

rq3_path <- file.path(data_dir, "rq3_kfold_elpd.rds")
rq1_path <- file.path(data_dir, "rq1_kfold_elpd.rds")
missing_results <- c(rq3_path, rq1_path)[
  !file.exists(c(rq3_path, rq1_path))
]
if (length(missing_results)) {
  stop("Missing regenerated result files: ",
       paste(missing_results, collapse=", "))
}
rq3 <- readRDS(rq3_path)
rq1 <- readRDS(rq1_path)

rq3_required <- c(
  "outcome", "predictor", "elpd_diff", "se_cluster", "p",
  "n_observations", "n_sentence_ids", "n_inference_clusters",
  "exclude_contrastive"
)
rq3_missing <- setdiff(rq3_required, names(rq3))
if (length(rq3_missing)) {
  stop("RQ3 result is missing columns: ",
       paste(rq3_missing, collapse=", "))
}
rq1_required <- c("res", "N", "J", "G", "exclude_contrastive")
rq1_missing <- setdiff(rq1_required, names(rq1))
if (length(rq1_missing)) {
  stop("RQ1 result is missing metadata: ",
       paste(rq1_missing, collapse=", "))
}
rq1_result_required <- c("predictor", "elpd_diff", "se_cluster", "p")
rq1_result_missing <- setdiff(rq1_result_required, names(rq1$res))
if (length(rq1_result_missing)) {
  stop("RQ1 result table is missing columns: ",
       paste(rq1_result_missing, collapse=", "))
}

if (anyNA(rq3$exclude_contrastive) || any(rq3$exclude_contrastive)) {
  stop("The primary RQ3 figure must not use the leave-pair-out result.")
}
if (!identical(rq1$exclude_contrastive, FALSE)) {
  stop("The primary RQ3 figure must use the primary RQ1 result.")
}
if (!identical(as.integer(rq1$J), 200L) ||
    !identical(as.integer(rq1$G), 199L) ||
    length(rq1$N) != 1L || !is.finite(rq1$N) || rq1$N <= 0) {
  stop("RQ1 primary N/J/G metadata are invalid.")
}

rq3_cnmt <- rq3 %>%
  filter(predictor == "c_nmt", outcome %in% c("FFD", "GD", "Go-past", "RRT"))
stopifnot(
  nrow(rq3_cnmt) == 4L,
  !anyDuplicated(rq3_cnmt$outcome),
  setequal(rq3_cnmt$outcome, c("FFD", "GD", "Go-past", "RRT")),
  all(rq3_cnmt$n_sentence_ids == 200L),
  all(rq3_cnmt$n_inference_clusters == 199L),
  all(is.finite(rq3_cnmt$n_observations)),
  all(rq3_cnmt$n_observations > 0L),
  all(rq3_cnmt$n_observations[rq3_cnmt$outcome %in% c("FFD", "GD")] ==
        rq1$N)
)
rq1_cnmt <- rq1$res %>% filter(predictor == "c_nmt")
stopifnot(nrow(rq1_cnmt) == 1L)

df <- bind_rows(
  rq3_cnmt %>% transmute(
    outcome, dllh=elpd_diff, se=se_cluster, nominal_p=p,
    n_observations, n_inference_clusters
  ),
  rq1_cnmt %>% transmute(
    outcome="TFD", dllh=elpd_diff, se=se_cluster, nominal_p=p,
    n_observations=rq1$N, n_inference_clusters=rq1$G
  )
) %>%
  mutate(
    outcome=factor(outcome, levels=c("FFD", "GD", "Go-past", "RRT", "TFD")),
    significant=nominal_p < .05,
    lo=dllh - 1.96 * se,
    hi=dllh + 1.96 * se
  ) %>%
  arrange(outcome)

p <- ggplot(df, aes(x=outcome, y=dllh)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0, linewidth=0.9,
                colour="#0072B2") +
  geom_point(aes(fill=significant), shape=21, colour="#0072B2",
             size=3.8, stroke=1.2) +
  scale_fill_manual(
    values=c(`TRUE`="#0072B2", `FALSE`="white"),
    labels=c(`TRUE`="Nominal p < .05", `FALSE`="Nominal p >= .05"),
    name=NULL
  ) +
  labs(x = NULL, y = expression("elpd"[diff])) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black", linewidth = 0.3))

figure_path <- file.path(output_dir, "rq3_gd_rrt.pdf")
ggsave(figure_path, p, width = 6, height = 4.2, device = "pdf")
cat("Saved ", figure_path, "\n", sep="")
