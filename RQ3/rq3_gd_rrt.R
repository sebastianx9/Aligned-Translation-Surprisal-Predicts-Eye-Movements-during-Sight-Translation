# ── RQ3: GD vs RRT kfold predictive gains (vertical dot-and-error-bar) ─────
# Values read DIRECTLY from rq3_kfold_elpd.rds (no hardcoding — a stale
# hardcoded copy of this script silently re-plotted pre-bugfix numbers once).
# Blue = c_nmt, orange = c_mono; filled = p < .05, open = p >= .05.
# Error bars are +-1.96 sentence-clustered SE.

OUT_DIR <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(dplyr); library(ggplot2)})

res <- readRDS("/Users/sebastianx/Dissertation RQ3/rq3_kfold_elpd.rds")
df <- res %>%
  transmute(predictor = ifelse(predictor=="c_nmt","c[nmt]","c[mono]"),
            outcome, dllh = elpd_diff, se = se_cluster, sig = p < .05) %>%
  mutate(outcome = factor(outcome, levels = c("GD", "RRT")),
         predictor = factor(predictor, levels = c("c[nmt]", "c[mono]")),
         lo = dllh - 1.96 * se, hi = dllh + 1.96 * se)

pred_colours <- c("c[nmt]" = "#0072B2", "c[mono]" = "#D55E00")
df <- df %>% mutate(fillcol = ifelse(sig, pred_colours[as.character(predictor)], "white"))
dodge <- position_dodge(width = 0.45)

p <- ggplot(df, aes(x = outcome, y = dllh, colour = predictor, group = predictor)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.9, position = dodge) +
  geom_point(aes(shape = predictor, fill = fillcol), size = 3.8, stroke = 1.2, position = dodge) +
  scale_colour_manual(values = pred_colours, labels = scales::parse_format(), name = NULL) +
  scale_shape_manual(values = c("c[nmt]" = 21, "c[mono]" = 24),
                     labels = scales::parse_format(), name = NULL) +
  scale_fill_identity() +
  labs(x = NULL, y = expression("elpd"[diff])) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black", linewidth = 0.3)) +
  guides(colour = guide_legend(order = 1, override.aes = list(
            shape = c(21, 24), fill = c("#0072B2", "#D55E00"))),
         shape = guide_legend(order = 1, override.aes = list(
            fill = c("#0072B2", "#D55E00"))))

ggsave(file.path(OUT_DIR, "rq3_gd_rrt.pdf"), p, width = 6, height = 4.2, device = "pdf")
cat("Saved rq3_gd_rrt.pdf\n")
