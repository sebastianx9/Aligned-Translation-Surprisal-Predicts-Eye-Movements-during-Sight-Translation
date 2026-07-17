# ── RQ2: stage-specific slope crossover plot ──────────────────────────────
# Reconstructed from the reported joint-model posterior slopes (main.tex).
# Blue = c_nmt, green = c_mono; filled = 95% CI excludes zero, open = spans zero.

OUT_DIR <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(dplyr); library(ggplot2)})

df <- tibble::tribble(
  ~predictor, ~stage,        ~slope, ~lo,     ~hi,    ~sig,
  "c[nmt]",   "Reading",      0.017, -0.010,  0.043,  FALSE,
  "c[nmt]",   "Translation",  0.056,  0.029,  0.083,  TRUE,
  "c[mono]",  "Reading",      0.045,  0.020,  0.070,  TRUE,
  "c[mono]",  "Translation",  0.008, -0.018,  0.034,  FALSE
) %>%
  mutate(stage = factor(stage, levels = c("Reading", "Translation")),
         predictor = factor(predictor, levels = c("c[nmt]", "c[mono]")))

pred_colours <- c("c[nmt]" = "#0072B2", "c[mono]" = "#D55E00")

p <- ggplot(df, aes(x = stage, y = slope, colour = predictor, group = predictor)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.06, linewidth = 0.9) +
  geom_point(aes(shape = sig), size = 4, stroke = 1.2, fill = "white") +
  scale_colour_manual(values = pred_colours, labels = scales::parse_format(), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "CI excludes zero", `FALSE` = "CI spans zero"),
                     name = NULL) +
  labs(x = NULL, y = "Stage-specific slope (log ms per SD)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black", linewidth = 0.3)) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)),
         shape = guide_legend(order = 2))

ggsave(file.path(OUT_DIR, "rq2_crossover.pdf"), p, width = 6.5, height = 4.2, device = "pdf")
cat("Saved rq2_crossover.pdf\n")
