# RQ2 interaction decomposition figure: for each stage interaction, the pooled
# held-out gain and its translation/reading components (elpd_diff +-1.96 SE).
# Shows that the pooled null hides a translation-concentrated gain for c_nmt and
# a reading-concentrated gain for c_mono. From rq2_interaction_decomp.rds.
OUT <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(dplyr); library(ggplot2)})
res <- readRDS("/Users/sebastianx/Dissertation RQ2/rq2_interaction_decomp.rds")

res <- res %>% mutate(
  subset = factor(subset, levels = c("pooled","translation","reading"),
                  labels = c("Pooled","Translation","Reading")),
  interaction = factor(interaction, levels = c("condition:c_nmt","condition:c_mono"),
                       labels = c("condition%*%italic(c)[nmt]","condition%*%italic(c)[mono]")),
  sig = p < .05, lo = elpd_diff - 1.96*se, hi = elpd_diff + 1.96*se)

scol <- c("Pooled"="grey45","Translation"="#0072B2","Reading"="#D55E00")
p <- ggplot(res, aes(subset, elpd_diff, colour = subset)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.9) +
  geom_point(aes(shape = sig), size = 3.6, stroke = 1.2, fill = "white") +
  facet_wrap(~interaction, labeller = label_parsed) +
  scale_colour_manual(values = scol, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "p < .05", `FALSE` = "p >= .05"), name = NULL) +
  labs(x = NULL, y = expression("elpd"[diff]~"(held-out gain)")) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "top",
        axis.line = element_line(colour = "black", linewidth = 0.3))
ggsave(file.path(OUT, "rq2_decomp.pdf"), p, width = 6.6, height = 3.8, device = "pdf")
cat("saved rq2_decomp.pdf\n"); print(res[,c("interaction","subset","elpd_diff","p","sig")])
