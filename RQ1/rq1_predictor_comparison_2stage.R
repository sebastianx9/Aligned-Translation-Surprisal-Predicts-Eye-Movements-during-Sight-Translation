# RQ1 predictor comparison — TWO-PANEL (translate | reading), brms 10-fold elpd.
# Translate panel: published Table (tab:rq1_loo) values; reading panel: the new
# within-reading brms run (rq1_kfold_reading.rds). Within-stage only — the panels
# show where each predictor's predictive gain concentrates, not a cross-stage test.
OUT <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
RQ1 <- "/Users/sebastianx/Dissertation RQ1"
suppressMessages({library(dplyr); library(ggplot2)})

lv <- c("c[nmt]","c[mono]","H[e]","f[e]","f[eos]","f[recv]","f[cross]")
grp <- c(rep("Surprisal",2), rep("Attention",5))

# Translate — authoritative rq1_kfold_elpd.rds values (identical to tab:rq1_loo)
# c_mono row updated 2026-07-17 after the monolingual-surprisal alignment fix.
# Significance marker uses HOLM-corrected p (matching the formal inference
# standard of tab:rq1_loo and the text) — raw p made c_mono translation filled
# after the fix (raw .042, Holm .25), visually contradicting the Holm claims.
tr <- tibble::tibble(
  label = lv, group = grp,
  elpd = c(14.48, 4.80, -2.06, 1.12, -1.65, -1.25, 0.06),
  se   = c( 5.65, 2.77,  0.62, 2.22,  0.68,  1.12, 1.09),
  p_holm = c(.038, .25, 1, 1, 1, 1, 1),
  stage = "Translation")

# Reading — from the new within-reading brms kfold
rd_raw <- readRDS(file.path(RQ1, "rq1_kfold_reading.rds"))$res
rd <- tibble::tibble(
  label = lv, group = grp,
  elpd = rd_raw$elpd_diff[match(c("c_nmt","c_mono","H_e","f_e","f_eos","f_recv","f_cross"), rd_raw$predictor)],
  se   = rd_raw$se_cluster[match(c("c_nmt","c_mono","H_e","f_e","f_eos","f_recv","f_cross"), rd_raw$predictor)],
  p_holm = rd_raw$p_holm[match(c("c_nmt","c_mono","H_e","f_e","f_eos","f_recv","f_cross"), rd_raw$predictor)],
  stage = "Reading aloud")

df <- bind_rows(tr, rd) %>%
  mutate(label = factor(label, levels = lv),
         stage = factor(stage, levels = c("Translation","Reading aloud")),
         significant = p_holm < .05, lo = elpd - 1.96*se, hi = elpd + 1.96*se)

cols <- c(Surprisal="#0072B2", Attention="#D55E00")
p <- ggplot(df, aes(label, elpd, colour=group)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0, linewidth=0.7) +
  geom_point(aes(shape=significant), size=3, fill="white", stroke=1.1) +
  facet_wrap(~stage, ncol=1) +
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=21),
                     labels=c(`TRUE`="Holm p < .05",`FALSE`="Holm p >= .05"), name=NULL) +
  scale_colour_manual(values=cols, name=NULL) +
  scale_x_discrete(labels=function(x) parse(text=x)) +
  labs(x=NULL, y=expression("elpd"[diff]*"  (held-out gain over controls, within stage)")) +
  theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        legend.position="bottom", strip.text=element_text(face="bold"),
        axis.line.y=element_line(colour="black", linewidth=0.3)) +
  guides(colour=guide_legend(order=1), shape=guide_legend(order=2))
ggsave(file.path(OUT,"rq1_predictor_comparison_2stage.pdf"), p, width=7.0, height=6.2, device="pdf")
cat("saved rq1_predictor_comparison_2stage.pdf\n")
print(df %>% select(stage,label,elpd,se,p_holm,significant))
