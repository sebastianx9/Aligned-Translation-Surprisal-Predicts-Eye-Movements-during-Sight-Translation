# RQ1 predictor comparison figure — brms 10-fold elpd_diff (sentence-clustered SE).
OUT <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
suppressMessages({library(dplyr); library(ggplot2)})
df <- tibble::tribble(
  ~label,      ~group,      ~elpd,   ~se,    ~p,
  "c[nmt]",    "Surprisal",  14.90,  5.12,  .003,
  "c[mono]",   "Surprisal",   3.24,  3.42,  .179,
  "H[e]",      "Attention",  -1.80,  0.82,  .995,
  "f[e]",      "Attention",  -2.10,  1.05,  .982,
  "f[eos]",    "Attention",   2.85,  2.60,  .142,
  "f[recv]",   "Attention",  -2.55,  1.20,  .988,
  "f[cross]",  "Attention",  -1.95,  0.95,  .991
) %>% mutate(label=factor(label, levels=label), significant=p<.05,
             lo=elpd-1.96*se, hi=elpd+1.96*se)
cols <- c(Surprisal="#0072B2", Attention="#D55E00")
p <- ggplot(df, aes(label, elpd, colour=group)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0, linewidth=0.7) +
  geom_point(aes(shape=significant), size=3, fill="white", stroke=1.1) +
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=21),
                     labels=c(`TRUE`="p < .05",`FALSE`="p >= .05"), name=NULL) +
  scale_colour_manual(values=cols, name=NULL) +
  scale_x_discrete(labels=function(x) parse(text=x)) +
  labs(x=NULL, y=expression("elpd"[diff])) +
  theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        legend.position="bottom", axis.line.y=element_line(colour="black", linewidth=0.3)) +
  guides(colour=guide_legend(order=1), shape=guide_legend(order=2))
ggsave(file.path(OUT,"rq1_predictor_comparison.pdf"), p, width=6.5, height=4.2, device="pdf")
cat("saved rq1_predictor_comparison.pdf (brms elpd)\n")
