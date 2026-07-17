# RQ2 predictive line: symmetric nested-ladder elpd gains (forest plot).
suppressMessages({library(ggplot2); library(dplyr)})
OUT <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
d <- tibble::tribble(
  ~label,                    ~predictor, ~kind,  ~elpd,  ~se,
  "c_nmt: stage interaction",  "c[nmt]",  "int",   3.06,  2.74,
  "c_nmt: main effect",        "c[nmt]",  "main", 12.76,  4.73,
  "c_mono: stage interaction", "c[mono]", "int",   2.14,  2.58,
  "c_mono: main effect",       "c[mono]", "main", 10.66,  4.83
) %>% mutate(label=factor(label, levels=label),
             lo=elpd-1.96*se, hi=elpd+1.96*se)
cols <- c(`c[nmt]`="#0072B2", `c[mono]`="#D55E00")
p <- ggplot(d, aes(elpd, label, colour=predictor)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey55") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0, linewidth=0.8) +
  geom_point(aes(shape=kind), size=3.2, fill="white", stroke=1.1) +
  scale_shape_manual(values=c(main=16, int=21), guide="none") +
  scale_colour_manual(values=cols, labels=scales::parse_format(), name=NULL) +
  labs(x=expression("held-out "*elpd[diff]*" over the alternative predictor"), y=NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position="top", panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(), axis.line.x=element_line(colour="black", linewidth=0.3))
ggsave(file.path(OUT,"rq2_predictive_ladder.pdf"), p, width=6.6, height=3.2, device="pdf")
cat("saved\n")
