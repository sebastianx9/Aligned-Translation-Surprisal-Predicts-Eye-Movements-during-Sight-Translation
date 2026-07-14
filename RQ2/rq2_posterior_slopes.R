# RQ2 stage-specific slopes — posterior violins + median & 95% CrI (double dissociation).
suppressMessages({library(brms); library(ggplot2); library(dplyr); library(tidyr)})
OUT <- "/Users/sebastianx/Dissertation Writeup/MSc_and_BEng_Dissertation_Template_the_University_of_Manchester_EEE/images"
m <- readRDS("/Users/sebastianx/Dissertation_Data/brm_cache/rq2_joint_maximal.rds")
d <- as_draws_df(m)
post <- tibble(
  read_nmt=d$b_c_nmt, tran_nmt=d$b_c_nmt+d[["b_condition:c_nmt"]],
  read_mono=d$b_c_mono, tran_mono=d$b_c_mono+d[["b_condition:c_mono"]]
) %>% pivot_longer(everything(), names_to="k", values_to="v") %>%
  mutate(predictor=ifelse(grepl("nmt",k),"c[nmt]","c[mono]"),
         stage=ifelse(grepl("read",k),"reading","translation"))
p <- ggplot(post, aes(stage, v, fill=predictor)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey55") +
  geom_violin(position=position_dodge(0.8), alpha=0.45, colour=NA, width=0.8) +
  stat_summary(aes(group=predictor),
    fun=median, fun.min=function(x)quantile(x,.025), fun.max=function(x)quantile(x,.975),
    geom="pointrange", position=position_dodge(0.8), linewidth=0.6, size=0.4, colour="grey20") +
  scale_fill_manual(values=c(`c[nmt]`="#0072B2",`c[mono]`="#D55E00"),
                    labels=scales::parse_format()) +
  labs(x=NULL, y="stage-specific slope (log-ms per SD)", fill=NULL) +
  theme_minimal(base_size=12) + theme(legend.position="top",
    panel.grid.minor=element_blank())
ggsave(file.path(OUT, "rq2_posterior_slopes.pdf"), p, width=6, height=4, device="pdf")
cat("saved violin + CrI\n")
