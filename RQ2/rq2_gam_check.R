# ── RQ2: robustness check for the linear specification ────────────────────
#
# The formal RQ2 models (rq2_brm_joint.R, and the partial-residual figure
# in rq2_interaction_scatter.R) specify c_nmt and c_mono as LINEAR
# predictors of log(TFD), following the linear (mixed-effects) approach in
# Lim et al. (2024), whose delta-llh framework this thesis adopts. This
# script checks whether that linear assumption is reasonable, by comparing
# a linear fit against a GAM with a smooth (spline) term for the predictor
# of interest, on the translation-stage data (where c_nmt's effect is
# reliable).
#
# Result (see console output): a small but statistically detectable
# deviation from linearity for c_nmt in the translation stage (edf=2.18,
# likelihood-ratio test of linear-vs-GAM p=.007), unaffected by whether
# c_mono is included as a covariate, but negligible in practical size
# (AIC essentially unchanged) and uncertain at high surprisal values where
# data are sparse. Reported as a limitation in main.tex.
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR <- "/Users/sebastianx/Dissertation_Data"
OUT_DIR  <- "/Users/sebastianx/Dissertation RQ2"
suppressMessages({library(dplyr); library(ggplot2); library(mgcv)})

fix  <- read.csv(file.path(DATA_DIR, "fixation_durations_word.csv"),    stringsAsFactors=FALSE)
nmt  <- read.csv(file.path(DATA_DIR, "nmt_surprisal_soft_word.csv"),    stringsAsFactors=FALSE)
mono <- read.csv(file.path(DATA_DIR, "monolingual_surprisal_word.csv"), stringsAsFactors=FALSE)
freq <- read.table(file.path(DATA_DIR, "subtlex_us.csv"), sep="\t",
                   header=TRUE, stringsAsFactors=FALSE, quote="") %>%
  select(word=Word, log10_freq=Lg10WF) %>% mutate(word=tolower(trimws(word)))

sent_lengths <- nmt %>% group_by(sentence_id) %>%
  summarise(sent_len=max(word_index)+1, .groups="drop")

predictors <- nmt %>% left_join(sent_lengths, by="sentence_id") %>%
  mutate(word_length=nchar(word), word_position=word_index/(sent_len-1),
         word_lower=tolower(trimws(word))) %>%
  left_join(freq, by=c("word_lower"="word")) %>%
  left_join(mono %>% select(sentence_id, word_index, mono_surprisal=surprisal_sum),
            by=c("sentence_id","word_index")) %>%
  select(sentence_id, word_index, word_length, word_position,
         nmt_surprisal=surprisal_soft, mono_surprisal, log10_freq)

df <- fix %>% filter(stage %in% c("translate","read")) %>%
  left_join(predictors, by=c("sentence_id","word_index")) %>%
  mutate(log_tfd=log(total_fixation_duration_ms), ambiguity=factor(ambiguity),
         condition=factor(ifelse(stage=="translate","Translation","Reading"),
                           levels=c("Reading","Translation"))) %>%
  filter(!is.na(nmt_surprisal), !is.na(mono_surprisal), !is.na(log10_freq)) %>%
  anti_join(data.frame(sentence_id="S003", word_index=3L),
            by=c("sentence_id","word_index"))

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
df <- df %>% mutate(c_nmt=z(nmt_surprisal), c_mono=z(mono_surprisal),
                     c_wlen=z(word_length), c_wpos=z(word_position), c_freq=z(log10_freq))

# ── 1. Smooth-term significance per predictor x condition (univariate) ────
check_nonlinearity <- function(pred_var) {
  for (cond in c("Reading", "Translation")) {
    sub <- df %>% filter(condition == cond)
    f <- as.formula(paste0("log_tfd ~ s(", pred_var, ", k=6) + c_wlen + c_wpos + c_freq + ambiguity"))
    g <- gam(f, data = sub, method = "REML")
    s <- summary(g)
    cat(sprintf("%s | %s: edf=%.2f (1=linear), smooth p=%.4f\n",
                pred_var, cond, s$s.table[1, "edf"], s$s.table[1, "p-value"]))
  }
}
cat("=== Univariate smooth-term check ===\n")
check_nonlinearity("c_nmt")
check_nonlinearity("c_mono")

# ── 2. Does controlling for c_mono change c_nmt's curvature in Translation? ──
sub <- df %>% filter(condition == "Translation")
cat("\n=== c_nmt curvature in Translation, with/without c_mono control ===\n")
g1 <- gam(log_tfd ~ s(c_nmt, k=6) + c_wlen + c_wpos + c_freq + ambiguity, data=sub, method="REML")
g2 <- gam(log_tfd ~ s(c_nmt, k=6) + c_mono + c_wlen + c_wpos + c_freq + ambiguity, data=sub, method="REML")
cat(sprintf("Without c_mono: edf=%.2f, p=%.4f\n", summary(g1)$s.table[1,"edf"], summary(g1)$s.table[1,"p-value"]))
cat(sprintf("With c_mono:    edf=%.2f, p=%.4f\n", summary(g2)$s.table[1,"edf"], summary(g2)$s.table[1,"p-value"]))

# ── 3. Formal test: linear vs GAM nested model comparison (the correct test
# for "is the deviation from linearity itself significant", as opposed to
# the smooth-term p-value above, which tests overall smooth significance) ──
cat("\n=== Linear vs GAM nested comparison (both with c_mono controlled) ===\n")
g_lin <- gam(log_tfd ~ c_nmt + c_mono + c_wlen + c_wpos + c_freq + ambiguity, data=sub, method="ML")
g_gam <- gam(log_tfd ~ s(c_nmt, k=6) + c_mono + c_wlen + c_wpos + c_freq + ambiguity, data=sub, method="ML")
cat("Linear AIC:", AIC(g_lin), " | GAM AIC:", AIC(g_gam), "\n")
print(anova(g_lin, g_gam, test="Chisq"))

# ── 4. Plot: GAM smooth vs linear fit for c_nmt in Translation ────────────
g_plot <- gam(log_tfd ~ s(c_nmt, k=6) + c_wlen + c_wpos + c_freq + ambiguity, data=sub, method="REML")
newx <- data.frame(c_nmt = seq(min(sub$c_nmt), max(sub$c_nmt), length.out=100),
                    c_wlen=0, c_wpos=0, c_freq=0, ambiguity=factor("A", levels=levels(sub$ambiguity)))
pred <- predict(g_plot, newx, se.fit=TRUE, terms="s(c_nmt)")
newx$fit <- pred$fit; newx$se <- pred$se.fit
lm_fit <- lm(log_tfd ~ c_nmt + c_wlen + c_wpos + c_freq + ambiguity, data=sub)
newx$lm_pred <- predict(lm_fit, newx)

p <- ggplot() +
  geom_ribbon(data=newx, aes(x=c_nmt, ymin=fit-1.96*se, ymax=fit+1.96*se), alpha=0.2, fill="#0072B2") +
  geom_line(data=newx, aes(x=c_nmt, y=fit), colour="#0072B2", linewidth=1) +
  geom_line(data=newx, aes(x=c_nmt, y=lm_pred - mean(lm_pred) + mean(fit)), colour="black", linetype="dashed") +
  labs(title="c_nmt effect in Translation: GAM (blue/shaded) vs linear (dashed)",
       x="c_nmt (SD units)", y="Smooth term contribution to log(TFD)") +
  theme_minimal(base_size=12)
ggsave(file.path(OUT_DIR, "gam_check_cnmt_translation.pdf"), p, width=6, height=4.5, device="pdf")
cat("\nSaved gam_check_cnmt_translation.pdf\n")
