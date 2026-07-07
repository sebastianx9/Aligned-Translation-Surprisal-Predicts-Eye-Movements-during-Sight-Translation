# Aligned Translation Surprisal Predicts Source-Text Eye Movements During Sight Translation

Code for the MSc dissertation extending the course paper of the same title (LELA70502) into a
three-part investigation (RQ1–RQ3) of whether NMT-derived surprisal predicts source-text fixation
durations during sight translation, and whether that effect is translation-specific.

Student ID: 11479116

## Overview

*Aligned translation surprisal* ($c_\text{nmt}$) is NMT decoder surprisal soft-aligned to source
words via cross-attention weights. Eye-tracking data come from the EMMT corpus
(Bhattacharya et al., 2022).

- **RQ1** — Does $c_\text{nmt}$ predict source-text fixation durations, and how does it compare to
  monolingual surprisal ($c_\text{mono}$) and NMT encoder attention features?
- **RQ2** — Does $c_\text{nmt}$ trace translation-specific difficulty beyond the general reading
  difficulty captured by $c_\text{mono}$?
- **RQ3** — Which aspect of source-text processing (first-pass gaze vs. re-reading) does
  $c_\text{nmt}$'s contribution reflect?

## Requirements

**Python** (feature extraction):
```bash
pip install -r requirements.txt
```

**R** (analysis):
```r
install.packages(c("lme4", "lmerTest", "dplyr", "ggplot2", "brms"))
```

## Data

Eye-tracking data are from the [EMMT corpus](https://github.com/rgu-iit-bt/emmt)
(Bhattacharya et al., 2022) and are not redistributed here. Word frequency norms are from
SUBTLEX-US (Brysbaert & New, 2009).

Run the extraction scripts in `data-extraction/` against your local copy of the EMMT corpus to
produce the word-level CSV files (`fixation_durations_word.csv`, `nmt_surprisal_soft_word.csv`,
`monolingual_surprisal_word.csv`) that every script in `RQ1/`, `RQ2/`, and `RQ3/` reads via a
`DATA_DIR` variable at the top of each file.

## Repository structure

```
├── data-extraction/                     # raw EMMT -> word-level feature CSVs
│   ├── extract_fixation_duration.py     # gaze files -> word-level TFD
│   ├── extract_nmt_surprisal_soft.py    # aligned translation surprisal (c_nmt)
│   ├── extract_monolingual_surprisal.py # GPT-2 monolingual surprisal (c_mono)
│   ├── plot_alignment.py                # cross-attention alignment heatmap
│   └── result-analysis.ipynb            # original course-paper (RQ1-only) analysis notebook
│
├── RQ1/                                  # predictive-power comparison (200-fold LOO-CV)
│   ├── rq1_loo_data_inspect_copy_failure.R   # diagnoses the "stoplight" copy-failure artefact
│   ├── copy_failures_full_scan.csv           # output of the above: all 26 flagged word positions
│   ├── rq1_variance_decomp.R                 # variance decomposition behind the Limitations
│   │                                          #   sentence-level precision discussion
│   ├── rq1_maximal_re_check.R                # robustness check: maximal random effects +
│   └── rq1_maximal_re_output.txt             #   Holm-corrected multiple comparisons (200x7 fits)
│
├── RQ2/                                  # translation-specificity test
│   ├── rq2_norm_loo.R / rq2_norm_output.txt      # RQ1 Table: normalised attention features,
│   │                                              #   200-fold LOO-CV vs baseline / vs baseline+c_nmt
│   ├── rq2_read_stage.R / rq2_read_output.txt    # reading-stage 200-fold LOO-CV (cross-stage table)
│   ├── rq2_brm_joint.R / rq2_joint_output.txt    # M_joint: the Bayesian mixed-effects model behind
│   │                                              #   all reported condition x predictor coefficients
│   ├── rq2_stage_slopes_check.R                  # in-sample lmer diagnostic: per-stage, solo vs.
│   │                                              #   joint coefficients (non-redundancy check)
│   ├── rq2_interaction_scatter.R                 # Figure: partial-residual binned scatterplot
│   ├── rq2_loo_results.csv / rq2_6features_results.csv   # supporting result tables
│
├── RQ3/                                  # processing-phase decomposition (GD vs RRT)
│   └── rq3_joint_coef_check.R            # joint-model coefficient check: does each predictor's
│                                          #   GD/RRT effect survive controlling for the other
│                                          #   (mirrors RQ2's non-redundancy logic)
│
└── figures/                              # the three figures actually included in the dissertation
    ├── rq2_crossover.pdf                 # stage-specific slopes (posterior mean + 95% CI)
    ├── rq2_interaction_scatter.pdf       # partial-residual binned scatterplot (see RQ2 script)
    └── rq3_gd_rrt.pdf                    # GD/RRT predictive gains
```

## A known gap

The scripts that produced the exact headline $\Delta\text{llh}$ figures for RQ1's primary
$c_\text{nmt}$ result and RQ3's GD/RRT decomposition could not be located among the working
files at the time this repository was assembled — only reference copies of their output values,
hardcoded into other diagnostic scripts, survive. `rq1_maximal_re_check.R` and
`rq3_joint_coef_check.R` reproduce the same data pipeline and give matching or near-identical
values as a robustness check, but are not verbatim the original computation. If you are trying to
exactly reproduce Table 1 (RQ1) or the GD/RRT $\Delta\text{llh}$ table (RQ3) and hit a discrepancy,
this is why.

## References

- Bhattacharya et al. (2022). EMMT corpus. *Scientific Data*.
- Lim et al. (2024). Predicting Human Translation Difficulty with Neural Machine Translation. *TACL*.
- Wilcox et al. (2023). Testing the Predictions of Surprisal Theory in 11 Languages. *TACL*.
- Tiedemann, J. and Thottingal, S. (2020). OPUS-MT. *EAMT*.
- Brysbaert & New (2009). SUBTLEX-US. *Behavior Research Methods*.
- Bürkner, P.-C. (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. *JSS*.
