# Aligned Translation Surprisal Predicts Source-Text Eye Movements During Sight Translation

**Student ID:** 11479116  
**Supervisor:** Dr Colin Bannard  
**Institution:** University of Manchester  
**Preregistration:** [OSF] *(link to be added upon submission)*

---

## Overview

Lim et al. (2024) showed that NMT **translation surprisal** — surprisal derived from
a neural machine translation decoder — outperforms monolingual source-language surprisal
as a predictor of human translation difficulty across 13 language pairs. Yet Lim et al.
identified a barrier to extending this advantage to source-text eye movements: decoder
surprisal is conditioned on the full source sequence simultaneously, so it cannot be
attributed to individual source words.

This project bridges that gap by introducing **aligned translation surprisal** — NMT
decoder surprisal distributed onto individual source words via **soft cross-attention
alignment** — and asks whether the same pattern Lim et al. observed on the target side
holds on the **source-text perceptual stage**: does a decoder-derived measure outperform
a pure source-language model at predicting source-text fixation durations during sight
translation?

Data come from the **EMMT corpus** (Bhattacharya et al., 2022): monocular eye-tracking
recordings of 43 Czech–English bilinguals performing reading aloud and sight translation
on 200 English sentences.

---

## Central Finding

| Predictor | β | Δllh over baseline | p (perm) |
|---|---|---|---|
| Aligned translation surprisal (`c_nmt`, soft) | **0.058** | **0.0025** | **.002** |
| Monolingual surprisal (`c_mono`, GPT-2) | — | 0.0001 | .379 |

The NMT decoder-derived measure yields a significantly larger held-out log-likelihood gain
than monolingual surprisal, mirroring Lim et al.'s target-side result on the source-text
perceptual stage. Results are based on 10-fold sentence-level cross-validation with a
paired permutation test (1,000 permutations).

---

## Repository Structure

```
.
├── src/                                    # Feature extraction (Python)
│   ├── extract_nmt_surprisal.py            # Argmax-aligned NMT surprisal (reference)
│   ├── extract_nmt_surprisal_soft.py       # Soft-aligned NMT surprisal (main method)
│   ├── plot_alignment.py                   # Cross-attention alignment heatmap figure
│   ├── extract_attention_features.py       # Encoder entropy and attention-to-context
│   └── extract_mono_surprisal.py           # GPT-2 monolingual surprisal
├── paper/
│   ├── essay.tex                           # Paper (ACL format)
│   ├── essay.bib                           # Bibliography
│   ├── acl.sty                             # ACL style file
│   ├── acl_natbib.bst                      # ACL natbib style
│   ├── alignment_map.pdf                   # Figure 1: soft alignment heatmap
│   └── delta_llh_comparison.pdf            # Figure 2: Δllh comparison
├── preregistration/
│   ├── osf_preregistration.tex
│   └── preregistration.bib
├── data/
│   └── README.md                           # Data sources and expected file layout
├── requirements.txt                        # Python dependencies
└── README.md
```

---

## Hypotheses and Models

| # | Hypothesis | Model |
|---|---|---|
| H2 | Aligned translation surprisal positively predicts TFD during sight translation | LME: `log(TFD) ~ c_nmt + controls + RE` (translate, content words) |
| H2a | `c_nmt` outperforms monolingual surprisal on source-text TFD | 10-fold CV Δllh + paired permutation test |

H2a directly mirrors Lim et al.'s (2024) target-side comparison on the source-text stage.

---

## Reproducing the Analysis

### 1. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 2. Extract NMT features

```bash
python src/extract_nmt_surprisal_soft.py   # soft-aligned surprisal (main)
python src/extract_mono_surprisal.py        # GPT-2 monolingual surprisal
python src/plot_alignment.py                # alignment heatmap figure
```

Output CSVs are written to `data/`.

### 3. Run R analysis

```r
install.packages(c("lme4", "lmerTest", "dplyr", "ggplot2"))
source("Dissertation_Data/lme_soft.R")
```

---

## Key Predictors

| Variable | Description |
|---|---|
| `c_nmt` | Soft-aligned translation surprisal: $\sum_t s_t \cdot \bar\alpha_{tw} / \sum_j \bar\alpha_{tj}$ |
| `c_mono` | Monolingual surprisal from GPT-2 |
| `c_freq` | Log word frequency from SUBTLEX-US (Lg10WF) |
| `c_wlen` | Word length (character count) |
| `c_wpos` | Normalised within-sentence position |

All continuous predictors are *z*-scored using the translate-stage content-word subset
as the reference distribution.

---

## NMT Model

All NMT features are extracted from
[`Helsinki-NLP/opus-mt-en-cs`](https://huggingface.co/Helsinki-NLP/opus-mt-en-cs)
(Tiedemann & Thottingal, 2020) via HuggingFace Transformers.

Monolingual surprisal is extracted from [`gpt2`](https://huggingface.co/gpt2).

---

## References

- Bahdanau, D., Cho, K., & Bengio, Y. (2015). Neural machine translation by jointly learning to align and translate. *ICLR 2015*.
- Bhattacharya, S., Kloudová, V., Zouhar, V., & Bojar, O. (2022). EMMT: A simultaneous eye-tracking, 4-electrode EEG and audio corpus for multi-modal reading and translation scenarios. *arXiv:2204.02905*.
- Lim, Z. W., Vylomova, E., Kemp, C., & Cohn, T. (2024). Predicting human translation difficulty with neural machine translation. *Transactions of the Association for Computational Linguistics*, 12, 1479–1496.
- Lijewska, A., Chmiel, A., & Inhoff, A. W. (2022). Stages of sight translation: Evidence from eye movements. *Applied Psycholinguistics*, 43(4), 997–1018.
- Tiedemann, J., & Thottingal, S. (2020). OPUS-MT — Building open translation services for the World. *EAMT 2020*.
- Brysbaert, M., & New, B. (2009). Moving beyond Kučera and Francis. *Behavior Research Methods*, 41(4), 977–990.
