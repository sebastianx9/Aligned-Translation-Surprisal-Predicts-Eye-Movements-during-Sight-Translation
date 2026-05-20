# Aligned Translation Surprisal Predicts Source-Text Eye Movements During Sight Translation

Code for the paper: **Aligned Translation Surprisal Predicts Source-Text Eye Movements During Sight Translation** (2025).

## Overview

We introduce *aligned translation surprisal*, NMT decoder surprisal soft-aligned to source words via cross-attention weights, and test whether it predicts source-text fixation durations during sight translation. Eye-tracking data come from the EMMT corpus (Bhattacharya et al., 2022).

## Requirements

**Python** (feature extraction):
```bash
pip install -r requirements.txt
```

**R** (analysis): The notebook uses an R kernel via [IRkernel](https://irkernel.github.io/).
```r
install.packages(c("lme4", "lmerTest", "dplyr", "ggplot2"))
```

## Data

Eye-tracking data are from the [EMMT corpus](https://github.com/rgu-iit-bt/emmt) (Bhattacharya et al., 2022) and are not redistributed here. Word frequency norms are from SUBTLEX-US (Brysbaert & New, 2009).

Run the extraction scripts in `src/` against your local copy of the EMMT corpus, then set `DATA_DIR` in `result-analysis.ipynb` to the folder containing the output CSV files.

## Reproducing the Analysis

### Prerequisites

1. Clone this repository:
   ```bash
   git clone https://github.com/<your-username>/nmt-surprisal-sight-translation.git
   cd nmt-surprisal-sight-translation
   ```

2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Download the [EMMT corpus](https://github.com/rgu-iit-bt/emmt) and note the path to its `preprocessed-data/gaze/` directory.

4. Download [SUBTLEX-US](https://www.ugent.be/pp/experimentele-psychologie/en/research/subtlex) and save as `subtlex_us.csv`.

### Step 1 — Extract features

```bash
# Word-level total fixation duration from raw EMMT gaze files
python src/extract_fixation_duration.py \
  --read_dir      /path/to/emmt/preprocessed-data/gaze/Read \
  --translate_dir /path/to/emmt/preprocessed-data/gaze/Translate \
  --sentences     /path/to/emmt/probes/Sentences.csv \
  --output        fixation_durations_word.csv

# Aligned translation surprisal (c_nmt)
# Uses Helsinki-NLP/opus-mt-en-cs (Tiedemann & Thottingal, 2020);
# model is downloaded automatically on first run.
python src/extract_nmt_surprisal_soft.py \
  --sentences /path/to/emmt/probes/Sentences.csv \
  --output    nmt_surprisal_soft_word.csv

# Monolingual surprisal (c_mono), uses GPT-2
python src/extract_monolingual_surprisal.py \
  --sentences /path/to/emmt/probes/Sentences.csv \
  --output    monolingual_surprisal_word.csv

# Optional: cross-attention alignment heatmap
python src/plot_alignment.py \
  --sentence "A man is blowing into a plastic ball." \
  --output   alignment_map.pdf
```

### Step 2 — Run the analysis notebook

Install R dependencies (once):
```r
install.packages(c("lme4", "lmerTest", "dplyr", "ggplot2"))
```

Open the notebook, set `DATA_DIR` to the folder containing the CSV files from Step 1, and run all cells:
```bash
jupyter notebook result-analysis.ipynb
```

Figures are saved to the `figures/` directory.

## Repository Structure

| Path | Description |
|------|-------------|
| `src/extract_fixation_duration.py` | Raw EMMT gaze files → word-level TFD |
| `src/extract_nmt_surprisal_soft.py` | Aligned translation surprisal via soft cross-attention alignment |
| `src/extract_monolingual_surprisal.py` | GPT-2 monolingual surprisal |
| `src/plot_alignment.py` | Cross-attention alignment visualisation |
| `result-analysis.ipynb` | Mixed-effects models and cross-validation (R kernel) |

## References

- Bhattacharya et al. (2022). EMMT corpus. *LREC*.
- Lim et al. (2024). Predicting Translation Difficulty with NMT. *ACL*.
- Tiedemann, J. and Thottingal, S. (2020). OPUS-MT. *EAMT*.
- Brysbaert & New (2009). SUBTLEX-US. *Behavior Research Methods*.
