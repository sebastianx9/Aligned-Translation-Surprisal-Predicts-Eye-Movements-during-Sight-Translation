# Data

## Source corpus

All eye-tracking data come from the **EMMT corpus** (Bhattacharya et al., 2022):

> Bhattacharya, S., Kloudová, V., Zouhar, V., & Bojar, O. (2022).
> EMMT: A simultaneous eye-tracking, 4-electrode EEG and audio corpus
> for multi-modal reading and translation scenarios.
> arXiv:2204.02905. https://arxiv.org/abs/2204.02905

The corpus is publicly available at the link above. Place the
fixation CSV files in this directory before running the analysis.

## Word frequency

Word frequency data come from **SUBTLEX-US** (Brysbaert & New, 2009).
Download from: https://www.ugent.be/pp/experimentele-psychologie/en/research/documents/subtlexus

## Expected files

| File | Description |
|---|---|
| `fixation_durations_word.csv` | Per-word TFD per participant per stage |
| `nmt_surprisal_word.csv` | Aligned translation surprisal (output of `src/extract_nmt_surprisal.py`) |
| `monolingual_surprisal_word.csv` | GPT-2 surprisal (output of `src/extract_mono_surprisal.py`) |
| `attention_features_word.csv` | Encoder attention features (output of `src/extract_attention_features.py`) |
| `subtlex_us.csv` | SUBTLEX-US frequency norms |
